#!/usr/bin/env python3
"""
Pipeline alineado al flujo de alphafold_multimer_validation_deep_dive.md:

1) PDB complejo (target + binder backbone)
2) ProteinMPNN diseña SOLO la cadena binder
3) Construcción de FASTA por diseño: seq_target_nativa : seq_binder_mpnn
4) ColabFold AlphaFold-Multimer v3 (num_recycle=3) con initial_guess=PDB
5) Métricas por diseño (rank_001):
   - ipTM
   - iPAE (promedio intercadena, Å)
   - ipSAE (opcional de apoyo)
   - RMSD_CA_all con alineación por target CA

Nota importante:
- El deep_dive describe la implementación de ColabDesign binder protocol con
  templates asimétricos (target visible, binder enmascarado). Este script usa
  ColabFold CLI + initial_guess para aproximar ese comportamiento.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}


def run(cmd: List[str], cwd: str | None = None) -> None:
    print("\n[RUN]", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=cwd)


def parse_pdb_chain_sequences(pdb_path: Path) -> OrderedDict[str, str]:
    chains: OrderedDict[str, List[str]] = OrderedDict()
    seen = set()
    with pdb_path.open("r") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            resname = line[17:20].strip().upper()
            if resname not in AA3_TO_1:
                continue
            altloc = line[16].strip()
            if altloc not in ("", "A"):
                continue
            chain = line[21].strip() or "_"
            resseq = line[22:26].strip()
            icode = line[26].strip()
            resid = (chain, resseq, icode)
            if resid in seen:
                continue
            seen.add(resid)
            chains.setdefault(chain, [])
            chains[chain].append(AA3_TO_1[resname])
    return OrderedDict((c, "".join(seq)) for c, seq in chains.items())


def run_proteinmpnn(
    python_bin: str,
    mpnn_script: Path,
    pdb: Path,
    binder_chain: str,
    out_dir: Path,
    num_seq: int,
    temp: str,
    seed: int,
    batch_size: int,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        python_bin,
        str(mpnn_script),
        "--pdb_path",
        str(pdb),
        "--pdb_path_chains",
        binder_chain,
        "--num_seq_per_target",
        str(num_seq),
        "--sampling_temp",
        temp,
        "--seed",
        str(seed),
        "--batch_size",
        str(batch_size),
        "--out_folder",
        str(out_dir),
    ]
    run(cmd)
    seq_file = out_dir / "seqs" / f"{pdb.stem}.fa"
    if not seq_file.exists():
        raise FileNotFoundError(f"No se encontró salida MPNN: {seq_file}")
    return seq_file


def parse_mpnn_samples(seq_file: Path) -> List[Dict[str, float | str | int]]:
    """
    Devuelve muestras MPNN en el orden del archivo.
    Cada muestra tiene sample_id, score, global_score, seq_recovery y binder_seq.
    """
    samples: List[Dict[str, float | str | int]] = []
    current_meta: Dict[str, float | str | int] | None = None

    re_sample = re.compile(
        r"sample=(\d+),\s*score=([0-9eE+.\-]+),\s*global_score=([0-9eE+.\-]+),\s*seq_recovery=([0-9eE+.\-]+)"
    )

    with seq_file.open("r") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                m = re_sample.search(line)
                if m:
                    current_meta = {
                        "sample_id": int(m.group(1)),
                        "mpnn_score": float(m.group(2)),
                        "mpnn_global_score": float(m.group(3)),
                        "mpnn_seq_recovery": float(m.group(4)),
                    }
                else:
                    current_meta = None
                continue

            # Línea de secuencia
            if current_meta is not None:
                meta = dict(current_meta)
                meta["binder_seq"] = line
                samples.append(meta)
                current_meta = None

    if not samples:
        raise ValueError(f"No se encontraron muestras MPNN en {seq_file}")
    return samples


def build_multimer_fastas(
    target_seq: str,
    mpnn_samples: List[Dict[str, float | str | int]],
    out_dir: Path,
) -> List[Tuple[Path, Dict[str, float | str | int]]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    out: List[Tuple[Path, Dict[str, float | str | int]]] = []

    for i, sample in enumerate(mpnn_samples, start=1):
        binder_seq = str(sample["binder_seq"])
        design_name = f"design_{i:03d}"
        fasta_path = out_dir / f"{design_name}.fasta"
        with fasta_path.open("w") as f:
            f.write(f">{design_name}\n{target_seq}:{binder_seq}\n")
        out.append((fasta_path, sample))

    return out


def run_colabfold_for_fastas(
    colabfold_python: str,
    fasta_and_meta: List[Tuple[Path, Dict[str, float | str | int]]],
    out_root: Path,
    model_type: str,
    num_recycle: int,
    num_models: int,
    num_seeds: int,
    random_seed: int,
    save_all: bool,
    overwrite_existing: bool,
    initial_guess_pdb: Path | None,
) -> List[Tuple[Path, Dict[str, float | str | int]]]:
    runs: List[Tuple[Path, Dict[str, float | str | int]]] = []
    for fasta_path, meta in fasta_and_meta:
        design_name = fasta_path.stem
        run_dir = out_root / f"af_{design_name}"
        run_dir.mkdir(parents=True, exist_ok=True)
        cmd = [
            colabfold_python,
            "-m",
            "colabfold.batch",
            "--model-type",
            model_type,
            "--num-recycle",
            str(num_recycle),
            "--num-models",
            str(num_models),
            "--num-seeds",
            str(num_seeds),
            "--random-seed",
            str(random_seed),
        ]
        if initial_guess_pdb is not None:
            cmd.extend(["--initial-guess", str(initial_guess_pdb)])
        if save_all:
            cmd.append("--save-all")
        if overwrite_existing:
            cmd.append("--overwrite-existing-results")
        cmd.extend([str(fasta_path), str(run_dir)])
        run(cmd)
        runs.append((run_dir, meta))
    return runs


def d0_from_n(n: int) -> float:
    if n >= 27:
        d0 = 1.24 * ((n - 15) ** (1 / 3)) - 1.8
    else:
        d0 = 1.0
    return max(1.0, d0)


def asym_ipsae(pae: List[List[float]], aligned_idx: range, partner_idx: range, pae_cutoff: float) -> float:
    best = 0.0
    for i in aligned_idx:
        js = [j for j in partner_idx if pae[i][j] < pae_cutoff]
        if not js:
            continue
        d0 = d0_from_n(len(js))
        vals = [1.0 / (1.0 + (pae[i][j] / d0) ** 2) for j in js]
        s = sum(vals) / len(vals)
        if s > best:
            best = s
    return best


def compute_metrics_from_scores_json(scores_json: Path, len_a: int, len_b: int, pae_cutoff: float) -> Dict[str, float]:
    d = json.loads(scores_json.read_text())
    pae = d["pae"]
    total_len = len_a + len_b
    a_idx = range(0, len_a)
    b_idx = range(len_a, total_len)

    ab_vals = [pae[i][j] for i in a_idx for j in b_idx]
    ba_vals = [pae[i][j] for i in b_idx for j in a_idx]
    i_pae = (sum(ab_vals) + sum(ba_vals)) / (len(ab_vals) + len(ba_vals))

    ipsae_ab = asym_ipsae(pae, a_idx, b_idx, pae_cutoff)
    ipsae_ba = asym_ipsae(pae, b_idx, a_idx, pae_cutoff)
    ipsae = max(ipsae_ab, ipsae_ba)

    plddt = d.get("plddt", [])
    plddt_mean = float(sum(plddt) / len(plddt)) if plddt else math.nan

    return {
        "iptm": float(d.get("iptm", math.nan)),
        "ptm": float(d.get("ptm", math.nan)),
        "plddt_mean": plddt_mean,
        "i_pae": float(i_pae),
        "ipsae": float(ipsae),
        "ipsae_ab": float(ipsae_ab),
        "ipsae_ba": float(ipsae_ba),
    }


def find_rank1_scores_json(af_out_dir: Path) -> Path | None:
    rank1 = sorted(af_out_dir.glob("*scores_rank_001*.json"))
    if rank1:
        return rank1[0]
    any_scores = sorted(af_out_dir.glob("*scores_rank_*.json"))
    if any_scores:
        print(f"[WARN] No encontré rank_001; usando {any_scores[0].name}")
        return any_scores[0]
    return None


def find_rank1_pdb(af_out_dir: Path) -> Path | None:
    rank1 = sorted(af_out_dir.glob("*unrelaxed_rank_001*.pdb"))
    if rank1:
        return rank1[0]
    any_pdb = sorted(af_out_dir.glob("*unrelaxed_rank_*.pdb"))
    if any_pdb:
        print(f"[WARN] No encontré rank_001 PDB; usando {any_pdb[0].name}")
        return any_pdb[0]
    return None


def parse_ca_coords_by_chain(pdb_path: Path) -> Dict[str, np.ndarray]:
    by_chain: Dict[str, List[np.ndarray]] = {}
    seen = set()
    with pdb_path.open("r") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            atom = line[12:16].strip()
            if atom != "CA":
                continue
            altloc = line[16].strip()
            if altloc not in ("", "A"):
                continue
            chain = line[21].strip() or "_"
            resseq = line[22:26].strip()
            icode = line[26].strip()
            resid = (chain, resseq, icode)
            if resid in seen:
                continue
            seen.add(resid)
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            by_chain.setdefault(chain, []).append(np.array([x, y, z], dtype=float))
    return {c: np.vstack(v) for c, v in by_chain.items() if v}


def kabsch_align(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Devuelve R, t tal que Q_aligned = Q @ R + t se alinea con P.
    """
    Pc = P.mean(axis=0)
    Qc = Q.mean(axis=0)
    P0 = P - Pc
    Q0 = Q - Qc
    H = Q0.T @ P0
    U, _, Vt = np.linalg.svd(H)
    R = U @ Vt
    if np.linalg.det(R) < 0:
        U[:, -1] *= -1
        R = U @ Vt
    t = Pc - Qc @ R
    return R, t


def rmsd_all_aligned_on_target(
    true_pdb: Path,
    pred_pdb: Path,
    target_chain: str,
    binder_chain: str,
) -> float:
    true_ca = parse_ca_coords_by_chain(true_pdb)
    pred_ca = parse_ca_coords_by_chain(pred_pdb)
    if target_chain not in true_ca or binder_chain not in true_ca:
        return math.nan
    if target_chain not in pred_ca or binder_chain not in pred_ca:
        return math.nan

    t_true = true_ca[target_chain]
    t_pred = pred_ca[target_chain]
    n_t = min(len(t_true), len(t_pred))
    if n_t < 3:
        return math.nan

    # Ajuste sobre target
    R, t = kabsch_align(t_true[:n_t], t_pred[:n_t])

    # RMSD sobre target+binder (como deep_dive describe tras alineación por target)
    b_true = true_ca[binder_chain]
    b_pred = pred_ca[binder_chain]
    n_b = min(len(b_true), len(b_pred))
    if n_b == 0:
        return math.nan

    pred_all = np.vstack([t_pred[:n_t], b_pred[:n_b]])
    true_all = np.vstack([t_true[:n_t], b_true[:n_b]])
    pred_all_aln = pred_all @ R + t
    diff = true_all - pred_all_aln
    return float(np.sqrt((diff * diff).sum() / len(true_all)))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--pdb", required=True, help="PDB del complejo target+binder")
    p.add_argument("--target-chain", default="A")
    p.add_argument("--binder-chain", default="B")
    p.add_argument("--proteinmpnn-python", default="python")
    p.add_argument("--proteinmpnn-script", required=True)
    p.add_argument("--num-mpnn-seqs", type=int, default=16)
    p.add_argument("--mpnn-temp", default="0.2")
    p.add_argument("--mpnn-seed", type=int, default=42)
    p.add_argument("--mpnn-batch-size", type=int, default=8)
    p.add_argument("--colabfold-python", required=True)
    p.add_argument("--model-type", default="alphafold2_multimer_v3")
    p.add_argument("--num-recycle", type=int, default=3)
    p.add_argument("--num-models", type=int, default=1, help="Deep dive usa model_1_multimer_v3")
    p.add_argument("--num-seeds", type=int, default=1)
    p.add_argument("--random-seed", type=int, default=0)
    p.add_argument("--pae-cutoff", type=float, default=15.0)
    p.add_argument("--no-save-all", action="store_true")
    p.add_argument("--overwrite-existing-results", action="store_true")
    p.add_argument("--no-initial-guess", action="store_true")
    p.add_argument("--max-designs", type=int, default=0, help="0=usar todas las secuencias MPNN")
    p.add_argument("--out-root", required=True)
    args = p.parse_args()

    pdb = Path(args.pdb).resolve()
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    chain_seqs = parse_pdb_chain_sequences(pdb)
    if args.target_chain not in chain_seqs or args.binder_chain not in chain_seqs:
        raise ValueError(f"Cadenas no encontradas. Disponibles: {list(chain_seqs.keys())}")

    target_seq = chain_seqs[args.target_chain]
    print(f"[INFO] Len {args.target_chain}={len(chain_seqs[args.target_chain])}, Len {args.binder_chain}={len(chain_seqs[args.binder_chain])}")

    # 1) ProteinMPNN
    mpnn_out = out_root / "mpnn_out"
    seq_file = run_proteinmpnn(
        python_bin=args.proteinmpnn_python,
        mpnn_script=Path(args.proteinmpnn_script).resolve(),
        pdb=pdb,
        binder_chain=args.binder_chain,
        out_dir=mpnn_out,
        num_seq=args.num_mpnn_seqs,
        temp=args.mpnn_temp,
        seed=args.mpnn_seed,
        batch_size=args.mpnn_batch_size,
    )
    samples = parse_mpnn_samples(seq_file)
    if args.max_designs > 0:
        samples = samples[: args.max_designs]

    # 2) FASTA target(native):binder(MPNN)
    fastas_out = out_root / "mpnn_fastas"
    fasta_and_meta = build_multimer_fastas(target_seq=target_seq, mpnn_samples=samples, out_dir=fastas_out)
    print(f"[INFO] FASTAs generados: {len(fasta_and_meta)}")

    # 3) AF-Multimer
    runs = run_colabfold_for_fastas(
        colabfold_python=args.colabfold_python,
        fasta_and_meta=fasta_and_meta,
        out_root=out_root,
        model_type=args.model_type,
        num_recycle=args.num_recycle,
        num_models=args.num_models,
        num_seeds=args.num_seeds,
        random_seed=args.random_seed,
        save_all=not args.no_save_all,
        overwrite_existing=args.overwrite_existing_results,
        initial_guess_pdb=None if args.no_initial_guess else pdb,
    )

    # 4) métricas
    len_a = len(chain_seqs[args.target_chain])
    len_b = len(chain_seqs[args.binder_chain])
    rows = []
    for run_dir, meta in runs:
        scores_json = find_rank1_scores_json(run_dir)
        if scores_json is None:
            print(f"[WARN] No encontré scores rank en {run_dir}")
            continue
        metrics = compute_metrics_from_scores_json(scores_json, len_a, len_b, args.pae_cutoff)

        pred_pdb = find_rank1_pdb(run_dir)
        rmsd_all = math.nan
        if pred_pdb is not None:
            rmsd_all = rmsd_all_aligned_on_target(
                true_pdb=pdb,
                pred_pdb=pred_pdb,
                target_chain=args.target_chain,
                binder_chain=args.binder_chain,
            )

        row = {
            "design": run_dir.name.replace("af_", ""),
            "run_dir": str(run_dir),
            "scores_json_used": str(scores_json),
            "pred_pdb_used": str(pred_pdb) if pred_pdb else "",
            "mpnn_sample_id": meta.get("sample_id", ""),
            "mpnn_score": meta.get("mpnn_score", math.nan),
            "mpnn_global_score": meta.get("mpnn_global_score", math.nan),
            "mpnn_seq_recovery": meta.get("mpnn_seq_recovery", math.nan),
            "rmsd_ca_all_align_target": rmsd_all,
            **metrics,
        }
        rows.append(row)
        print(
            f"[METRIC] {row['design']} "
            f"ipTM={row['iptm']:.3f} iPAE={row['i_pae']:.3f} "
            f"ipSAE={row['ipsae']:.3f} RMSD={row['rmsd_ca_all_align_target']:.3f}"
        )

    summary_csv = out_root / "summary_metrics.csv"
    with summary_csv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "design",
                "run_dir",
                "scores_json_used",
                "pred_pdb_used",
                "mpnn_sample_id",
                "mpnn_score",
                "mpnn_global_score",
                "mpnn_seq_recovery",
                "iptm",
                "ptm",
                "plddt_mean",
                "i_pae",
                "ipsae",
                "ipsae_ab",
                "ipsae_ba",
                "rmsd_ca_all_align_target",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n[OK] Resumen guardado en: {summary_csv}")


if __name__ == "__main__":
    main()

