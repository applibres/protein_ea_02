#!/usr/bin/env python3
"""
Descarga complejos de SKEMPI v2.0 y extrae automáticamente las posiciones
de interfaz donde se han hecho mutaciones experimentales.

Útil para validar tu GA con datos reales.

Usage:
    python skempi_downloader.py --pdb 1B27 --output ./test_case/
    python skempi_downloader.py --list-all  # Ver todos los complejos disponibles
"""

import pandas as pd
import requests
import argparse
import re
import math
import time
import difflib
from pathlib import Path

SKEMPI_CSV_URLS = [
    "https://life.bsc.es/pid/skempi2/database/download/skempi_v2.csv",
    "https://life.bsc.es/pid/skempi2/data/skempi_v2.csv",
]
PDB_URL = "https://files.rcsb.org/download/"

class SKEMPIDownloader:
    def __init__(self, cache_dir="./skempi_cache"):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        self.csv_path = self.cache_dir / "skempi_v2.csv"
        self.df = None
        self.colmap = None

    def _find_column(self, df, *names):
        for name in names:
            if name in df.columns:
                return name
        return None

    def _detect_columns(self, df):
        colmap = {
            "pdb": self._find_column(df, "PDB", "#Pdb", "Pdb"),
            "mutation": self._find_column(df, "Mutation(s) in PDB", "Mutation(s)_PDB"),
            "ddg": self._find_column(df, "DDG(kcal/mol)", "DDG", "ddG"),
            "aff_mut": self._find_column(df, "Affinity_Kd_mut(nM)", "Affinity_mut_parsed", "Affinity_mut (M)"),
            "aff_wt": self._find_column(df, "Affinity_Kd_wt(nM)", "Affinity_wt_parsed", "Affinity_wt (M)"),
        }
        missing = [k for k in ("pdb", "mutation") if colmap[k] is None]
        if missing:
            raise ValueError(f"CSV SKEMPI no reconocido, faltan columnas: {missing}")
        return colmap

    def _to_nM(self, value, source_col):
        if pd.isna(value):
            return float("nan")
        # En el export actual las columnas parsed y (M) están en molar.
        if source_col and ("parsed" in source_col or "(M)" in source_col):
            return float(value) * 1e9
        return float(value)

    def _parse_mutation(self, mut):
        mut = mut.strip()
        # Formato SKEMPI frecuente: KA27A -> WT K, cadena A, posición 27, mutante A
        m = re.match(r"^([A-Za-z])([A-Za-z])(-?\d+)([A-Za-z])$", mut)
        if m:
            return m.group(2), int(m.group(3)), m.group(1), m.group(4)
        # Formato legado simple: A123B -> cadena A, posición 123, mutante B (WT desconocido)
        m = re.match(r"^([A-Za-z])(-?\d+)([A-Za-z])$", mut)
        if m:
            return m.group(1), int(m.group(2)), "?", m.group(3)
        return None
        
    def download_csv(self):
        """Descarga el CSV de SKEMPI si no existe"""
        if self.csv_path.exists():
            print(f"CSV encontrado en {self.csv_path}")
        else:
            last_error = None
            for url in SKEMPI_CSV_URLS:
                for attempt in range(1, 4):
                    try:
                        print(f"Descargando SKEMPI v2.0 desde {url} (intento {attempt}/3)...")
                        response = requests.get(
                            url,
                            timeout=30,
                            headers={"User-Agent": "SKEMPI-Downloader/1.0"},
                        )
                        response.raise_for_status()
                        with open(self.csv_path, "w") as f:
                            f.write(response.text)
                        print(f"Guardado en {self.csv_path}")
                        break
                    except requests.RequestException as exc:
                        last_error = exc
                        print(f"  Aviso: fallo de descarga ({exc})")
                        if attempt < 3:
                            time.sleep(1.5)
                else:
                    print(f"  Aviso: no se pudo descargar desde {url}")
                    continue
                break
            else:
                raise RuntimeError(
                    "No se pudo descargar SKEMPI v2.0 desde ninguna URL conocida. "
                    "Puedes descargar manualmente el CSV y guardarlo en "
                    f"{self.csv_path}"
                ) from last_error
        
        self.df = pd.read_csv(self.csv_path, sep=";")
        self.colmap = self._detect_columns(self.df)
        raw_pdb_col = self.colmap["pdb"]
        self.df["pdb_complex"] = self.df[raw_pdb_col].astype(str).str.upper()
        self.df["pdb_code"] = self.df["pdb_complex"].str.slice(0, 4)
        return self.df
    
    def list_all_complexes(self):
        """Lista todos los complejos únicos en SKEMPI"""
        if self.df is None:
            self.download_csv()
        
        complexes = self.df["pdb_code"].dropna().unique()
        print(f"\nTotal de complejos en SKEMPI v2.0: {len(complexes)}\n")
        
        # Agrupar por complejidad
        for pdb_id in sorted(complexes):
            subset = self.df[self.df["pdb_code"] == pdb_id]
            n_mutations = len(subset)
            print(f"  {pdb_id}: {n_mutations} mutations")
    
    def get_interface_residues_for_complex(self, pdb_id):
        """
        Extrae todas las posiciones de interfaz donde se han hecho mutaciones
        en un complejo SKEMPI dado
        """
        if self.df is None:
            self.download_csv()
        
        pdb_id = pdb_id.upper()
        subset = self.df[self.df["pdb_code"] == pdb_id]
        
        if len(subset) == 0:
            print(f"ERROR: PDB {pdb_id} no encontrado en SKEMPI")
            available = sorted(self.df["pdb_code"].dropna().unique())
            suggestions = difflib.get_close_matches(pdb_id, available, n=5, cutoff=0.5)
            if len(pdb_id) >= 2:
                same_prefix = [x for x in available if x.startswith(pdb_id[:2])]
                suggestions = list(dict.fromkeys(suggestions + same_prefix[:10]))
            if suggestions:
                print(f"Sugerencias cercanas: {', '.join(suggestions)}")
            return None
        
        mutations_info = []
        interface_positions = {'chain1': set(), 'chain2': set()}
        
        for _, row in subset.iterrows():
            mutation = row[self.colmap["mutation"]]
            aff_mut = float("nan")
            aff_wt = float("nan")
            if self.colmap["aff_mut"] is not None:
                aff_mut = self._to_nM(row[self.colmap["aff_mut"]], self.colmap["aff_mut"])
            if self.colmap["aff_wt"] is not None:
                aff_wt = self._to_nM(row[self.colmap["aff_wt"]], self.colmap["aff_wt"])
            affinity_change = aff_mut - aff_wt if not (pd.isna(aff_mut) or pd.isna(aff_wt)) else float("nan")

            if self.colmap["ddg"] is not None and pd.notna(row[self.colmap["ddg"]]):
                ddg = float(row[self.colmap["ddg"]])
            elif not (pd.isna(aff_mut) or pd.isna(aff_wt)) and aff_mut > 0 and aff_wt > 0:
                # ddG aprox. desde Kd_mut/Kd_wt (asumiendo 298.15 K)
                ddg = 0.0019872041 * 298.15 * math.log(aff_mut / aff_wt)
            else:
                ddg = float("nan")

            chain_tokens = str(row["pdb_complex"]).split("_")[1:]
            chain1_ids = set(chain_tokens[0]) if len(chain_tokens) >= 1 else set()
            chain2_ids = set(chain_tokens[1]) if len(chain_tokens) >= 2 else set()
            
            # Parse mutation format: "A123B" means Chain A, position 123, to residue B
            if pd.notna(mutation) and mutation != '-':
                mutations = str(mutation).split(',')
                
                for mut in mutations:
                    parsed = self._parse_mutation(mut)
                    if not parsed:
                        continue
                    chain, pos, wt_aa, mt_aa = parsed

                    mutations_info.append({
                        'chain': chain,
                        'position': pos,
                        'wt': wt_aa,
                        'mutant': mt_aa,
                        'ddG': ddg,
                        'affinity_change': affinity_change
                    })

                    if chain in chain1_ids:
                        interface_positions['chain1'].add(pos)
                    elif chain in chain2_ids:
                        interface_positions['chain2'].add(pos)
                    elif chain.upper() == 'A':
                        interface_positions['chain1'].add(pos)
                    else:
                        interface_positions['chain2'].add(pos)
        
        return {
            'pdb_id': pdb_id,
            'n_mutations': len(subset),
            'mutations': mutations_info,
            'interface_chain1': sorted(interface_positions['chain1']),
            'interface_chain2': sorted(interface_positions['chain2'])
        }
    
    def download_pdb(self, pdb_id, output_dir="./"):
        """Descarga un PDB desde RCSB"""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        pdb_path = output_dir / f"{pdb_id}.pdb"
        
        if pdb_path.exists():
            print(f"PDB {pdb_id} ya descargado")
            return pdb_path
        
        url = f"{PDB_URL}{pdb_id}.pdb"
        print(f"Descargando {url}...")
        
        response = requests.get(url)
        response.raise_for_status()
        
        with open(pdb_path, 'w') as f:
            f.write(response.text)
        
        print(f"Guardado en {pdb_path}")
        return pdb_path
    
    def export_for_rosetta(self, complex_info, output_dir="./"):
        """
        Exporta la información de interfaz en formatos útiles para Rosetta
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        pdb_id = complex_info['pdb_id']
        
        # 1. Lista de posiciones de interfaz (simple text)
        interface_file = output_dir / f"{pdb_id}_interface_positions.txt"
        with open(interface_file, 'w') as f:
            f.write("# Interface residues for " + pdb_id + "\n")
            f.write(f"# Chain A: {len(complex_info['interface_chain1'])} residues\n")
            f.write(f"# Chain D: {len(complex_info['interface_chain2'])} residues\n\n")
            
            f.write("Chain A interface positions:\n")
            f.write(" ".join(map(str, complex_info['interface_chain1'])) + "\n\n")
            
            f.write("Chain D interface positions:\n")
            f.write(" ".join(map(str, complex_info['interface_chain2'])) + "\n")
        
        # 2. Archivo con mutaciones para alanine scanning (Rosetta format)
        mut_file = output_dir / f"{pdb_id}_mutations.txt"
        with open(mut_file, 'w') as f:
            f.write("# Mutation file for Rosetta ddG protocols\n")
            f.write("# Format: total_residues\n")
            f.write("# position chain_id wt_aa mutant_aa\n\n")
            
            total_interface = len(complex_info['interface_chain1']) + len(complex_info['interface_chain2'])
            f.write(f"{total_interface}\n")
            
            for mut in complex_info['mutations']:
                f.write(f"{mut['position']} {mut['chain']} {mut['wt']} {mut['mutant']}\n")
        
        # 3. CSV con datos detallados
        mut_df = pd.DataFrame(complex_info['mutations'])
        mut_csv = output_dir / f"{pdb_id}_detailed_mutations.csv"
        mut_df.to_csv(mut_csv, index=False)
        
        print(f"\nArchivos exportados a {output_dir}:")
        print(f"  - {interface_file.name}")
        print(f"  - {mut_file.name}")
        print(f"  - {mut_csv.name}")
        
        return {
            'interface_file': interface_file,
            'mutation_file': mut_file,
            'csv_file': mut_csv
        }

def main():
    parser = argparse.ArgumentParser(description='SKEMPI Downloader and Interface Extractor')
    parser.add_argument('--pdb', help='PDB ID to download and extract (e.g., 1B27)')
    parser.add_argument('--list-all', action='store_true', help='List all SKEMPI complexes')
    parser.add_argument('--output', default='./', help='Output directory')
    parser.add_argument('--cache', default='./skempi_cache', help='SKEMPI cache directory')
    
    args = parser.parse_args()
    
    downloader = SKEMPIDownloader(cache_dir=args.cache)
    downloader.download_csv()
    
    if args.list_all:
        downloader.list_all_complexes()
    
    if args.pdb:
        pdb_id = args.pdb.upper()
        print(f"\n{'='*60}")
        print(f"Processing complex: {pdb_id}")
        print(f"{'='*60}\n")
        
        # Get interface info
        complex_info = downloader.get_interface_residues_for_complex(pdb_id)
        
        if complex_info:
            print(f"Found {complex_info['n_mutations']} mutations in SKEMPI")
            print(f"\nInterface positions - Chain A: {complex_info['interface_chain1']}")
            print(f"Interface positions - Chain D: {complex_info['interface_chain2']}")
            
            # Download PDB
            downloader.download_pdb(pdb_id, args.output)
            
            # Export for Rosetta
            downloader.export_for_rosetta(complex_info, args.output)
            
            print(f"\n✓ Ready for Rosetta ddG calculations!")
            print(f"✓ Use the files in {args.output} as input for your GA")

if __name__ == "__main__":
    main()
