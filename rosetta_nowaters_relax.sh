#!/usr/bin/env bash
set -euo pipefail

# Uso:
#   ./rosetta_nowaters_relax.sh <input_pdb> [output_dir] [docker_image]
#
# Ejemplos:
#   ./rosetta_nowaters_relax.sh scenarios/test07/1BRS.pdb
#   ./rosetta_nowaters_relax.sh scenarios/test07/1BRS.pdb scenarios/test07
#   ./rosetta_nowaters_relax.sh scenarios/test07/1BRS.pdb scenarios/test07 protein_ea_02

if [[ $# -lt 1 || $# -gt 3 ]]; then
  echo "Uso: $0 <input_pdb> [output_dir] [docker_image]"
  exit 1
fi

INPUT_PDB="$(realpath "$1")"
OUT_DIR="${2:-$(dirname "$INPUT_PDB")}"
DOCKER_IMAGE="${3:-protein_ea_02}"

if [[ ! -f "$INPUT_PDB" ]]; then
  echo "Error: no existe el PDB de entrada: $INPUT_PDB"
  exit 1
fi

mkdir -p "$OUT_DIR"
OUT_DIR="$(realpath "$OUT_DIR")"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_STEM="$(basename "${INPUT_PDB%.pdb}")"
NOWATERS_PDB="${OUT_DIR}/${INPUT_STEM}_nowaters.pdb"
RELAXED_PDB="${OUT_DIR}/${INPUT_STEM}_nowaters_relaxed.pdb"

echo "[1/2] Quitando moléculas de agua -> ${NOWATERS_PDB}"
# Elimina registros ATOM/HETATM cuyo residuo sea agua (HOH/WAT/DOD/TIP*)
awk '
{
  rec = substr($0, 1, 6)
  res = substr($0, 18, 3)
  if ((rec == "ATOM  " || rec == "HETATM") &&
      (res == "HOH" || res == "WAT" || res == "DOD" || res == "TIP")) {
    next
  }
  print
}' "$INPUT_PDB" > "$NOWATERS_PDB"

if command -v /usr/local/bin/relax.default.linuxgccrelease >/dev/null 2>&1; then
  echo "[2/2] Relajando con Rosetta local (/usr/local/bin/relax.default.linuxgccrelease)"
  /usr/local/bin/relax.default.linuxgccrelease \
    -database /usr/local/database \
    -s "$NOWATERS_PDB" \
    -relax:fast \
    -nstruct 1 \
    -out:path:all "$OUT_DIR" \
    -out:suffix "_relaxed" \
    -out:no_nstruct_label true >/dev/null
else
  echo "[2/2] Relajando con Rosetta en Docker (${DOCKER_IMAGE})"

  if [[ "$NOWATERS_PDB" != "$ROOT_DIR"* || "$OUT_DIR" != "$ROOT_DIR"* ]]; then
    echo "Error: para modo Docker, input/output deben estar dentro del proyecto:"
    echo "       $ROOT_DIR"
    exit 1
  fi

  REL_NOWATERS="${NOWATERS_PDB#${ROOT_DIR}/}"
  REL_OUTDIR="${OUT_DIR#${ROOT_DIR}/}"

  docker run --rm \
    -v "${ROOT_DIR}:/app" \
    -w /app \
    "${DOCKER_IMAGE}" \
    /usr/local/bin/relax.default.linuxgccrelease \
      -database /usr/local/database \
      -s "/app/${REL_NOWATERS}" \
      -relax:fast \
      -nstruct 1 \
      -out:path:all "/app/${REL_OUTDIR}" \
      -out:suffix "_relaxed" \
      -out:no_nstruct_label true >/dev/null
fi

if [[ -f "$RELAXED_PDB" ]]; then
  echo "Listo:"
  echo "  - No waters: $NOWATERS_PDB"
  echo "  - Relaxed:   $RELAXED_PDB"
  exit 0
fi

ALT_RELAXED="$(ls -1 "${OUT_DIR}/${INPUT_STEM}_nowaters_relaxed"*.pdb 2>/dev/null | head -n 1 || true)"
if [[ -n "$ALT_RELAXED" ]]; then
  echo "Listo:"
  echo "  - No waters: $NOWATERS_PDB"
  echo "  - Relaxed:   $ALT_RELAXED"
  exit 0
fi

echo "Error: no se encontró el PDB relajado esperado."
exit 1
