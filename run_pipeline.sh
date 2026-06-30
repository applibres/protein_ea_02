#!/usr/bin/env bash
set -euo pipefail

# Uso:
# ./run_pipeline.sh /ruta/estructura.pdb A B
# ./run_pipeline.sh /ruta/estructura.pdb A B /ruta/out_dir

if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "Uso: $0 <pdb_path> <target_chain> <binder_chain> [out_dir]"
  exit 1
fi

PDB_PATH="$1"
TARGET_CHAIN="$2"
BINDER_CHAIN="$3"
OUT_DIR="${4:-/home/ariel/Desktop/protea02/repro_run}"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_PY="${ROOT_DIR}/reproduce_deep_dive_pipeline.py"
PROTEINMPNN_SCRIPT="${ROOT_DIR}/ProteinMPNN/protein_mpnn_run.py"

# Defaults (sobrescribibles por variables de entorno)
PROTEINMPNN_PYTHON="${PROTEINMPNN_PYTHON:-/home/ariel/venvs/mpnn_cpu/bin/python}"
COLABFOLD_PYTHON="${COLABFOLD_PYTHON:-/home/ariel/venvs/colabfold/bin/python}"

if [[ ! -f "${PDB_PATH}" ]]; then
  echo "Error: no existe PDB: ${PDB_PATH}"
  exit 1
fi

if [[ ! -f "${PIPELINE_PY}" ]]; then
  echo "Error: no existe script pipeline: ${PIPELINE_PY}"
  exit 1
fi

if [[ ! -f "${PROTEINMPNN_SCRIPT}" ]]; then
  echo "Error: no existe script ProteinMPNN: ${PROTEINMPNN_SCRIPT}"
  exit 1
fi

if [[ ! -x "${PROTEINMPNN_PYTHON}" ]]; then
  echo "Error: no existe python de ProteinMPNN: ${PROTEINMPNN_PYTHON}"
  echo "Tip: export PROTEINMPNN_PYTHON=/ruta/a/python"
  exit 1
fi

if [[ ! -x "${COLABFOLD_PYTHON}" ]]; then
  echo "Error: no existe python de ColabFold: ${COLABFOLD_PYTHON}"
  echo "Tip: export COLABFOLD_PYTHON=/ruta/a/python"
  exit 1
fi

echo "[INFO] PDB: ${PDB_PATH}"
echo "[INFO] Target chain: ${TARGET_CHAIN}"
echo "[INFO] Binder chain: ${BINDER_CHAIN}"
echo "[INFO] Output: ${OUT_DIR}"
echo "[INFO] ProteinMPNN python: ${PROTEINMPNN_PYTHON}"
echo "[INFO] ColabFold python: ${COLABFOLD_PYTHON}"

python "${PIPELINE_PY}" \
  --pdb "${PDB_PATH}" \
  --target-chain "${TARGET_CHAIN}" \
  --binder-chain "${BINDER_CHAIN}" \
  --proteinmpnn-script "${PROTEINMPNN_SCRIPT}" \
  --proteinmpnn-python "${PROTEINMPNN_PYTHON}" \
  --colabfold-python "${COLABFOLD_PYTHON}" \
  --out-root "${OUT_DIR}"
