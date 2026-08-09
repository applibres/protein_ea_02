#!/usr/bin/env bash
set -euo pipefail

# 32 runs total: 8 experiments x 4 replicas (seeds 1..4)
# Experiments vary:
#   - random_mutation: 0/1
#   - bo_enabled: 0/1
#   - llm_crossover: 0/1
#
# Outputs are isolated by replicates-id:
#   exp1..exp8 -> replicates3..replicates10
#
# Usage:
#   ./run_replicates.sh [moead|moea|nsgaiii|nsga3]

SCENARIO="test04"
# Choose the multi-objective backend for all experiments in this batch.
# Argument default: moead
if [[ $# -gt 1 ]]; then
  echo "Usage: $0 [moead|moea|nsgaiii|nsga3]"
  exit 1
fi
ALGORITHM_BACKEND_INPUT="${1:-moead}"
SIM_PARAMS="pdbfile=protein01.pdb,partners=A_C,ligand_chain=C"
FITNESS_IDXS="fitness_idxs=7"
INIT_SEED=1
END_SEED=4
GEN=20
POPSIZE=28
N_NEIGHBORS=5
N_PARTITIONS=6
BO_CANDIDATES=8
BO_BETA=1.0
BO_MIN_TRAIN=100
REPLICATES_ID_START=1

ALGORITHM_BACKEND="$(echo "${ALGORITHM_BACKEND_INPUT}" | tr '[:upper:]' '[:lower:]')"
case "${ALGORITHM_BACKEND}" in
  moea|moead)
    ALGO="moea"
    ;;
  nsga3|nsgaiii)
    ALGO="nsga3"
    ;;
  *)
    echo "Invalid algorithm argument '${ALGORITHM_BACKEND_INPUT}'. Use: moead/moea or nsgaiii/nsga3."
    exit 1
    ;;
esac

run_experiment() {
  local random_mutation="$1"
  local bo_enabled="$2"
  local llm_crossover="$3"
  local replicates_id="$4"

  local algo_params
  algo_params="gen=${GEN},popsize=${POPSIZE},n_neighbors=${N_NEIGHBORS},n_partitions=${N_PARTITIONS},random_mutation=${random_mutation},llm_crossover=${llm_crossover},bo_enabled=${bo_enabled},bo_candidates_per_parent=${BO_CANDIDATES},bo_beta=${BO_BETA},bo_min_train=${BO_MIN_TRAIN}"

  echo "============================================================"
  echo "Running experiment: algo=${ALGO} random_mutation=${random_mutation}, bo_enabled=${bo_enabled}, llm_crossover=${llm_crossover}, replicates-id=${replicates_id}"
  echo "============================================================"

  python replicates.py \
    --init "${INIT_SEED}" --end "${END_SEED}" \
    --scenario "${SCENARIO}" \
    --algo "${ALGO}" \
    --sim-params "${SIM_PARAMS}" \
    --algo-params "${algo_params}" \
    --fitness-idxs "${FITNESS_IDXS}" \
    --checkpoint false \
    --freq 2 \
    --mobj true \
    --replicates-id "${replicates_id}"
}

replicates_id="${REPLICATES_ID_START}"
echo "Selected backend: ALGORITHM_BACKEND=${ALGORITHM_BACKEND} -> --algo ${ALGO}"

for llm_crossover in 0 1; do
  for random_mutation in 0 1; do
    for bo_enabled in 0 1; do
      run_experiment "${random_mutation}" "${bo_enabled}" "${llm_crossover}" "${replicates_id}"
      replicates_id="$((replicates_id + 1))"
    done
  done
done

num_experiments="$((replicates_id - REPLICATES_ID_START))"
num_replicates="$((END_SEED - INIT_SEED + 1))"
echo "All $((num_experiments * num_replicates)) runs completed (${num_experiments} experiments x ${num_replicates} replicates)."
