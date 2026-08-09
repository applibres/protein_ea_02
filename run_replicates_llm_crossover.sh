#!/usr/bin/env bash
set -euo pipefail

# 8 runs total: 4 MOEA/D + 4 NSGA-III experiments x 1 replicate each
# Experiments:
#   MOEA/D: llm_crossover=1, all combos random_mutation(0/1) x bo_enabled(0/1)
#   NSGA-III: llm_crossover=1, all combos random_mutation(0/1) x bo_enabled(0/1)
#
# Outputs isolated by replicates-id (1..8).

SCENARIO="test06"
SIM_PARAMS="pdbfile=9Q1V_prepared_clean_relaxed.pdb,partners=A_B,ligand_chain=B"
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

run_experiment() {
  local algo="$1"
  local random_mutation="$2"
  local bo_enabled="$3"
  local llm_crossover="$4"
  local replicates_id="$5"

  local algo_params
  algo_params="gen=${GEN},popsize=${POPSIZE},n_neighbors=${N_NEIGHBORS},n_partitions=${N_PARTITIONS},random_mutation=${random_mutation},llm_crossover=${llm_crossover},bo_enabled=${bo_enabled},bo_candidates_per_parent=${BO_CANDIDATES},bo_beta=${BO_BETA},bo_min_train=${BO_MIN_TRAIN}"

  echo "============================================================"
  echo "Running experiment: algo=${algo} random_mutation=${random_mutation}, bo_enabled=${bo_enabled}, llm_crossover=${llm_crossover}, replicates-id=${replicates_id}"
  echo "============================================================"

  python replicates.py \
    --init "${INIT_SEED}" --end "${END_SEED}" \
    --scenario "${SCENARIO}" \
    --algo "${algo}" \
    --sim-params "${SIM_PARAMS}" \
    --algo-params "${algo_params}" \
    --fitness-idxs "${FITNESS_IDXS}" \
    --checkpoint false \
    --freq 2 \
    --mobj true \
    --replicates-id "${replicates_id}"
}

replicates_id=1

# --- MOEA/D experiments (replicates-id 1..4) ---
for random_mutation in 0 1; do
  for bo_enabled in 0 1; do
    run_experiment "moea" "${random_mutation}" "${bo_enabled}" "1" "${replicates_id}"
    replicates_id="$((replicates_id + 1))"
  done
done

# --- NSGA-III experiments (replicates-id 5..8) ---
for random_mutation in 0 1; do
  for bo_enabled in 0 1; do
    run_experiment "nsga3" "${random_mutation}" "${bo_enabled}" "1" "${replicates_id}"
    replicates_id="$((replicates_id + 1))"
  done
done

num_experiments="$((replicates_id - 1))"
num_replicates="$((END_SEED - INIT_SEED + 1))"
echo "All $((num_experiments * num_replicates)) runs completed (${num_experiments} experiments x ${num_replicates} replicates)."
