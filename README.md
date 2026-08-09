# protein_ea_02

Evolutionary Algorithms for Protein Design

Phage Therapy Group Yachay Tech

## Introduction
This project implements **multi-objective evolutionary algorithms (MOEA)** to design protein-protein interfaces. Given a wild-type protein complex, the algorithm searches for point mutations on the interface of a target (ligand) chain that improve the binding interaction. Two MOEA backends are provided:

- **MOEA/D** (`prot_mob_pymoo.py`) — neighbourhood-decomposition-based search.
- **NSGA-III** (`prot_mob_nsga3_pymoo.py`) — reference-direction niching + non-dominated sorting.

Both backends share the same genetic operators (crossover, mutation, evaluation, relaxation) and differ only in how parents are selected and how the population survives across generations.

## Installation
You need to install the software listed below:

1) Create docker image from Dockerfile using

```bash
docker build -t protein_ea_02 .
```

2) Run the docker image

```bash
docker run -it --rm protein_ea_02
```

MultiObjective mode makes use of ESM2 model, so if your device has GPU, you can use the next command:

```bash
docker run -it --rm --gpus all protein_ea_02
```

Alternatively, `docker-compose.yml` starts an **ollama** service (used for LLM-guided crossover) alongside the `protein-rosetta` container.

3) Change the respective parameters values inside file prot_interface/prot_settingsI.py. The following parameters need to set according to your local repository and rosetta-commons instalation:

```
ROSETTA_BIN = Main directory path wher is rosetta-commons binary files. ("/rosetta.binary.m1.release-371/main/source/bin/")
CONFIG_PATH = Main directory where are located the scenarios to test ("/scenarios/")
INTERF_EN = Name of interface_energy binary file ("interface_energy.static.macosclangrelease -s")
FACE1_FILE_NAME = Face A name file ("faceA.txt")
FACE2_FILE_NAME = Face C name file ("faceC.txt")
MSA_MATRIX = Name of MSA Matrix ("MSA_matrix.tsv")
INTERF_AN = Name of InterfaceAnalyzer binary file ("InterfaceAnalyzer.static.macosclangrelease -s")
SCORE_INDEXES = [4,7,5,11,23] #score positions for dG_separated, dSASA_int, dG_separated/dSASAx100, hbonds_int
FLAGS = "-compute_packstat true -tracer_data_print false -pack_input true " \
        "-pack_separated true -add_regular_scores_to_scorefile true " \
        "-atomic_burial_cutoff 0.01 -sasa_calculator_probe_radius 1.4 " \
        "-pose_metrics::interface_cutoff 8.0 -use_input_sc " \
        "-out:file:score_only "
```

## Basic Use
You can download the project's content as a zip file. Uncompressed it in the corresponding location that you decide to copy the project's content. The directory contains the next structure:

- `prot_interface/`: Python classes that interface with PyRosetta, Rosetta binaries, ESM2 and the BO surrogate.
- `genetic_operators/`: Crossover operators and the MOEA/D / NSGA-III selection-survival operators.
- `scenarios/`: Folder with the protein file in pdb format plus per-scenario `faceA.txt`/`faceB.txt` files.
- `prot_mob_pymoo.py`: Multi-objective evolutionary algorithm (MOEA/D backend).
- `prot_mob_nsga3_pymoo.py`: Multi-objective evolutionary algorithm (NSGA-III backend).
- `eaprot_call.py`: Main python script that parses the evolutionary algorithm parameters and calls the selected backend.
- `replicates.py`: Sequential replicate runner that isolates outputs per experiment/replicate set.
- `run_replicates.sh` / `run_replicates_llm_crossover.sh`: Batch scripts that sweep `random_mutation` × `bo_enabled` × `llm_crossover`.

Example for multiobjective NSGA-III backend with `replicates.py`:

```bash
python replicates.py \
  --init 1 --end 1 \
  --scenario test06 \
  --algo nsga3 \
  --sim-params "pdbfile=9Q1V_prepared_clean_relaxed.pdb,partners=A_B,ligand_chain=B" \
  --algo-params "gen=5,popsize=28,n_partitions=6,random_mutation=1,llm_crossover=0,bo_enabled=1,bo_candidates_per_parent=8,bo_beta=1.0,bo_min_train=100" \
  --fitness-idxs "fitness_idxs=7" \
  --checkpoint false \
  --freq 2 \
  --mobj true \
  --replicates-id 1
```

### Parameters description

- `scenario`: name of the scenario folder (contains `faceA.txt`, `faceB.txt`, `protein.pdb`, etc.).
- `algo`: evolutionary backend name. For multiobjective mode (`mobj=True`) supported values are:
  - `moea` / `moead` (MOEA/D backend)
  - `nsga3` / `nsgaiii` (NSGA-III backend)
- `pdbfile`: name of the pdb file.
- `partners`: the chains that conform the protein structure.
- `ligand_chain`: the chain that will be used for evolution process.
- `gen`: number of generations.
- `popsize`: population size.
- `n_partitions`: Das-Dennis partitions for reference directions in MOEA/D and NSGA-III. For NSGA-III, `popsize` must be less than or equal to the number of generated reference directions.
- `n_neighbors` (MOEA/D only): neighbourhood size (T).
- `nr` (MOEA/D only): maximum number of neighbour replacements per offspring.
- `random_mutation`: `0` uses Rosetta+ESM2 candidate proposal, `1` uses random position/AA proposal.
- `llm_crossover`: `0` uses uniform crossover, `1` enables LLM-guided crossover (fallback to uniform on failure).
- `bo_enabled`: `1` activates the Bayesian-optimization surrogate to pick the best mutation candidate per child.
- `bo_candidates_per_parent`: number of mutation candidates generated per child (upper bound; actual = `top_positions × top_aa`).
- `bo_beta`: exploration/exploitation trade-off of the LCB acquisition function.
- `bo_min_train`: minimum number of observations before the surrogate is fitted/used.
- `fitness_idxs`: score positions for dG_separated, dSASA_int, dG_separated/dSASAx100, hbonds_int (indexes into the filtered `SCORE_INDEXES` list).
- `fitness_weights`: weights for each fitness function (-1 for minimize or 1 for maximize).
- `checkpoint`: restore from a previous checkpoint.
- `freq`: checkpoint frequency.
- `mobj`: algorithm mode. If `mobj=True` uses `fitness_idxs` scores and the ll score from esm2.
- `randomseed`: random seed.

## Algorithm description

### Genotype / individual
An individual is the **interface sequence of the ligand chain**: the residues that appear in the pairwise short-range energy output of the Rosetta `interface_energy` binary (mapped from relative to absolute residue numbers and stored in `aa_pos_list`). Each individual stores the amino-acid sequence as an integer array, the PDB path of its (relaxed) structure, its fitness vector, the number of mutations with respect to the wild type (`nmut`) and the mutation labels (e.g. `A45V`).

### Objectives (3 by default)
With the default `fitness_idxs=7`, the algorithm optimizes **3 objectives** (`n_obj = len(fitness_idxs) + 2`):

1. **dG_separated/dSASAx100** (Rosetta `InterfaceAnalyzer`) — minimised. Interaction energy per buried surface area, the main proxy for binding affinity. This is also the target used to train the BO surrogate.
2. **ΔLL (ESM2)** — maximised. `LL(child) - LL(wild-type)` computed by the ESM2 masked language model on the complete ligand-chain sequence. Acts as a sequence-plausibility / evolutionary-fitness regulariser.
3. **nmut** — minimised. Number of point mutations relative to the wild type; penalises drift away from the native sequence.

All objectives are internally converted to a **minimisation space** `F = fitness × (-1 × weight)` before being handed to the selection operators.

### Generation 0 (initial population)
1. The wild-type PDB is copied, evaluated (relaxed structure is scored with Rosetta) and becomes `individual 0`.
2. `popsize - 1` mutants are produced from the wild type through the **mutation operator** (one single-point mutation each, ESM2-informed or random), each relaxed and evaluated.
3. The full population is moved to `g0/`, a Pareto front (hall-of-fame) is initialised, and statistics are logged.

### Main evolutionary loop (per generation)

**Step 1 — Parent selection + crossover.** One child is produced per population slot:

- **MOEA/D**: parents are sampled uniformly at random from the *neighbourhood* of the sub-problem's reference-direction weight vector.
- **NSGA-III**: parents are picked globally by **binary tournament on non-dominated rank** (whole population mating pool).

The crossover operator computes, per interface position, which alleles the child inherits from parent B:

- **Uniform crossover**: every position where the parents differ is taken from parent B with probability 0.5.
- **LLM-guided crossover** (`llm_crossover=1`): a local LLM (`gemma4:e4b` via ollama) receives both parent sequences, their `dG_separated/dSASAx100` fitness and the whole population, and returns a child interface sequence. On any failure (invalid length/AA, identical to a parent, timeout) it falls back to uniform crossover.

The child is built on top of parent A's structure; the resulting PDB is obtained either from `sequence_pdb_cache` (sequence already relaxed in a previous evaluation) or by applying the residue substitutions through Rosetta (`MutateResidue` + local `FastRelax`) in parallel.

**Step 2 — Mutation of every child.** For each child, the mutation operator proposes up to `top_positions × top_aa` **single-point mutation candidates** (default 8×8 = 64):

1. Candidate positions are extracted from the Rosetta `interface_energy` output (`aa_stab_nstab_list`), which ranks interface residues by their per-residue interaction energy. In informed mode (`random_mutation=0`) the top `top_positions` are taken; in random mode (`random_mutation=1`) positions are shuffled.
2. For each candidate position, candidate amino acids are proposed:
   - **Informed (ESM2)**: `top_k_replacements` uses the ESM2 masked-language model probability at that position, with a temperature that **anneals across generations** (`T = 3·(1 − (gen+1)/ngen)`, min 0.1), so exploration decreases over time. The wild-type residue is excluded.
   - **Random**: `top_aa` random non-wild-type residues.
3. Each candidate carries a `candidate_score` = `−(pos_rank + aa_rank)` used for softmax sampling.

**Step 3 — Candidate selection (one survivor per child).** Exactly one candidate per child becomes the actual mutant:

- **Bayesian optimisation path** (`bo_enabled=1` and the surrogate has been trained, i.e. `n_observations ≥ bo_min_train`): candidates are embedded with ESM2 (pooled hidden-state embeddings), reduced with PCA (64 components) and predicted by a **Gaussian-Process surrogate** (constant × RBF + white kernel) trained on real `dG_separated/dSASAx100` evaluations. The candidate with the best **LCB acquisition** `μ − β·σ` is chosen. Duplicate complete sequences are skipped globally.
- **Softmax path** (otherwise): candidates are sampled with a softmax over `candidate_score`, favouring the top-ranked position/AA proposals.

**Step 4 — Relaxation & real evaluation.** Each selected candidate is applied to the parent structure with Rosetta's `PM_Mutation_Relax_Local.xml` protocol (point mutation + cartesian FastRelax with coordinate constraints on a local neighbourhood around the mutated residue, `ref2015_cart`). Identical candidate sequences are deduplicated so each unique structure is relaxed only once; results are cached in `sequence_pdb_cache`. Real fitness is then computed with Rosetta `InterfaceAnalyzer` + ESM2 ΔLL, and the observations feed the BO surrogate (refitted when `≥ bo_min_train`).

**Step 5 — Survival selection.**

- **MOEA/D**: classical **Tchebycheff** neighbourhood replacement. The ideal point `z*` is updated incrementally; each offspring is compared against the individuals in its neighbourhood and replaces those whose Tchebycheff scalarized score is worse, up to `nr` replacements per offspring (duplicate winners are deep-copied to keep one PDB per slot).
- **NSGA-III**: the combined population + offspring (2·popsize) is filtered by pymoo's NSGA-III survival (non-dominated sorting + niching around the Das-Dennis reference directions), keeping exactly `popsize` individuals.

**Step 6 — Bookkeeping.** PDBs are moved to `g<gen>/`, the Pareto front / hall-of-fame is updated, per-generation statistics and population CSVs are written, the temp directory is cleaned, and a **checkpoint** (population, Pareto front, logbook, RNG states) is saved every `freq` generations.

### After the run
- `statistics/run<seed>_statistics.csv` — per-generation avg/min/max/std of every objective.
- `statistics/run<seed>_individuals.csv` + `.json` — per-generation individuals, fitness, mutations and sequences.
- `statistics/run<seed>_hallofame.csv` — final non-dominated (Pareto) front.
- `statistics/run<seed>_scfile.csv` — raw Rosetta `.sc` scores per individual.
- `g<gen>/` — relaxed PDBs of every surviving individual per generation plus `pop_g<gen>_AA.txt` with the AA sequences.

## Remarks / details

- **Multiprocessing**: fitness, crossover and mutation run in parallel worker pools with PyRosetta initialised per worker (`spawn` start method). ESM2 is lazily instantiated inside workers to avoid pickling heavy torch state.
- **Sequence cache**: any complete ligand-chain sequence that has already been relaxed (e.g. from a checkpoint or an earlier generation) is reused by copying its PDB, saving expensive relaxations.
- **ESM2 temperature schedule**: mutation exploration is annealed from `T=3.0` down to `T=0.1` over the generations, sharpening the sampling distribution as the run progresses.
- **Reference directions**: both backends generate Das-Dennis weight vectors and, when `popsize` is smaller than the generated set, select a well-spread subset via greedy farthest-point sampling.
