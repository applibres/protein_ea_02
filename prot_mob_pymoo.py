#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 16/04/2025

@author: Rolando Armas
Yachay Tech University
Phage Therapy Group

Simple Genetic Algorithm for Protein Mutation - MOEA/D
"""

from typing import List

from genetic_operators.crossoever_operators import llm_crossover, uniform_crossover
from genetic_operators.moead import SelectionOperatorMOEAD
from genetic_operators.individual import Individual
from prot_interface.prot_problemI import *
import prot_interface.prot_problemI as problem
import prot_interface.prot_settingsI as sets
import os
import shutil
import numpy as np
import random
from prot_interface.logging_config import setup_logging
from prot_interface.prot_esm2 import ESM2ProbMatrix
from utils import (
    csvToTree, read_scfiles, save_population_to_csv, save_HallofFame,
    mutation_labels, save_evolutions_statistics, save_population_aa, num2str
)
import logging
import pickle

setup_logging()
logger = logging.getLogger(__name__)


SCORE_OBJECTIVE_NAMES = {
    0: "packstat",
    1: "sc_value",
    2: "total_score",
    3: "delta_unsatHbonds",
    4: "fa_rep",
    5: "per_residue_energy_int",
    6: "dSASA_int",
    7: "dG_separated_per_dSASA",
    8: "hbonds_int",
}



class pymoo_sga_protein:

    def __init__(self, scenario, algoritm_params, sim_params, fitness_idxs, output, randomseed):
        self.algoritm_params = algoritm_params
        self.sim_params      = sim_params
        self.fitness_idxs    = fitness_idxs
        self.output          = output
        self.randomseed      = randomseed
        self.scenario        = scenario

        self.partners      = self.sim_params['partners']
        self.ligand_chain  = self.sim_params['ligand_chain']
        self.pdbfile       = self.sim_params['pdbfile']

        self.popsize      = self.algoritm_params['popsize']
        self.ngenerations = self.algoritm_params['gen']
        self.nobj         = len(self.fitness_idxs) + 1
        self.mutprob      = self.algoritm_params['mutp']

        self.n_obj = len(self.fitness_idxs) + 2

        # ── MOEA/D operators ─────────────────────────────────────────────────
        moead = SelectionOperatorMOEAD(self.n_obj, n_neighbors=5, n_partitions=6)
        self.moead          = moead                    # full object kept for crossover
        self.selection      = moead.selection          # survival update
        self.parent_select  = moead.parent_selection   # neighbourhood parent picker
        # ─────────────────────────────────────────────────────────────────────

        weights      = [sets.SCORE_OBJECTIVE[fidx] for fidx in self.fitness_idxs]
        self.weights = [*weights, 1, -1]

        self.my_protein_problem = problem.prot_problem(
            self.scenario, self.partners, self.ligand_chain
        )

        logging.info("Initializing ESM2 model...")
        self.esm2 = ESM2ProbMatrix()
        logging.info("ESM2 model initialized successfully")

        self.aa0 = self.my_protein_problem.create_individual0(self.pdbfile)
        self.aa0_complete = self.my_protein_problem.get_complete_interest_sequence(
            sets.CONFIG_PATH + self.scenario + "/" + self.pdbfile, self.ligand_chain
        )
        self.ll_father, _ = self.esm2.get_esm_ll(self.aa0_complete)

    # ------------------------------------------------------------------ #
    #  Path helpers                                                        #
    # ------------------------------------------------------------------ #

    def path_temp_file(self, generation: int, i: int):
        temp_dir = self.output + "/tmp"
        os.makedirs(temp_dir, exist_ok=True)
        return temp_dir + "/temp_g" + str(generation) + "_" + num2str(i) + ".pdb"

    def path_individual_file(self, generation: int, i: int):
        gen_dir = self.output + "/g" + str(generation)
        os.makedirs(gen_dir, exist_ok=True)
        return gen_dir + "/g" + str(generation) + "_" + num2str(i) + ".pdb"

    def mutation_labels_from_original(self, sequence) -> List[str]:
        return mutation_labels(sequence, self.aa0, self.my_protein_problem.aa_pos_list)

    def register_mutations_from_original(self, ind: Individual):
        ind.mutations = self.mutation_labels_from_original(ind.sequence())
        ind.nmut = len(ind.mutations)

    def count_mutations_from_original(self, sequence) -> int:
        return len(self.mutation_labels_from_original(sequence))

    # ------------------------------------------------------------------ #
    #  Core GA operations                                                  #
    # ------------------------------------------------------------------ #

    def create_initial_individual(self) -> Individual:
        dst = self.path_individual_file(0, 0)

        ind        = Individual(self.aa0)
        ind.pdb    = dst
        ind.nmut   = 0
        ind.mutations = []

        src_pdb = sets.CONFIG_PATH + self.scenario + "/" + self.pdbfile
        shutil.copy(src_pdb, dst)

        if not os.path.isfile(dst):
            exit(1)

        fitness_indv0  = self.my_protein_problem.fitness(dst)
        fitness_values = [*fitness_indv0, 0]
        ind.F = np.array([f * w for f, w in zip(fitness_values, self.weights)])

        return ind

    def evaluate(self, population: List[Individual]):
        population = sorted(
            population, key=lambda x: int(x.pdb.split("_")[-1].replace(".pdb", ""))
        )
        fitness = self.my_protein_problem.fitnessPop([ind.pdb for ind in population])

        for i, ind in enumerate(population):
            ind_complete_sequence = self.my_protein_problem.get_complete_interest_sequence(
                ind.pdb, self.ligand_chain
            )
            if len(population) == 1:
                delta_ll = 0
            else:
                ll, _    = self.esm2.get_esm_ll(ind_complete_sequence)
                delta_ll = ll - self.ll_father

            fitness_individual = [fitness[i][fitness_idx] for fitness_idx in self.fitness_idxs]
            fit   = [*fitness_individual, delta_ll, ind.nmut]
            ind.fitness = fit
            ind.F = np.array([f * (-1 * w) for f, w in zip(fit, self.weights)])

            logging.info(f"Individual fitness -> {ind.F}")
            logging.info(f"Rosetta scores     -> {fit}")

    def mutation(self, population: List[Individual], generation: int) -> List[Individual]:
        """Apply mutation to every individual in *population*."""
        arguments   = []
        individuals = []

        population = sorted(
            population, key=lambda x: int(x.pdb.split("_")[-1].replace(".pdb", ""))
        )

        for i, mutant in enumerate(population):
            output_file_path = self.path_temp_file(generation, i)
            arguments.append(
                (
                    mutant.id,
                    mutant.pdb,
                    output_file_path,
                    self.mutprob,
                    generation,
                    self.ngenerations,
                    self.aa0_complete,
                )
            )

        mutation_results = self.my_protein_problem.mutate_population(arguments)

        for result in mutation_results:
            ind        = Individual(result["sequence"])
            ind.pdb    = result["pdb"]
            ind.parent_id = result["parent_id"]
            self.register_mutations_from_original(ind)
            individuals.append(ind)

        logging.info(f"PDB Files: {arguments}")
        return individuals

    # ------------------------------------------------------------------ #
    #  MOEA/D crossover                                                    #
    # ------------------------------------------------------------------ #

    def crossover(self, population: List[Individual], generation: int, mut_start_offset: int = 0) -> List[Individual]:
        """
        Produce one child per sub-problem via neighbourhood crossover.

        For each sub-problem *i*:
          1. ``parent_selection`` picks (parent_a, parent_b) from neighbour hood *i*.
          2. ``CrossoverOperatorMOEAD.get_crossover_args`` computes which
             positions differ and which alleles the child inherits from parent_b.
          3. The sequence-level indices are mapped to absolute Rosetta positions
             through ``aa_pos_list``.
          4. ``prot_problem.crossover_population`` applies all structural changes
             in parallel via Rosetta.
          5. ``Individual`` shells are assembled with lineage metadata.

        Parameters
        ----------
        population       : current population
        generation       : current generation (used for temp-file naming)
        mut_start_offset : offset added to the file index to avoid name clashes
                           when crossover and mutation temps coexist.

        Returns
        -------
        List[Individual]  – one child per sub-problem (same length as population)
        """
        aa_pos_list = self.my_protein_problem.aa_pos_list  # absolute Rosetta positions

        crossover_args: List[tuple] = []   # args for prot_problem.crossover_population
        meta: List[dict] = []              # lineage info per child

        for i in range(len(population)):
            parent_a, parent_b = self.parent_select(population, i)[0]

            # Sequence-level indices + alleles for the child
            pdb_base = parent_a.pdb
            #try:
            #    seq_indices, child_aas = llm_crossover(
            #        parent_a=parent_a,
            #        parent_b=parent_b,
            #        population=population,
            #        sequence_initial="".join(self.aa0),
            #    )
            #except Exception as exc:
            #    logging.warning(
            #        "[Crossover] LLM crossover failed for sub-problem %s; ",
            #        "falling back to uniform crossover: %s",
            #        i, exc,
            #    )
            seq_indices, child_aas = uniform_crossover(parent_a, parent_b)

            
            # Build the actual output path
            output_file = self.path_temp_file(generation, i + mut_start_offset + len(population))

            # Map sequence-level (0-based) indices → absolute Rosetta positions
            rosetta_positions = [aa_pos_list[idx] for idx in seq_indices]

            crossover_args.append((pdb_base, output_file, rosetta_positions, child_aas))

            meta.append({"output_file": output_file})

            logging.info(
                f"[Crossover] sub-problem {i}: "
                f"parents=({parent_a.id}, {parent_b.id}), "
                f"positions_changed={len(seq_indices)}"
            )

        # Apply structural crossovers in parallel (Rosetta)
        new_seqs = self.my_protein_problem.crossover_population(crossover_args)

        # Assemble Individual objects
        children: List[Individual] = []
        for m, aa in zip(meta, new_seqs):
            child        = Individual(aa)
            child.pdb    = m["output_file"]
            self.register_mutations_from_original(child)
            children.append(child)

        return children

    def move_population(self, population: List[Individual], generation: int):
        population = sorted(
            population, key=lambda x: int(x.pdb.split("_")[-1].replace(".pdb", ""))
        )
        for i, ind in enumerate(population):
            src = ind.pdb
            dst = self.path_individual_file(generation, i)
            logging.info(f"Moving PDB file: {src} -> {dst}")

            if src == dst:
                continue

            shutil.copyfile(src, dst)
            ind.pdb = dst

            src_sc = src + ".sc"
            dst_sc = dst + ".sc"
            if os.path.isfile(src_sc):
                shutil.copyfile(src_sc, dst_sc)
                logging.info(f"Copying sc file: {src_sc} -> {dst_sc}")
            else:
                logging.warning(f"Expected .sc file not found: {src_sc}")

    # ------------------------------------------------------------------ #
    #  Pareto front & statistics                                           #
    # ------------------------------------------------------------------ #

    def update_pareto_front(self, pareto_front: List[Individual], population: List[Individual]) -> List[Individual]:
        combined = pareto_front + population
        F = np.array([ind.F for ind in combined])

        is_pareto = np.ones(len(combined), dtype=bool)
        for i, fi in enumerate(F):
            if not is_pareto[i]:
                continue
            dominated = np.all(F <= fi, axis=1) & np.any(F < fi, axis=1)
            dominated[i] = False
            is_pareto[dominated] = False
            if np.any(np.all(F[dominated] <= fi, axis=1)):
                is_pareto[i] = False

        return [ind for ind, keep in zip(combined, is_pareto) if keep]

    def _compile_stats(self, population: List[Individual]) -> dict:
        F = np.array([ind.F for ind in population])
        return {
            "avg": np.mean(F, axis=0).tolist(),
            "std": np.std(F, axis=0).tolist(),
            "min": np.min(F, axis=0).tolist(),
            "max": np.max(F, axis=0).tolist(),
        }

    # ------------------------------------------------------------------ #
    #  Main evolutionary loop                                              #
    # ------------------------------------------------------------------ #

    def run(self, checkpoint: bool = False, freq: int = 2):

        PATH_STATISTICS     = os.path.join(os.path.dirname(self.output), "statistics")
        SCFILE_CSV          = os.path.join(PATH_STATISTICS, f"run{self.randomseed}_scfile.csv")
        SAVEPOPGEN_CSV      = os.path.join(PATH_STATISTICS, f"run{self.randomseed}_individuals.csv")
        SAVEHALLOFFAME_CSV  = os.path.join(PATH_STATISTICS, f"run{self.randomseed}_hallofame.csv")
        SAVE_STATISTICS_CSV = os.path.join(PATH_STATISTICS, f"run{self.randomseed}_statistics.csv")
        CHECKPOINT_FILE     = os.path.join(self.output, "checkpoint.pkl")

        os.makedirs(PATH_STATISTICS, exist_ok=True)
        os.makedirs(self.output + "/tmp", exist_ok=True)

        random.seed(self.randomseed)

        # ── Checkpoint restore ──────────────────────────────────────────
        if checkpoint and os.path.isfile(CHECKPOINT_FILE):
            try:
                logging.info("Loading from checkpoint...")
                with open(CHECKPOINT_FILE, 'rb') as cp_file:
                    cp = pickle.load(cp_file)
                population   = cp['population']
                pareto_front = cp['pareto_front']
                ngen_start   = cp['generation'] + 1
                logbook      = cp['logbook']
                random.setstate(cp['rndstate'])
                logging.info(f"Resumed from generation {cp['generation']}")
            except (FileNotFoundError, EOFError) as e:
                logging.critical(f"Checkpoint error: {e}")
                exit(1)

        else:
            # ── Fresh start ─────────────────────────────────────────────
            logging.info("Starting new run")

            if os.path.exists(SAVEHALLOFFAME_CSV): os.remove(SAVEHALLOFFAME_CSV)
            if os.path.exists(SAVEPOPGEN_CSV):     os.remove(SAVEPOPGEN_CSV)

            ind0 = self.create_initial_individual()
            population: List[Individual] = [ind0]

            population += self.mutation(population * self.popsize, 0)

            self.evaluate(population)
            self.move_population(population, 0)

            pareto_front: List[Individual] = self.update_pareto_front([], population)

            record  = self._compile_stats(population)
            logbook = [{"gen": 0, **record}]

            save_population_to_csv(population, 0, SAVEPOPGEN_CSV)
            save_population_aa(
                self.output + "/g0", 0, population, [ind.pdb for ind in population]
            )
            logging.info(f"Generation 0 stats: {record}")

            ngen_start = 1

        # ── Evolutionary loop ───────────────────────────────────────────
        for generation in range(ngen_start, self.ngenerations):
            logging.info(f"-- Generation {generation} --")

            # Step 1: Parent selection + crossover (one child per sub-problem)
            logging.info(f"[Gen {generation}] Crossover...")
            children = self.crossover(population, generation, mut_start_offset=0)

            # Step 2: Mutation applied to the crossover children
            logging.info(f"[Gen {generation}] Mutation...")
            mutants = self._mutate_children(children, population, generation)

            # Step 3: Evaluate mutants
            logging.info(f"[Gen {generation}] Evaluation...")
            self.evaluate(mutants)

            # Step 4: MOEA/D survival update (Tchebycheff neighbourhood replacement)
            logging.info(f"[Gen {generation}] Selection (MOEA/D survival)...")
            population = self.selection(population, mutants)

            # Step 5: Move PDB files to permanent generation directory
            self.move_population(population, generation)

            pareto_front = self.update_pareto_front(pareto_front, population)

            record = self._compile_stats(population)
            logbook.append({"gen": generation, **record})
            logging.info(f"Generation {generation} stats: {record}")

            save_population_to_csv(population, generation, SAVEPOPGEN_CSV)
            save_population_aa(
                self.output + f"/g{generation}", generation,
                population, [ind.pdb for ind in population]
            )

            # Checkpoint
            if generation % freq == 0 or generation == self.ngenerations - 1:
                cp = dict(
                    population=population,
                    pareto_front=pareto_front,
                    generation=generation,
                    logbook=logbook,
                    rndstate=random.getstate(),
                )
                with open(CHECKPOINT_FILE, 'wb') as cp_file:
                    pickle.dump(cp, cp_file)
                logging.info(f"Checkpoint saved at generation {generation}")

            # Clean temp directory
            temp_dir = os.path.join(self.output, "tmp")
            for f in os.listdir(temp_dir):
                entry_path = os.path.join(temp_dir, f)
                if os.path.islink(entry_path) or os.path.isfile(entry_path):
                    os.remove(entry_path)
                elif os.path.isdir(entry_path):
                    shutil.rmtree(entry_path)
                else:
                    logging.warning(f"Skipping unknown temp entry type: {entry_path}")

            logging.info(f"-- End Generation {generation} --")

        # ── Post-evolution outputs ──────────────────────────────────────
        logging.info("-- End of Evolution --")

        csvToTree(SAVEPOPGEN_CSV, SAVEPOPGEN_CSV.replace('.csv', '.json'))
        save_HallofFame(pareto_front, SAVEHALLOFFAME_CSV)
        read_scfiles(self.output, SCFILE_CSV)

        logging.info("-- Saving Evolution Statistics --")
        save_evolutions_statistics(logbook, SAVE_STATISTICS_CSV, self.nobj)

    # ------------------------------------------------------------------ #
    #  Private helpers                                                     #
    # ------------------------------------------------------------------ #

    def _mutate_children(self, children: List[Individual], population: List[Individual], generation: int) -> List[Individual]:
        """
        Apply point mutation to a list of child individuals.

        Reuses the existing ``mutation`` logic but assigns file paths that
        avoid collisions with the crossover temp files (offset by popsize).
        """
        arguments   = []
        individuals = []

        for i, child in enumerate(children):
            # Offset by popsize so temp names don't clash with crossover files
            output_file_path = self.path_temp_file(generation, i + 2 * self.popsize)
            arguments.append(
                (
                    child.id,
                    child.pdb,
                    output_file_path,
                    self.mutprob,
                    generation,
                    self.ngenerations,
                    self.aa0_complete,
                )
            )

        mutation_results = self.my_protein_problem.mutate_population(arguments)

        for result in mutation_results:
            ind        = Individual(result["sequence"])
            ind.pdb    = result["pdb"]
            ind.parent_id = result["parent_id"]
            self.register_mutations_from_original(ind)
            individuals.append(ind)

        logging.info(f"[Mutation of children] PDB Files: {arguments}")
        return individuals
