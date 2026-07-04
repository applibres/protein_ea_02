#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 22/03/2025

@author: Rolando Armas
Yachay Tech University
Phage Therapy Group

Problems classes for protein design problem

Single-objective
- fitness: dGSeparated/dSASA * 100


Multi-objective
- fitness: Not defined yet


"""
import os
import re
import shutil
import random
import numpy as np
import pyrosetta
from multiprocessing import Pool
import multiprocessing as mp
import prot_interface.prot_pyrosettaI as ppyrst
import prot_interface.prot_energyInterfI as prot_en_intf
import prot_interface.prot_mutI as p_mut
import prot_interface.prot_pyrosettaI as ppyrst
import prot_interface.prot_aa_stI as prot_aa
import prot_interface.prot_settingsI as sets
from prot_interface.logging_config import setup_logging
import logging
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1


# Initialize logging before anything else
setup_logging()
logger = logging.getLogger(__name__)

def _init_pyrosetta_worker():
    """Initialize PyRosetta in each worker process"""
    import pyrosetta
    pyrosetta.init(
        "-nstruct 1 "
        "-ignore_zero_occupancy false "
        "-ex1 -ex2 "
        "-use_input_sc "
        "-flip_HNQ "
        "-no_optH false "
        "-mute all"  # Suppress output in workers
    )

class prot_problem:
    """ Protein Design Problem - Class to mutate and evaluate protein chain 
 
    """
    def __init__(self, scenario, partners, ligand_chain):
        """Constructor
        
        Parameters
        ----------
        sname: string
            scenario name where the relaxed protein chain base is located
        """
        
        self.scenario = scenario #scenario name

        #Instantiate Protein Energy Class
        self.protEn = prot_en_intf.prot_energyInterf(scenario)

        #Instantiate Mutation Operator Class
        self.mut = p_mut.prot_mut(scenario)

        #Instantiate Protein Get AA stabilized/no stabilized list
        self.list_aa = prot_aa.prot_aa_extract(scenario,ligand_chain)

        #Partners String
        self.partners = partners

        #Partners String
        self.ligand_chain = ligand_chain        

        #Config Path
        self.config_path = sets.CONFIG_PATH + self.scenario + "/"

        #AA absolute position list
        self.aa_pos_list = []

        #Output Path

    def extract_mappings(self, energy_filepath):
        """
        Maps a relative position from chain C to absolute position
        Extracts a dictionary mapping Number1 to Number2 from a text file with
        the format 'AAA###_C(###)', where:
        - 'AAA' is a 3-letter amino acid code.
        - '###' (Number1) is the residue number.
        - '_C(###)' (Number2) is the mapped residue number.

        Args:
            energy_filepath (str): Path to the input file.

        Returns:
            dict: A dictionary mapping Number1 to Number2.
        """
        mapping = {}

        with open(energy_filepath, 'r') as file:
            for line in file:
                match = re.search(r"([A-Z]{3})(\d+)_" + re.escape(self.ligand_chain) + "\((\d+)\)", line)
                if match:
                    number1 = int(match.group(2))  # Extract Number1
                    number2 = int(match.group(3))  # Extract Number2
                    mapping[number1] = number2  # Store in dictionary

        return dict(sorted(mapping.items()))

    def absolute_AA_positions(self, AA_map_positions):
        # Extract values as a list
        values_list = list(AA_map_positions.values())

        # Create new dictionary where key = value from list, value = index position
        indexed_dict = {value: index for index, value in enumerate(values_list)}
        return indexed_dict

   
    #Given a pdbfile and absolute positions list - return a string sequence    
    def get_individual_seq(self, pdb_file):
        # Generate the sequence and initialize fixed values based on it
        sequence = self.sequence(pdb_file)
        # Si no hay posiciones definidas (posible después de checkpoint), reconstruirlas
        if not self.aa_pos_list:
            logging.warning("aa_pos_list is empty. Recomputing from energy file...")
            energy_filepath = self.list_aa.energy_interact_file(pdb_file)
            aa_pos_dict = self.extract_mappings(energy_filepath)
            self.aa_pos_list = list(aa_pos_dict.values())
            logging.debug(f"Recomputed aa_pos_list: {self.aa_pos_list}")
    
        # Extraer aminoácidos en las posiciones definidas
        aa = []
        for aa_pos in self.aa_pos_list:
            if aa_pos <= len(sequence):
                aa.append(sequence[aa_pos - 1])
            else:
                logging.error(f"Invalid position {aa_pos} for sequence of length {len(sequence)}")
                aa.append('X')  # placeholder in caso de error

        return aa




    #Initialize an individual with aa list taking ligand chain (face) absolute position
    def create_individual0(self, pdbfile):        
        #1 Generate Energy Interaction File
        pdb_file_path = self.config_path + pdbfile
        energy_filepath = self.list_aa.energy_interact_file(pdb_file_path)

        logging.info("energy_file = %s", energy_filepath)

        #2 Map Positions
        try:
            aa_pos_dict=self.extract_mappings(energy_filepath)
        finally:
            if os.path.isfile(energy_filepath):
                os.remove(energy_filepath)
        logging.debug(aa_pos_dict)

        # Extract absolute position values as a list
        self.aa_pos_list = list(aa_pos_dict.values())
        logging.debug(f"List absolute position values: {self.aa_pos_list}")


        # Create new dictionary where key = value from list, value = index position
        indexed_dict = {value: index for index, value in enumerate(self.aa_pos_list)}
        logging.debug(indexed_dict)

        # Generate the sequence and initialize fixed values based on it
        sequence = self.sequence(pdb_file_path)
        
        #aminoacid list
        aa = []
        for aa_pos in self.aa_pos_list:
            aa.append(sequence[aa_pos-1])

        return aa       
        

    def fitness(self, pdb_file):
        """Individual fitness function 
        
        Parameters
        ----------
        pose_file: mutated chain pdb file name 
            
        
        Returns
        -------
        Tuple 
            Fitness values
        """
        
        #Evaluate the solution
        # f[0]: packstat -> MAXIMIZAR (Calidad del empaquetamiento; valor ideal > 0.65)
        # f[1]: sc_value -> MAXIMIZAR (Complementariedad de formas; valor ideal > 0.60) 
        # f[2]: total_score -> MINIMIZAR (Estabilidad global del complejo) [3, 4]
        # f[3]: delta_unsatHbonds -> MINIMIZAR (Penalización por polares no satisfechos; ideal cercano a 0) [5]
        # f[4]: fa_rep -> MINIMIZAR (Repulsión estérica/choques; debe mantenerse bajo para ser físicamente posible) 
        # f[5]: per_residue_energy_int -> MINIMIZAR (Energía promedio por residuo en la interfaz) 
        # f[6]: dSASA_int -> MAXIMIZAR (Área enterrada; valor ideal entre 1200-2000 A^2)
        # f[7]: dG_separated/dSASAx100 -> MINIMIZAR
        f = self.protEn.getEnergyInterf(pdb_file)
        return f

    ## Evaluate Population in parallel
    def fitnessPop(self, pop_pdb_files):
        """Fitness function 
        
        Parameters
        ----------
        pop_pdb_files: Tuple
            pdb files population 
        
        Returns
        -------
        Tuple 
            Evaluated population
        """  
 
        # Create pool with PyRosetta initializer
        with mp.Pool(processes=mp.cpu_count(), initializer=_init_pyrosetta_worker) as pool:
            evaluation = []
            
            # Calling pool for run in parallel
            # Fitness:
            # f[0]: packstat -> MAXIMIZAR (Calidad del empaquetamiento; valor ideal > 0.65)
            # f[1]: sc_value -> MAXIMIZAR (Complementariedad de formas; valor ideal > 0.60) 
            # f[2]: total_score -> MINIMIZAR (Estabilidad global del complejo) [3, 4]
            # f[3]: delta_unsatHbonds -> MINIMIZAR (Penalización por polares no satisfechos; ideal cercano a 0) [5]
            # f[4]: fa_rep -> MINIMIZAR (Repulsión estérica/choques; debe mantenerse bajo para ser físicamente posible) 
            # f[5]: per_residue_energy_int -> MINIMIZAR (Energía promedio por residuo en la interfaz) 
            # f[6]: dSASA_int -> MAXIMIZAR (Área enterrada; valor ideal entre 1200-2000 A^2)
            # f[7]: dG_separated/dSASAx100 -> MINIMIZAR
            
            for result_f in pool.map(self.fitness, pop_pdb_files):
                # return the fitness solution
                evaluation.append(result_f)
        
        return evaluation 


    ## Mutate Population in parallel
    def mutate_population(self, args, bo_context=None):
        """Mutate population in parallel
        
        Parameters
        ----------
        args: list of tuples
            Each tuple contains
            (parent_id, pdb_file, output_file, generation, ngen, original_sequence)
        
        Returns
        -------
        list
            Population of mutated individuals as dicts:
            {"sequence": list[str], "pdb": str, "parent_id": str}
        """
        plans = self.plan_mutation_population(args, bo_context=bo_context)
        selected_plans = self.select_mutation_plans(plans, bo_context=bo_context)
        sequence_pdb_cache = {}
        if bo_context is not None:
            sequence_pdb_cache = bo_context.get("sequence_pdb_cache", {}) or {}

        resolved_results = [None] * len(selected_plans)
        unresolved_plans = []
        unresolved_indices = []

        # Global reuse: if this complete sequence was already relaxed before, just copy its PDB.
        for idx, plan in enumerate(selected_plans):
            key = plan.get("mutated_complete_sequence", plan["mutated_interface_sequence"])
            cached_pdb = sequence_pdb_cache.get(key)
            target_pdb = plan["pdb_out"]

            if cached_pdb and os.path.isfile(cached_pdb):
                if os.path.abspath(cached_pdb) != os.path.abspath(target_pdb):
                    shutil.copy2(cached_pdb, target_pdb)
                sequence_pdb_cache[key] = target_pdb
                resolved_results[idx] = {
                    "sequence": self.get_individual_seq(target_pdb),
                    "pdb": target_pdb,
                    "parent_id": plan["parent_id"],
                }
            else:
                unresolved_indices.append(idx)
                unresolved_plans.append(plan)

        if unresolved_plans:
            unique_plans, index_map = self.deduplicate_mutation_plans(unresolved_plans)
            unique_results = self.apply_mutation_population(unique_plans)

            for local_idx, unique_idx in enumerate(index_map):
                global_idx = unresolved_indices[local_idx]
                plan = unresolved_plans[local_idx]
                canonical = unique_results[unique_idx]

                source_pdb = canonical["pdb"]
                target_pdb = plan["pdb_out"]

                if os.path.abspath(source_pdb) != os.path.abspath(target_pdb):
                    shutil.copy2(source_pdb, target_pdb)

                key = plan.get("mutated_complete_sequence", plan["mutated_interface_sequence"])
                sequence_pdb_cache[key] = target_pdb
                resolved_results[global_idx] = {
                    "sequence": canonical["sequence"],
                    "pdb": target_pdb,
                    "parent_id": plan["parent_id"],
                }

        return resolved_results

    def plan_mutation_population(self, args, bo_context=None):
        """Plan mutation jobs in parallel using deterministic top-position/top-aa candidates."""
        top_positions = 8 if bo_context is None else int(bo_context.get("top_positions", 8))
        top_aa = 8 if bo_context is None else int(bo_context.get("top_aa", 8))
        expanded_args = []
        for parent_slot, arg in enumerate(args):
            expanded_args.append((parent_slot, *arg, top_positions, top_aa))

        with mp.Pool(processes=mp.cpu_count(), initializer=_init_pyrosetta_worker) as pool:
            plans_per_slot = pool.starmap(self._plan_mutation_worker_with_slot, expanded_args)

        plans = []
        for slot_plans in plans_per_slot:
            plans.extend(slot_plans)
        return plans

    @staticmethod
    def _softmax_sample(candidates):
        logits = np.array([float(c.get("candidate_score", 0.0)) for c in candidates], dtype=np.float64)
        logits = np.nan_to_num(logits, nan=0.0, posinf=0.0, neginf=0.0)
        logits = logits - np.max(logits)
        probs = np.exp(logits)
        probs_sum = probs.sum()
        if probs_sum <= 0:
            probs = np.ones(len(candidates), dtype=np.float64) / len(candidates)
        else:
            probs = probs / probs_sum
        idx = np.random.choice(len(candidates), p=probs)
        return candidates[int(idx)]

    def select_mutation_plans(self, plans, bo_context=None):
        """Select exactly one candidate plan per parent slot."""
        if not plans:
            return []

        grouped = {}
        for plan in plans:
            slot = plan["parent_slot"]
            grouped.setdefault(slot, []).append(plan)

        selected_plans = []
        bo_enabled = bo_context is not None and bo_context.get("enabled", False)
        surrogate = None if bo_context is None else bo_context.get("surrogate")
        beta = 1.0 if bo_context is None else float(bo_context.get("beta", 1.0))
        min_train = 24 if bo_context is None else int(bo_context.get("min_train", 24))
        bo_ready = bo_enabled and surrogate is not None and surrogate.is_ready(min_train)
        used_sequences = set()

        for slot in sorted(grouped.keys()):
            candidates = grouped[slot]
            available_candidates = [
                c for c in candidates
                if c.get("mutated_complete_sequence", c.get("mutated_interface_sequence")) not in used_sequences
            ]
            pool = available_candidates if available_candidates else candidates

            if bo_ready:
                sequences = [
                    plan.get("mutated_complete_sequence", plan["mutated_interface_sequence"])
                    for plan in pool
                ]
                acquisition = surrogate.acquisition(sequences, beta=beta)
                best_idx = min(range(len(pool)), key=lambda i: float(acquisition[i]))
                selected = pool[best_idx]
            else:
                selected = self._softmax_sample(pool)

            selected_plans.append(selected)
            used_sequences.add(selected.get("mutated_complete_sequence", selected.get("mutated_interface_sequence", "")))

        return selected_plans

    def deduplicate_mutation_plans(self, plans):
        """Deduplicate plans globally by mutated complete sequence."""
        unique_plans = []
        index_map = []
        key_to_unique_idx = {}

        for plan in plans:
            key = plan.get("mutated_complete_sequence", plan["mutated_interface_sequence"])
            if key in key_to_unique_idx:
                idx = key_to_unique_idx[key]
            else:
                idx = len(unique_plans)
                key_to_unique_idx[key] = idx
                unique_plans.append(plan)
            index_map.append(idx)

        return unique_plans, index_map

    def apply_mutation_population(self, unique_plans):
        """Apply unique mutation plans in parallel."""
        with mp.Pool(processes=mp.cpu_count(), initializer=_init_pyrosetta_worker) as pool:
            results = pool.map(self._apply_mutation_plan_worker, unique_plans)
        return results

    def _plan_mutation_worker(self, parent_id, pdb_file, output_file, generation, ngen, original_sequence,
                              top_positions, top_aa):
        sequence = "".join(self.get_complete_interest_sequence(pdb_file, self.ligand_chain))
        return self.mut.plan_mutation_candidates(
            scenario=self.scenario,
            ligand_chain=self.ligand_chain,
            pdb_file=pdb_file,
            output_file=output_file,
            sequence=sequence,
            generation=generation,
            ngen=ngen,
            original_sequence=original_sequence,
            parent_id=parent_id,
            top_positions=top_positions,
            top_aa=top_aa,
        )

    def _plan_mutation_worker_with_slot(self, parent_slot, parent_id, pdb_file, output_file, generation, ngen,
                                        original_sequence, top_positions, top_aa):
        plans = self._plan_mutation_worker(
            parent_id, pdb_file, output_file, generation, ngen, original_sequence, top_positions, top_aa
        )
        for plan in plans:
            plan["parent_slot"] = parent_slot
        return plans

    def _apply_mutation_plan_worker(self, plan):
        result = self.mut.apply_mutation_plan(plan)
        aa = self.get_individual_seq(result["pdb"])
        result["sequence"] = aa
        return result
    

    def _crossover_worker(self, pdb_file, output_file, positions_to_mutate, aminoacids_to_place):
        """
        Worker wrapper for parallel crossover.
        Calls prot_mut.crossover() and returns the resulting sequence.
    
        Parameters
        ----------
        pdb_file : str
            Structural base PDB (parent_a).
        output_file : str
            Path for the output crossover PDB.
        positions_to_mutate : List[int]
            Rosetta residue positions to apply from parent_b.
        aminoacids_to_place : List[str]
            Single-letter AA codes to place at each position.
    
        Returns
        -------
        list
            Amino acid sequence of the resulting child at the interface positions.
        """
        if positions_to_mutate:
            self.mut.crossover(
                pdb_file=pdb_file,
                output_file=output_file,
                positions_to_mutate=positions_to_mutate,
                aminoacids_to_place=aminoacids_to_place,
            )
        else:
            import shutil
            shutil.copy(pdb_file, output_file)
    
        return self.get_individual_seq(output_file)
    
    
    def crossover_population(self, args):
        """
        Apply crossover to a list of (parent_a_pdb, output_file, positions, aas) tuples in parallel.
    
        Each arg tuple:
            (pdb_file, output_file, positions_to_mutate, aminoacids_to_place)
    
        Parameters
        ----------
        args : list of tuples
            Each tuple: (pdb_file, output_file, positions_to_mutate, aminoacids_to_place)
    
        Returns
        -------
        list of list
            One AA sequence per child, in the same order as args.
        """
        with mp.Pool(processes=mp.cpu_count(), initializer=_init_pyrosetta_worker) as pool:
            results = pool.starmap(self._crossover_worker, args)
        return results


    def get_complete_interest_sequence(self, pdbfile, chain_id):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdbfile)
    
        sequences = {}
    
        for model in structure:
            seq = ""
            for chain in model:
                if chain.id != chain_id:
                    continue
    
                for residue in chain:
                    if residue.id[0] != " ":
                        continue
                    try:
                        seq += seq1(residue.resname)
                    except Exception:
                        seq += "X"
    
            if seq:
                sequences[model.id] = seq

        logging.debug(f"SEQUENCES: \n{sequences}\n")
    
        return sequences[0]



    ## return the pose sequence    
    def sequence(self, pdb_file):
        """Get sequence from PDB file
        
        Parameters
        ----------
        pdb_file: str
            PDB file path
            
        Returns
        -------
        str
            Protein sequence
        """
        # Call the object constructor
        seqTemp = ppyrst.prot_mut_py_rosetta(pdb_file, self.partners)
        # Get the sequence
        return seqTemp.sequence()
