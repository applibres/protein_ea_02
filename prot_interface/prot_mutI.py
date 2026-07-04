#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 18/03/2025

@author: Rolando Armas
Yachay Tech University
Phage Therapy Group
"""

''' Class Mutation Operator'''

##Pyrosetta Libs
import os
import shutil
import subprocess
import pyrosetta
from pyrosetta.rosetta.core.pack.task import *
from pyrosetta.toolbox import *

# Implemented Libs
import prot_interface.prot_aa_stI as prot_aa

## Python Libs
import re
import random
import pandas as pd
import numpy as np

import prot_interface.prot_settingsI as sets
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.pack.task import TaskFactory, operation
from pyrosetta.rosetta.core.select import residue_selector as rs
from pyrosetta.rosetta.core.select.movemap import MoveMapFactory, mm_enable
from prot_interface.logging_config import setup_logging
from prot_interface.prot_esm2 import ESM2ProbMatrix
import logging

# Initialize logging before anything else
setup_logging()
logger = logging.getLogger(__name__)

class prot_mut:

    def __init__(self,scenario):
        
    	#Config Files

        self.scenario = scenario
        self.config_path = sets.CONFIG_PATH + self.scenario + "/"
        self.scenarios_path = sets.CONFIG_PATH
        self.matrix_file_name = sets.MSA_MATRIX
        self.rosetta_bin = sets.ROSETTA_BIN

        # Lazy initialization to avoid pickling heavy torch storages when
        # multiprocessing serializes bound methods.
        self.esm2_prob_matrix = None
        self.positions = self.getPositions()

        """
        Load the MSA matrix from a TSV (tab-separated values) file.
        Assumes the first column contains amino acids, and the rest are probabilities.
        """
        #self.msa_df = pd.read_csv(self.config_path + self.matrix_file_name, sep='\t', index_col=0)
    
		#Initialize pyrosetta
        
        
        pyrosetta.init(
            "-nstruct 1 "
            "-ignore_zero_occupancy false "
            "-ex1 -ex2 "
            "-use_input_sc "
            "-flip_HNQ "
            "-no_optH false"
        )
        #pyrosetta.init()


    def getPositions(self):
        with open(self.config_path + sets.FACE2_FILE_NAME, 'r') as file:
            faceE_data = file.read()
        positions = list(map(lambda x: re.findall(r'\d+', x), faceE_data.split(" _")))
        positions = [int(num) for sublist in positions for num in sublist]
        return positions

    def _get_esm2_prob_matrix(self):
        if self.esm2_prob_matrix is None:
            self.esm2_prob_matrix = ESM2ProbMatrix()
        return self.esm2_prob_matrix

    def get_probabilities(self, amino_acid, threshold=0.0):
        """
        Given an amino acid, return the probability row from the MSA matrix.
        """
        if amino_acid not in self.msa_df.index:
            raise ValueError(f"Amino acid {amino_acid} not found in MSA matrix.")
        
        row = self.msa_df.loc[amino_acid]

        return [[col, float(row[col])] for col in self.msa_df.columns if row[col] >= threshold]


    def wildtype(self,aatype):
        """This function returns the wild type amino acids"""

        AA = ['G','A','L','M','F','W','K','Q','E','S','P'
            ,'V','I','C','Y','H','R','N','D','T']

        AA_3 = ['GLY','ALA','LEU','MET','PHE','TRP'
            ,'LYS','GLN','GLU', 'SER','PRO','VAL'
            ,'ILE','CYS','TYR','HIS','ARG','ASN'
            ,'ASP','THR']

        for i in range(0, len(AA_3)):
            if(aatype == AA_3[i]):
                return AA[i]
            
    def aatype(self,wildtype):
        """This function returns the aa type amino acids"""

        AA = ['G','A','L','M','F','W','K','Q','E','S','P'
            ,'V','I','C','Y','H','R','N','D','T']

        AA_3 = ['GLY','ALA','LEU','MET','PHE','TRP'
            ,'LYS','GLN','GLU', 'SER','PRO','VAL'
            ,'ILE','CYS','TYR','HIS','ARG','ASN'
            ,'ASP','THR']

        for i in range(0, len(AA)):
            if(wildtype == AA[i]):
                return AA_3[i]
    
    def mutate_local_relax(self, input_pdb, output_pdb, pos, aa):
        """
        Ejecuta Rosetta Scripts y renombra la salida al nombre deseado.
        
        Args:
            input_pdb (str): Ruta al PDB de entrada (ej: /ruta/input.pdb)
            output_pdb (str): Ruta y nombre exacto deseado (ej: /ruta/final_mutado.pdb)
            pos (int): Posición a mutar.
            aa (str): Aminoácido destino (3 letras).
            i (int/str): Identificador único para el sufijo temporal.
        """
        input_filename = os.path.basename(input_pdb)
        input_stem = os.path.splitext(input_filename)[0]
        output_dir = os.path.dirname(output_pdb)
        
        temp_suffix = f"_{random.random()}"
        temp_filename = f"{input_stem}{temp_suffix}.pdb"
        temp_path = os.path.join(output_dir, temp_filename)
        
        command = (
            f"{self.rosetta_bin}rosetta_scripts.default.linuxgccrelease "
            f"-database /usr/local/database "
            f"-s {input_pdb} "
            f"-parser:protocol {self.scenarios_path}/PM_Mutation_Relax_Local.xml "
            f"-parser:script_vars pos1={pos} aa={aa} "
            f"-ignore_unrecognized_res "
            f"-out:path:all {output_dir} "
            f"-out:suffix {temp_suffix} "
            f"-out:no_nstruct_label true "
            f"&& mv {temp_path} {output_pdb}"
        )

        try:
            subprocess.run(command, shell=True, check=True, executable='/bin/bash')
        except subprocess.CalledProcessError as e:
            print(f'Error: {e}')

    def _interface_from_complete(self, complete_sequence: str) -> str:
        # faceB residue numbers are 1-based for the ligand chain.
        # Keep ascending residue order to match create_individual0/get_individual_seq.
        interface = []
        for residue_number in sorted(self.positions):
            idx = int(residue_number) - 1
            if 0 <= idx < len(complete_sequence):
                interface.append(complete_sequence[idx])
        return "".join(interface)

    def _empty_plan(self, scenario, ligand_chain, pdb_file, output_file,
                    generation, ngen, parent_id=None, sequence=""):
        interface_sequence = self._interface_from_complete(sequence)
        return {
            "parent_id": parent_id,
            "scenario": scenario,
            "ligand_chain": ligand_chain,
            "pdb_in": pdb_file,
            "pdb_out": output_file,
            "generation": generation,
            "ngen": ngen,
            "mutations": [],
            "mutated_interface_sequence": interface_sequence,
            "mutated_complete_sequence": sequence,
            "position_rank": -1,
            "aa_rank": -1,
            "candidate_score": 0.0,
        }

    def plan_mutation_candidates(self, scenario, ligand_chain, pdb_file, output_file,
                                 sequence, generation, ngen, original_sequence, parent_id=None,
                                 top_positions: int = 8, top_aa: int = 8):
        """Build deterministic 8x8 single-point mutation candidates for one parent slot."""
        list_aa = prot_aa.prot_aa_extract(scenario, ligand_chain)
        aans, aas = list_aa.aa_stab_nstab_list(pdb_file)
        aminoacids = aans + aas

        logging.debug(f"AANS aminoacids: {aans}")
        logging.debug(f"AAS aminoacids: {aas}")
        logging.debug(f"Sequence for ESM2: {sequence}")

        if not aminoacids:
            return [
                self._empty_plan(
                    scenario, ligand_chain, pdb_file, output_file,
                    generation, ngen, parent_id, sequence=sequence
                )
            ]

        top_positions = max(1, int(top_positions))
        top_aa = max(1, int(top_aa))

        # Higher interaction-energy entries have higher priority in ranking.
        ranked_positions = sorted(aminoacids, key=lambda row: float(row[1]), reverse=True)[:top_positions]
        candidates = []

        for pos_rank, aa2mut in enumerate(ranked_positions):
            aa_pos = re.findall(r'\d+', aa2mut[0])
            if len(aa_pos) < 2:
                continue

            residue_number = int(aa_pos[0])
            rosetta_position = int(aa_pos[1])
            sequence_index = residue_number - 1

            if sequence_index < 0 or sequence_index >= len(sequence):
                continue

            aa_options = self._get_esm2_prob_matrix().top_k_replacements(
                seq=sequence,
                position=sequence_index,
                generation=generation,
                ngen=ngen,
                k=top_aa,
                exclude_wt=True,
            )

            for aa_rank, (aa_mut, _) in enumerate(aa_options):
                res = self.aatype(aa_mut)
                if res is None:
                    continue

                mutations = [(rosetta_position, res, residue_number)]

                mut_complete_sequence = list(sequence)
                mut_complete_sequence[sequence_index] = aa_mut
                mut_complete_sequence = "".join(mut_complete_sequence)
                mutated_interface_sequence = self._interface_from_complete(mut_complete_sequence)

                candidates.append(
                    {
                        "parent_id": parent_id,
                        "scenario": scenario,
                        "ligand_chain": ligand_chain,
                        "pdb_in": pdb_file,
                        "pdb_out": output_file,
                        "generation": generation,
                        "ngen": ngen,
                        "mutations": mutations,
                        "mutated_interface_sequence": mutated_interface_sequence,
                        "mutated_complete_sequence": mut_complete_sequence,
                        "position_rank": int(pos_rank),
                        "aa_rank": int(aa_rank),
                        "candidate_score": float(-(pos_rank + aa_rank)),
                    }
                )

        if not candidates:
            return [
                self._empty_plan(
                    scenario, ligand_chain, pdb_file, output_file,
                    generation, ngen, parent_id, sequence=sequence
                )
            ]

        return candidates

    def plan_mutation(self, scenario, ligand_chain, pdb_file, output_file,
                      sequence, generation, ngen, original_sequence, parent_id=None):
        """
        Backwards-compatible wrapper: returns one candidate sampled from the
        deterministic candidate list.
        """
        candidates = self.plan_mutation_candidates(
            scenario=scenario,
            ligand_chain=ligand_chain,
            pdb_file=pdb_file,
            output_file=output_file,
            sequence=sequence,
            generation=generation,
            ngen=ngen,
            original_sequence=original_sequence,
            parent_id=parent_id,
        )
        if not candidates:
            return self._empty_plan(
                scenario, ligand_chain, pdb_file, output_file, generation, ngen, parent_id, sequence=sequence
            )
        return random.choice(candidates)

    def apply_mutation_plan(self, plan):
        """Apply a precomputed mutation plan by executing mutate_local_relax."""
        pdb_file = plan["pdb_in"]
        output_file = plan["pdb_out"]
        mutations = plan["mutations"]

        path_pdb = os.path.dirname(os.path.abspath(output_file))
        base = os.path.basename(output_file)
        temp = os.path.join(path_pdb, f"tmp_{str(base).split('.pdb')[0]}")
        os.makedirs(temp, exist_ok=True)

        current_pdb = pdb_file

        for i, (posi, res, _) in enumerate(mutations):
            temp_file = os.path.join(temp, f"{base.split('.pdb')[0]}_{i}.pdb")
            self.mutate_local_relax(current_pdb, temp_file, posi, res)
            current_pdb = temp_file

            if os.path.isfile(current_pdb):
                logging.info(f"Temporal mutant: {current_pdb}")
            else:
                logging.info(f"Not found: {current_pdb}")

        if os.path.abspath(current_pdb) == os.path.abspath(output_file):
            pass
        elif os.path.abspath(current_pdb) == os.path.abspath(pdb_file):
            logging.info(f"Copying Mutant pdb: {current_pdb} -> {output_file}")
            shutil.copy2(current_pdb, output_file)
        else:
            logging.info(f"Moving Mutant pdb: {current_pdb} -> {output_file}")
            shutil.move(current_pdb, output_file)

        if os.path.isdir(temp):
            shutil.rmtree(temp)

        return {
            "parent_id": plan.get("parent_id"),
            "pdb": output_file,
            "mutated_interface_sequence": plan.get("mutated_interface_sequence", ""),
            "mutations": mutations,
        }
    
    def crossover(self, pdb_file, output_file, positions_to_mutate, aminoacids_to_place):
        """
        Perform crossover on a protein structure by applying a predefined list
        of mutations (position → amino acid), without probabilistic selection.
    
        Parameters
        ----------
        pdb_file : str
            Input PDB file path.
        output_file : str
            Output PDB file path.
        positions_to_mutate : list[int]
            List of residue positions (Rosetta numbering) to mutate.
            e.g. [45, 78, 112]
        aminoacids_to_place : list[str]
            List of single-letter amino acid codes to place at each position.
            Must have the same length as positions_to_mutate.
            e.g. ['A', 'V', 'L']
    
        Raises
        ------
        ValueError
            If positions_to_mutate and aminoacids_to_place have different lengths.
    
        Example
        -------
        prot.crossover(
            pdb_file="input.pdb",
            output_file="crossover.pdb",
            positions_to_mutate=[45, 78, 112],
            aminoacids_to_place=['A', 'V', 'L']
        )
        """
        if len(positions_to_mutate) != len(aminoacids_to_place):
            raise ValueError(
                f"positions_to_mutate (len={len(positions_to_mutate)}) and "
                f"aminoacids_to_place (len={len(aminoacids_to_place)}) must have the same length."
            )
    
        path_pdb = os.path.dirname(os.path.abspath(output_file))
        base = os.path.basename(output_file)
    
        temp = os.path.join(path_pdb, f"tmp_crossover_{str(base).split('.pdb')[0]}")
        os.makedirs(temp, exist_ok=True)
    
        current_pdb = pdb_file
    
        for i, (posi, aa_single) in enumerate(zip(positions_to_mutate, aminoacids_to_place)):
    
            res = self.aatype(aa_single)   # single-letter → 3-letter (e.g. 'A' → 'ALA')
    
            if res is None:
                logging.warning(f"[Crossover] Unknown amino acid '{aa_single}' at position {posi}. Skipping.")
                continue
    
            logging.info(f"[Crossover] Step {i+1}/{len(positions_to_mutate)}: "
                         f"pos={posi}  aa={aa_single} ({res})")
    
            temp_file = os.path.join(temp, f"{base.split('.pdb')[0]}_xover_{i}.pdb")
            self.mutate_local_relax(current_pdb, temp_file, posi, res)
    
            if os.path.isfile(temp_file):
                logging.info(f"[Crossover] Temporal crossover file: {temp_file}")
                current_pdb = temp_file
            else:
                logging.warning(f"[Crossover] File not found after step {i+1}: {temp_file}. "
                                "Keeping previous structure.")
    
        logging.info(f"[Crossover] Moving final structure: {current_pdb} -> {output_file}")
        shutil.move(current_pdb, output_file)
        shutil.rmtree(temp)
