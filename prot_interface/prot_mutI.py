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

        self.esm2_prob_matrix = ESM2ProbMatrix()
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

    def plan_mutation(self, scenario, ligand_chain, pdb_file, output_file, mut_rate,
                      sequence, generation, ngen, original_sequence, parent_id=None):
        """Plan mutation steps without touching PDB files."""
        list_aa = prot_aa.prot_aa_extract(scenario, ligand_chain)
        aans, aas = list_aa.aa_stab_nstab_list(pdb_file)
        logging.debug(f"AANS aminoacids: {aans}")
        logging.debug(f"AAS aminoacids: {aas}")

        mutations = []

        if len(aans) > 0 or len(aas) > 0:
            aminoacids = aans + aas
            energies = np.array([energ[1] for energ in aminoacids])
            max_e = np.max(energies)

            weights_aas = np.exp(energies - max_e) / np.sum(np.exp(energies - max_e))

            logging.debug(f"Sequence for ESM2: {sequence}")

            max_mut = len(aminoacids) * mut_rate
            logging.debug(f"Max mut: {max_mut}")

            if int(max_mut) <= 1:
                num_of_mut = 1
            else:
                num_of_mut = np.random.randint(1, int(max_mut))

            logging.info("Number of Mutations =%s", num_of_mut)

            for _ in range(num_of_mut):
                aa2mut = random.choices(aminoacids, weights=weights_aas, k=1)[0]

                logging.info("AA to Mutate: %s", aa2mut)

                aa = aa2mut[0][0:3]
                logging.info("Amino to replace: %s %s %s", aa, "-", self.wildtype(aa))

                aa_pos = re.findall(r'\d+', aa2mut[0])
                logging.info("In Position: %s", aa_pos[1])

                position_in_seq = self.positions.index(int(aa_pos[0]))
                posi = int(aa_pos[1])

                aa_mut = self.esm2_prob_matrix.most_probable_replacement(
                    sequence, position_in_seq, generation, ngen
                )
                logging.info("Decision: %s %s %s", aa2mut, " --> ", aa_mut[0])
                res = self.aatype(aa_mut[0])

                mutations.append((posi, res, int(aa_pos[0])))

        mut_sequence = list(sequence[:])
        for _, res, posi in mutations:
            mut_sequence[posi] = self.wildtype(res)

        interface_sequence = [mut_sequence[a] for a in self.positions]
        mutated_interface_sequence = "".join(interface_sequence)
        logging.info(f"Sequence mutated: {mut_sequence}")

        return {
            "parent_id": parent_id,
            "scenario": scenario,
            "ligand_chain": ligand_chain,
            "pdb_in": pdb_file,
            "pdb_out": output_file,
            "generation": generation,
            "ngen": ngen,
            "mut_rate": mut_rate,
            "mutations": mutations,
            "mutated_interface_sequence": mutated_interface_sequence,
        }

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
