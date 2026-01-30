#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 15/03/2025

@author: Rolando Armas
Yachay Tech University
Phage Therapy Group
"""

import subprocess
import prot_interface.prot_settingsI as sets
from prot_interface.logging_config import setup_logging
import logging

# Initialize logging before anything else
setup_logging()
logger = logging.getLogger(__name__)

''' Class for get the list of stabilized/unstabilized amino-acids'''

class prot_aa_extract:

	def __init__(self, scenario, ligand_chain):



		# #Config parameters
		self.ligand_chain = ligand_chain
		self.start_pattern = sets.START_PATTERN
		self.rosetta_bin = sets.ROSETTA_BIN
		self.config_path = sets.CONFIG_PATH
		self.interf_en = sets.INTERF_EN
		self.face1_file_name = sets.FACE1_FILE_NAME
		self.face2_file_name = sets.FACE2_FILE_NAME
		self.scenario = scenario

	def extract_text(self,file_name):
		extracted_lines = []
		capturing = False  # Flag to indicate when to start capturing

		with open(file_name, 'r', encoding='utf-8') as file:
			for line in file:
				if self.start_pattern in line:  # Detect start pattern
					capturing = True

				if capturing:
					if line.strip() == "":  # Stop at the first blank line
						break
					extracted_lines.append(line.rstrip())  # Append non-empty lines
		return extracted_lines


	###Version 3
	def extract_interactions(self, energy_interact, ligand_chain):
		"""
    	Extracts pairwise interactions with positive and negative energies from a file, focusing on a specific ligand chain.

    	Args:
       	filename (str): Name of the file containing interaction data.
       	ligand_chain (str): The chain identifier for the ligand (e.g., 'C').

    	Returns:
       	tuple: Two lists - positive_interactions and negative_interactions, two sums - sum_positive and sum_negative,
              and a dictionary of ligand entities and their interactions.
    	"""
		sum_positive = 0.0  # Initialize sum for positive energies
		sum_negative = 0.0  # Initialize sum for negative energies
		ligand_interactions = {}  # Dictionary to store interactions for ligand entities

		for line in energy_interact:
			try:
                # Skip lines that don't contain interaction data
				if "---" not in line or ":" not in line:
					continue

				# Split the line into interaction and energy part
				interaction, energy_str = line.strip().split(':')
				energy = float(energy_str)  # Convert energy to float

				# Split the interaction into two entities
				entity1, entity2 = interaction.split(' --- ')

				# Check if either entity belongs to the ligand chain
				if f"_{ligand_chain}(" in entity1:
					ligand_entity = entity1
					other_entity = entity2
				elif f"_{ligand_chain}(" in entity2:
					ligand_entity = entity2
					other_entity = entity1
				else:
					# Skip interactions not involving the ligand chain
					continue

				# Add interaction to ligand entity's list
				if ligand_entity not in ligand_interactions:
					ligand_interactions[ligand_entity] = {"interactions": {}, "total_energy": 0.0, "positive_interactions":0,"negative_interactions":0}
				ligand_interactions[ligand_entity]["interactions"][other_entity] = energy
				ligand_interactions[ligand_entity]["total_energy"] += energy  # Update group sum

				if energy > 0:
					sum_positive += energy  # Add to positive sum
					ligand_interactions[ligand_entity]["positive_interactions"] += 1
				elif energy < 0:
					sum_negative += energy  # Add to negative sum
					ligand_interactions[ligand_entity]["negative_interactions"] += 1

			except ValueError as e:
				logging.error(f"Skipping line: {line.strip()} - Error: {e}")
		
		return sum_positive, sum_negative, ligand_interactions


	def energy_interact_file(self, pdb_file_name):
		##Execute the command and reeturn the file name
		# #Build the command
		output_file_name = pdb_file_name + ".txt"
		command = self.rosetta_bin + self.interf_en + \
          " " + pdb_file_name + \
          " " + "-face1 " + self.config_path + self.scenario + "/" + self.face1_file_name + \
          " " + "-face2 " + self.config_path + self.scenario + "/" + self.face2_file_name + \
          " " + "-score:hbond_bb_per_residue_energy >" + \
          " " + output_file_name


		try:
			return_code = subprocess.call(command, shell=True)
			subprocess.call("sleep 5", shell=True)			

		except subprocess.CalledProcessError as e:
			logging.critical(f"Unexpected error trying to run command: {command}, {return_code}: return_code")
			logging.critical(e.output)
			exit(1)

		energy_file_path = output_file_name

		return energy_file_path
	

	def aa_stab_nstab_list (self, pdb_file_name):
		
		# #Build the command
		
		output_file_name = self.energy_interact_file(pdb_file_name)
		extracted_text = self.extract_text(output_file_name)

		sum_positive, sum_negative, ligand_interactions = self.extract_interactions(extracted_text, self.ligand_chain)


		list_of_aa_s = []
		list_of_aa_ns = []
		for ligand_entity, data in ligand_interactions.items():
			total_energy = data["total_energy"]
			n_positive = data["positive_interactions"]
			n_negative = data["negative_interactions"]
			if total_energy >= 0:
				list_of_aa_ns.append([ligand_entity, total_energy, n_positive, n_negative])
			else:
				list_of_aa_s.append([ligand_entity, total_energy, n_positive, n_negative])

		return list_of_aa_ns, list_of_aa_s
