#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 15/03/2025

@author: Rolando Armas
Yachay Tech University
Phage Therapy Group
"""

import subprocess
import re
from pathlib import Path
import prot_interface.prot_settingsI as sets

''' Class for get the list of stabilized/unstabilized amino-acids'''

class prot_aa_extract:

	def __init__(self, scenario, ligand_chain):



		# #Config parameters
		# self.ligand_chain = 'C'  # Specify the ligand chain (e.g., 'C')
		# self.start_pattern = "##### PAIRWISE SHORT-RANGE ENERGIES #####"
		# self.rosetta_bin = "PATH/rosetta.binary.m1.release-371/main/source/bin/"
		# self.config_path = "PATH/tests/"
		# self.interf_en = "interface_energy.static.macosclangrelease -s"
		# self.face1_file_name = "faceA.txt"
		# self.face2_file_name = "faceC.txt"


		self.ligand_chain = ligand_chain
		self.start_pattern = sets.START_PATTERN
		self.rosetta_bin = sets.ROSETTA_BIN
		self.config_path = sets.CONFIG_PATH
		self.interf_en = sets.INTERF_EN
		self.face1_file_name = sets.FACE1_FILE_NAME
		self.face2_file_name = sets.FACE2_FILE_NAME
		self.scenario = scenario
		#self.output_path = sets.OUTPUT_PATH


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
    
		#return "\n".join(extracted_lines)
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
		positive_interactions = []  # Initialize empty list for positive interactions
		negative_interactions = []  # Initialize empty list for negative interactions
		sum_positive = 0.0  # Initialize sum for positive energies
		sum_negative = 0.0  # Initialize sum for negative energies
		ligand_interactions = {}  # Dictionary to store interactions for ligand entities

		for line in energy_interact:
			try:
                # Skip lines that don't contain interaction data
				if "---" not in line or ":" not in line:
					#print(f"Skipping line: {line.strip()} - Invalid format")
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
					#positive_interactions.append((interaction, energy))
					sum_positive += energy  # Add to positive sum
					ligand_interactions[ligand_entity]["positive_interactions"] += 1
				elif energy < 0:
					#negative_interactions.append((interaction, energy))
					sum_negative += energy  # Add to negative sum
					ligand_interactions[ligand_entity]["negative_interactions"] += 1

			except ValueError as e:
				print(f"Skipping line: {line.strip()} - Error: {e}")
		
		return sum_positive, sum_negative, ligand_interactions


	def energy_interact_file (self, pdb_file_name):
		##Execute the command and reeturn the file name

		#Build the command
		#pattern = '[\w-]+?(?=\.)'

		# get the name of pdb_file_name
		#output_file_name = re.search(pattern, pdb_file_name).group() + ".txt"
		#output_file_name = "output04.txt"
		#output_path = "PATH/Research/Bio/tests/"
		output_file_name = pdb_file_name + ".txt"

		# command = self.rosetta_bin + self.interf_en + \
        #   " " + self.config_path + pdb_file_name + \
        #   " " + "-face1 " + self.config_path + self.face1_file_name + \
        #   " " + "-face2 " + self.config_path + self.face2_file_name + \
        #   " " + "-score:hbond_bb_per_residue_energy >" + \
        #   " " + self.output_path + output_file_name

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
			print("Unexpected error trying to run command: ", command, "return_code: ", return_code )
			print(e.output)

		#energy_file_path = self.output_path + output_file_name
		energy_file_path = output_file_name

		return energy_file_path
	

	def aa_stab_nstab_list (self, pdb_file_name):
		
		# #Build the command
		# pattern = '[\w-]+?(?=\.)'

		# # get the name of pdb_file_name
		# output_file_name = re.search(pattern, pdb_file_name).group() + ".txt"
		# #output_file_name = "output04.txt"
		# #output_path = "PATH/tests/"


		# command = self.rosetta_bin + self.interf_en + \
        #   " " + self.config_path + pdb_file_name + \
        #   " " + "-face1 " + self.config_path + self.face1_file_name + \
        #   " " + "-face2 " + self.config_path + self.face2_file_name + \
        #   " " + "-score:hbond_bb_per_residue_energy >" + \
        #   " " + self.output_path + output_file_name

		# try:
		# 	return_code = subprocess.call(command, shell=True)

		# except subprocess.CalledProcessError as e:
		# 	print("Unexpected error trying to run command: ", command, "return_code: ", return_code )
		# 	print(e.output)

  		# #list_aa = self.prot_aa_extract(output_path + output_file_name) 
		
		output_file_name = self.energy_interact_file(pdb_file_name)
		#extracted_text = self.extract_text(self.output_path + output_file_name)
		extracted_text = self.extract_text(output_file_name)

		sum_positive, sum_negative, ligand_interactions = self.extract_interactions(extracted_text, self.ligand_chain)


		list_of_aa_s = []
		list_of_aa_ns = []
		positive_interactions_ranked = []
		negative_interactions_ranked = []
		for ligand_entity, data in ligand_interactions.items():
    		#interactions = data["interactions"]
			total_energy = data["total_energy"]
			n_positive = data["positive_interactions"]
			n_negative = data["negative_interactions"]
    		#print(f"{ligand_entity} (Total Energy: {total_energy:.6f},+I:{n_positive},-I:{n_negative})")
			if total_energy >= 0:
				list_of_aa_ns.append([ligand_entity, total_energy, n_positive, n_negative])
        		# Rank positive interactions (most positive first)
				positive_interactions_ranked = sorted(list_of_aa_ns, key=lambda x: x[1], reverse=True)
			else:
				list_of_aa_s.append([ligand_entity, total_energy, n_positive, n_negative])
        		# Rank negative interactions (most negative first)
				negative_interactions_ranked = sorted(list_of_aa_s, key=lambda x: x[1])

		return positive_interactions_ranked, negative_interactions_ranked


	# def extract_interactions(self, filename, ligand_chain):
	# 	"""
    # 	Extracts pairwise interactions with positive and negative energies from a file, focusing on a specific ligand chain.

    # 	Args:
    #    	filename (str): Name of the file containing interaction data.
    #    	ligand_chain (str): The chain identifier for the ligand (e.g., 'C').

    # 	Returns:
    #    	tuple: Two lists - positive_interactions and negative_interactions, two sums - sum_positive and sum_negative,
    #           and a dictionary of ligand entities and their interactions.
    # 	"""
	# 	positive_interactions = []  # Initialize empty list for positive interactions
	# 	negative_interactions = []  # Initialize empty list for negative interactions
	# 	sum_positive = 0.0  # Initialize sum for positive energies
	# 	sum_negative = 0.0  # Initialize sum for negative energies
	# 	ligand_interactions = {}  # Dictionary to store interactions for ligand entities

	# 	try:
	# 		with open(filename, 'r') as f:
	# 			for line in f:
	# 				try:
    #                 	# Skip lines that don't contain interaction data
	# 					if "---" not in line or ":" not in line:
	# 						print(f"Skipping line: {line.strip()} - Invalid format")
	# 						continue

	# 					# Split the line into interaction and energy part
	# 					interaction, energy_str = line.strip().split(':')
	# 					energy = float(energy_str)  # Convert energy to float

	# 					# Split the interaction into two entities
	# 					entity1, entity2 = interaction.split(' --- ')

	# 					# Check if either entity belongs to the ligand chain
	# 					if f"_{ligand_chain}(" in entity1:
	# 						ligand_entity = entity1
	# 						other_entity = entity2
	# 					elif f"_{ligand_chain}(" in entity2:
	# 						ligand_entity = entity2
	# 						other_entity = entity1
	# 					else:
	# 						# Skip interactions not involving the ligand chain
	# 						continue

	# 					# Add interaction to ligand entity's list
	# 					if ligand_entity not in ligand_interactions:
	# 						ligand_interactions[ligand_entity] = {"interactions": {}, "total_energy": 0.0}
	# 					ligand_interactions[ligand_entity]["interactions"][other_entity] = energy
	# 					ligand_interactions[ligand_entity]["total_energy"] += energy  # Update group sum

	# 					if energy > 0:
	# 						positive_interactions.append((interaction, energy))
	# 						sum_positive += energy  # Add to positive sum
	# 					elif energy < 0:
	# 						negative_interactions.append((interaction, energy))
	# 						sum_negative += energy  # Add to negative sum
	# 				except ValueError as e:
	# 					print(f"Skipping line: {line.strip()} - Error: {e}")
	# 	except FileNotFoundError:
	# 		print(f"Error: The file '{filename}' was not found.")
	# 	except Exception as e:
	# 		print(f"An error occurred: {e}")

	# 	return positive_interactions, negative_interactions, sum_positive, sum_negative, ligand_interactions

	### Version 2
	# def extract_interactions(self, energy_interact, ligand_chain):
	# 	"""
    # 	Extracts pairwise interactions with positive and negative energies from a file, focusing on a specific ligand chain.

    # 	Args:
    #    	filename (str): Name of the file containing interaction data.
    #    	ligand_chain (str): The chain identifier for the ligand (e.g., 'C').

    # 	Returns:
    #    	tuple: Two lists - positive_interactions and negative_interactions, two sums - sum_positive and sum_negative,
    #           and a dictionary of ligand entities and their interactions.
    # 	"""
	# 	positive_interactions = []  # Initialize empty list for positive interactions
	# 	negative_interactions = []  # Initialize empty list for negative interactions
	# 	sum_positive = 0.0  # Initialize sum for positive energies
	# 	sum_negative = 0.0  # Initialize sum for negative energies
	# 	ligand_interactions = {}  # Dictionary to store interactions for ligand entities

	# 	for line in energy_interact:
	# 		try:
    #             # Skip lines that don't contain interaction data
	# 			if "---" not in line or ":" not in line:
	# 				#print(f"Skipping line: {line.strip()} - Invalid format")
	# 				continue

	# 			# Split the line into interaction and energy part
	# 			interaction, energy_str = line.strip().split(':')
	# 			energy = float(energy_str)  # Convert energy to float

	# 			# Split the interaction into two entities
	# 			entity1, entity2 = interaction.split(' --- ')

	# 			# Check if either entity belongs to the ligand chain
	# 			if f"_{ligand_chain}(" in entity1:
	# 				ligand_entity = entity1
	# 				other_entity = entity2
	# 			elif f"_{ligand_chain}(" in entity2:
	# 				ligand_entity = entity2
	# 				other_entity = entity1
	# 			else:
	# 				# Skip interactions not involving the ligand chain
	# 				continue

	# 			# Add interaction to ligand entity's list
	# 			if ligand_entity not in ligand_interactions:
	# 				ligand_interactions[ligand_entity] = {"interactions": {}, "total_energy": 0.0, "positive_interactions":0,"negative_interactions":0}
	# 			ligand_interactions[ligand_entity]["interactions"][other_entity] = energy
	# 			ligand_interactions[ligand_entity]["total_energy"] += energy  # Update group sum

	# 			if energy > 0:
	# 				positive_interactions.append((interaction, energy))
	# 				sum_positive += energy  # Add to positive sum
	# 				ligand_interactions[ligand_entity]["positive_interactions"] += 1
	# 			elif energy < 0:
	# 				negative_interactions.append((interaction, energy))
	# 				sum_negative += energy  # Add to negative sum
	# 				ligand_interactions[ligand_entity]["negative_interactions"] += 1

	# 		except ValueError as e:
	# 			print(f"Skipping line: {line.strip()} - Error: {e}")
		
	# 	return positive_interactions, negative_interactions, sum_positive, sum_negative, ligand_interactions
