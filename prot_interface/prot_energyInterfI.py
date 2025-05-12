#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 22/03/2025

@author: Rolando Armas
Yachay Tech University
Phage Therapy Group
"""

from . import prot_settingsI as sets
import subprocess
import re



''' Class for get the interface energy between protein chains'''

class prot_energyInterf:

	def __init__(self, scenario):

		self.rosetta_bin = sets.ROSETTA_BIN
		self.config_path = sets.CONFIG_PATH
		self.interf_an = sets.INTERF_AN
		self.score_indexes = sets.SCORE_INDEXES
		self.flags = sets.FLAGS
		self.scenario = scenario
		#self.output_path = sets.OUTPUT_PATH

	def getEnergyInterf(self, pdb_file_name):

		#Build the command
		#pattern = '[\w-]+?(?=\.)'

		# get the name of pdb_file_name
		#output_file_name = re.search(pattern, pdb_file_name).group() + ".sc"
		output_file_name = pdb_file_name + ".sc"

		# printing the match
		#print(a.group())

		# command = self.rosetta_bin + self.interf_an + \
		# 	" " + self.config_path + pdb_file_name + \
		# 	" " + self.flags + \
		# 	" " + self.output_path + output_file_name

		# command = self.rosetta_bin + self.interf_an + \
		# 	" " + self.config_path + pdb_file_name + \
		# 	" " + self.flags + \
		# 	" " + output_file_name
		
		command = self.rosetta_bin + self.interf_an + \
			" " + pdb_file_name + \
			" " + self.flags + \
			" " + output_file_name


		try:
			return_code = subprocess.call(command, shell=True)

		except subprocess.CalledProcessError as e:
			print("Unexpected error trying to run command: ", command, "return_code: ", return_code )
			print(e.output)

		# open the output file 
		#file = open(self.output_path + output_file_name) 
		file = open(output_file_name) 

		# read the content of the file opened 
		content = file.readlines() 

		#read the third scores line 
		scores = content[2]

		#We use a regular expression to find all numeric fields in line

		#[-+]?: Matches an optional sign (either - or +).
		#\d*: Matches zero or more digits before the decimal point (to allow for numbers like .5).
		#\.: Matches the decimal point.
		#\d+: Matches one or more digits after the decimal point.
		#|: Acts as an OR operator to allow for matching integers as well.
		#[-+]?\d+: Matches positive or negative integers.

		result = re.findall(r'[-+]?\d*\.\d+|[-+]?\d+', scores)
		
		# Convert strings to number
		result_numbers = [float(num_str) for num_str in result]

		energy_intf = []
		for i in self.score_indexes:
			energy_intf.append(result_numbers[i])

		return energy_intf

