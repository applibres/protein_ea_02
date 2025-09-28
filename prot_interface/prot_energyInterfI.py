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
from prot_interface.logging_config import setup_logging
import logging

# Initialize logging before anything else
setup_logging()
logger = logging.getLogger(__name__)


''' Class for get the interface energy between protein chains'''

class prot_energyInterf:

	def __init__(self, scenario):

		self.rosetta_bin = sets.ROSETTA_BIN
		self.config_path = sets.CONFIG_PATH
		self.interf_an = sets.INTERF_AN
		self.score_indexes = sets.SCORE_INDEXES
		self.flags = sets.FLAGS
		self.scenario = scenario

	def getEnergyInterf(self, pdb_file_name):

		#Build the command

		# get the name of pdb_file_name
		output_file_name = pdb_file_name + ".sc"
		
		command = self.rosetta_bin + self.interf_an + \
			" " + pdb_file_name + \
			" " + self.flags + \
			" " + output_file_name


		try:
			return_code = subprocess.call(command, shell=True)

		except subprocess.CalledProcessError as e:
			logging.critical(f"Unexpected error trying to run command: {command}, {return_code}: return_code")
			logging.critical(e.output)
			exit(1)

		# open the output file 
		file = open(output_file_name) 

		# read the content of the file opened 
		content = file.readlines() 

		#read the third scores line 
		scores = content[2]

		#We use a regular expression to find all numeric fields in line

		result = re.findall(r'[-+]?\d*\.\d+|[-+]?\d+', scores)
		
		# Convert strings to number
		result_numbers = [float(num_str) for num_str in result]

		energy_intf = []
		for i in self.score_indexes:
			energy_intf.append(result_numbers[i])

		return energy_intf

