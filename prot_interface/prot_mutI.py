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
import pyrosetta
from pyrosetta.rosetta.core.pack.task import *
from pyrosetta.toolbox import *

# Implemented Libs
import prot_interface.prot_aa_stI as prot_aa

## Python Libs
import re
import time
import random
import pandas as pd
import numpy as np

import prot_interface.prot_settingsI as sets
from prot_interface.logging_config import setup_logging
import logging

# Initialize logging before anything else
setup_logging()
logger = logging.getLogger(__name__)

class prot_mut:

    def __init__(self,scenario):
        
    	#Config Files

        self.scenario = scenario
        self.config_path = sets.CONFIG_PATH + self.scenario + "/"
        self.matrix_file_name = sets.MSA_MATRIX

        """
        Load the MSA matrix from a TSV (tab-separated values) file.
        Assumes the first column contains amino acids, and the rest are probabilities.
        """
        self.msa_df = pd.read_csv(self.config_path + self.matrix_file_name, sep='\t', index_col=0)
    
		#Initialize pyrosetta
        pyrosetta.init()


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

    def mutate(self, scenario, ligand_chain ,pdb_file, output_file, mut_rate):

		#Instantiate Objects

		## Step 1: Choose Position to Mutate ####
		# Get list of stable and unstable aa
        list_aa = prot_aa.prot_aa_extract(scenario, ligand_chain)
        aans, aas = list_aa.aa_stab_nstab_list(pdb_file) 


        pyrosetta.init()
        init_pose = pyrosetta.io.pose_from_pdb(pdb_file)
        mut_pose = pyrosetta.io.Pose()
        mut_pose.assign(init_pose)



        if (len(aans) > 0 or len(aas) > 0):
            min_mut = int(1)
            max_mut = 0

            if (len(aans) > 0):
                max_mut = len(aans)*mut_rate
            elif (len(aas) > 0):
                max_mut = len(aas)*mut_rate
            
            if (int(max_mut) <= 1):
                num_of_mut = 1
            else:    
                num_of_mut = np.random.randint(min_mut,int(max_mut))

            logging.info("Number of Mutations =%s", num_of_mut)


            mut_locations = []
            for i in range(num_of_mut):
                if (len(aans)>0):
                    aa2mut = random.choice(aans)
                elif (len(aas)>0):
                    aa2mut = random.choice(aas)
                
                logging.info("AA to Mutate: %s", aa2mut)

                #Get AA
                aa = aa2mut[0][0:3]
                logging.info("Amino to replace: %s %s %s", aa,"-", self.wildtype(aa))

                #Get position to mutate
                aa_pos = re.findall(r'\d+', aa2mut[0])
                logging.info("In Position: %s", aa_pos[1])

                ## Step 2: Choose replace residue from MSA ###
                aa_prob = self.get_probabilities(self.wildtype(aa), 0.2)  # Get probabilities for Alanine
                logging.info("List of Substitute Aminoacids:")
                logging.info(aa_prob)

                aa_mut = random.choice(aa_prob)
                logging.info("Decision: %s %s %s", aa2mut," --> ",aa_mut)
                mut_locations.append([aa_pos[1],aa_mut[0]])

                #Step 3: mutate 
                res = aa_mut[0]
                posi = int(aa_pos[1])
                pyrosetta.toolbox.mutants.mutate_residue(mut_pose,posi,res)

        #Create the mutant in pdb file
        mut_pose.dump_pdb(output_file)




