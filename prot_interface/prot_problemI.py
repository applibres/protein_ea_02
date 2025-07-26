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
import re
import pyrosetta
from pyrosetta.rosetta.protocols.relax import FastRelax
from multiprocessing import Pool
import prot_interface.prot_pyrosettaI as ppyrst
import prot_interface.prot_energyInterfI as prot_en_intf
import prot_interface.prot_mutI as p_mut
import prot_interface.prot_pyrosettaI as ppyrst
import prot_interface.prot_aa_stI as prot_aa
import prot_interface.prot_settingsI as sets
from prot_interface.logging_config import setup_logging
import logging
import os
import shutil

# Initialize logging before anything else
setup_logging()
logger = logging.getLogger(__name__)

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
                match = re.search(r'([A-Z]{3})(\d+)_C\((\d+)\)', line)
                if match:
                    number1 = int(match.group(2))  # Extract Number1
                    number2 = int(match.group(3))  # Extract Number2
                    mapping[number1] = number2  # Store in dictionary

        return dict(sorted(mapping.items()))

    def absolute_AA_positions(AA_map_positions):
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
        aa_pos_dict=self.extract_mappings(energy_filepath)
        logging.debug(aa_pos_dict)

        # Extract absolute position values as a list
        self.aa_pos_list = list(aa_pos_dict.values())
        logging.debug(self.aa_pos_list)


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
        f = [0.0,0.0,0.0,0.0]
        
        #Evaluate the solution
        #f[0]: dG_separated
        #f[1]: dSASA_int
        #f[2]: dG_separated/dSASAx100
        #f[3]: hbonds_int]
        f = self.protEn.getEnergyInterf(pdb_file)
        return f

    ## Evaluate Population in parallel
    def fitnessPop(self,pop_pdb_files):
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
 
        #Create pool to run in parallel
        pool = Pool()
        evaluation=[]
           
            
        #calling pool for run in parallel
            #Fitness:
            #result_f[0]:dG_separated
            #result_f[1]:dSASA_int
            #result_f[2]:dG_separated/dSASAx100
            #result_f[3]:hbonds_int
        

        for result_f in pool.map(self.fitness,pop_pdb_files):
            #return the fitness solution
            evaluation.append(result_f)
        
        return evaluation 


    ##Mutate Population in parallel
    def mutate_population(self,args):
        #Create pool to run in parallel
        pool = Pool()
        population=[]

        for result_f in pool.starmap(self.mutate, args):
            ind = result_f
            population.append(ind)
        
        return population


    ##Mutation Operator##
    def mutate(self,pdb_file,output_file,mut_rate):
        #Mutate
        self.mut.mutate(self.scenario, self.ligand_chain, pdb_file, output_file, mut_rate)
        aa = self.get_individual_seq(output_file)
        return aa



     ##return the pose sequence    
    def sequence(self,pdb_file):
       #Call the object constructor
       
       seqTemp = ppyrst.prot_mut_py_rosetta(pdb_file,self.partners)
       #Get the sequence
       return seqTemp.sequence()
    
    def relax_population(self, pdb_files):
        with Pool() as pool:
            pool.map(self.relax, pdb_files)
        return list(map(lambda file: file.replace(".pdb", "_relaxed.pdb"), pdb_files))
    
    def relax(self, pdb_file):
        #pyrosetta.init() <-- For Mac
        logging.debug(f"Relaxing: {pdb_file}")
        scorefxn = pyrosetta.get_fa_scorefxn()
        relax = FastRelax()
        relax.set_scorefxn(scorefxn)
        relax.constrain_relax_to_start_coords(False)
    
        pose = pyrosetta.pose_from_pdb(pdb_file)
        relax.apply(pose)
        os.remove(pdb_file)
        output_file = pdb_file.replace(".pdb", "_relaxed.pdb")
        pose.dump_pdb(output_file)
