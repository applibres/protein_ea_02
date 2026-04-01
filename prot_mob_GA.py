#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 16/04/2025

@author: Rolando Armas
Yachay Tech University
Phage Therapy Group


Simple Genetic Algorithm for Protein Mutation
"""

from deap import base, creator, tools
from deap.benchmarks.tools import hypervolume
from deap.tools.emo import sortNondominated

from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

from prot_interface.prot_problemI import *
import prot_interface.prot_problemI as problem
import prot_interface.prot_settingsI as sets
import pyrosetta
import os
import shutil
import numpy as np
import csv
import random
from prot_interface.logging_config import setup_logging
from prot_interface.prot_esm2 import ESM2ProbMatrix
from utils import csvToTree, read_scfiles, save_population_to_csv, save_HallofFame, hamming_distance
import logging
import pickle

# Initialize logging before anything else
setup_logging()
logger = logging.getLogger(__name__)

class CustomIndividual(list):
    counter = 0
    def __init__(self, *args):
        super().__init__(*args)
        self.id = f'{CustomIndividual.counter}'
        CustomIndividual.counter += 1
        self.father = None
        self.fitness = creator.FitnessMin()
        self.pdb = None
        self.nmut = 0


class deap_mob_sga_protein:

    def __init__(self, scenario, algoritm_params, sim_params, fitness_idxs, fitness_weights, output, randomseed):
        """Constructor
        
        Parameters
        ----------     
        algorithm_params:string
            algorithm parameters

        pdbfile: string
            pdbfile path and name

        partners: string
            pdbfile chains

        scenario: string
            scenario name

        """
        self.algoritm_params = algoritm_params 
        self.sim_params = sim_params
        self.fitness_idxs = fitness_idxs
        self.fitness_weights = fitness_weights
        self.output = output
        self.randomseed = randomseed
        self.scenario = scenario
        
        #Simulation Params
        self.partners = self.sim_params['partners']
        self.ligand_chain = self.sim_params['ligand_chain']
        self.pdbfile = self.sim_params['pdbfile']
        
        #Initialize pyrosetta FIRST (before problem instantiation)
        pyrosetta.init()
        
        #Instantiate problem
        self.my_protein_problem = problem.prot_problem(self.scenario,self.partners,self.ligand_chain)

        # Initialize ESM2 model
        logging.info("Initializing ESM2 model...")
        self.esm2 = ESM2ProbMatrix()
        logging.info("ESM2 model initialized successfully")


    def unique_offspring(self, population, selection_func, k):
        """
        Select k unique individuals using the given selection function, 
        returning both selected individuals and their original indices.

        :param population: List of individuals.
        :param selection_func: Selection function like toolbox.select.
        :param k: Number of individuals to select.
        :return: List of selected individuals (no duplicates), and their indices.
        """
        selected = set()
        unique_inds = []
        indices = []

        attempts = 0
        max_attempts = k * 10  # Prevent infinite loops

        while len(unique_inds) < k and attempts < max_attempts:
            candidate_inds = selection_func(population, 1)
            for sel in candidate_inds:
                idx = next(i for i, orig in enumerate(population) if sel is orig)
                if idx not in selected:
                    selected.add(idx)
                    unique_inds.append(sel)
                    indices.append(idx)
                    break
            attempts += 1

        if len(unique_inds) < k:
            logging.warning("Could not find enough unique individuals, filled with randoms.")
            while len(unique_inds) < k:
                rand_idx = random.randint(0, len(population) - 1)
                if rand_idx not in selected:
                    unique_inds.append(population[rand_idx])
                    indices.append(rand_idx)
                    selected.add(rand_idx)

        return unique_inds, indices

    def get_complete_sequence(self, pdbfile, chain_id):
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
    
        return sequences[0]


    # Custom simple evolutionary algorithm
    def run(self, checkpoint=False, freq=2):

        """Run the simple genetic algorithm
        Parameters
        ----------
        checkpoint: bool
            If True, the algorithm will save the state of the algorithm in a checkpoint file
        freq: int
            Frequency to save the checkpoint file, default is 2 generations
        """

        #Capture parameters

        #Algorithm Params

        PATH_STATISTICS = os.path.join(os.path.dirname(self.output), "statistics")
        SCFILE_CSV = os.path.join(PATH_STATISTICS, f"run{self.randomseed}_scfile.csv")
        SAVEPOPGEN_CSV = os.path.join(PATH_STATISTICS, f'run{self.randomseed}_individuals.csv')
        SAVEHALLOFFAME_CSV = os.path.join(PATH_STATISTICS, f'run{self.randomseed}_hallofame.csv')
        SAVEHYPERVOLUME_CSV = os.path.join(PATH_STATISTICS, f"run{self.randomseed}_hv.csv")
        SAVE_STATISTICS_CSV = os.path.join(PATH_STATISTICS, f"run{self.randomseed}_statistics.csv")

        os.makedirs(PATH_STATISTICS, exist_ok=True)
        
        popsize = self.algoritm_params['popsize']
        ngenerations = self.algoritm_params['gen']
        nobj = len(self.fitness_idxs)+1
        mutprob = self.algoritm_params['mutp']
        logging.info("###Algorithm Parameters###")
        logging.info(f"popsize = {popsize}")
        logging.info(f"ngenerations = {ngenerations}")
        logging.info(f"nobj = {nobj}")
        logging.info(f"mut prob = {mutprob}")
    
        
        #Create tmp directory
        output_path = self.output+"/tmp/"
        logging.debug(f"output_path = {output_path}")
        

        #Create directory if not exist
        if os.path.isdir(output_path):
            logging.debug(f"{output_path} directory already exist")
        else:
            os.makedirs(output_path)
            logging.debug(f"{output_path} directory created")
    


        output_path = self.output+"/g0/"
        logging.debug(f"output_path = {output_path}")
        

        #Create directory if not exist
        if os.path.isdir(output_path):
            logging.debug(f"{output_path} directory already exist")
        else:
            os.makedirs(output_path)
            logging.debug(f"{output_path} directory created")


        ##Declare FitnessMinimization and Individual
        creator.create("FitnessMin", base.Fitness, weights=(*self.fitness_weights, 1.0))  # Minimization problem
        creator.create("Individual", CustomIndividual)

        # Initialize toolbox
        toolbox = base.Toolbox()
        
        ##Selection operator implemented by framework##
        toolbox.register("select", tools.selNSGA2)
        
        # Statistics
        stats = tools.Statistics(key=lambda ind: ind.fitness.values)
        stats.register("avg", np.mean, axis=0)
        stats.register("std", np.std, axis=0)
        stats.register("min", np.min, axis=0)
        stats.register("max", np.max, axis=0)

        #Logbook to save statistics
        logbook = tools.Logbook()

        ## Set random seed ##
        random.seed(self.randomseed)

        ##Generation 0

        # # #Get the fitness value
        pdbfile_path = sets.CONFIG_PATH + self.scenario + "/" +self.pdbfile 

        src = pdbfile_path
        dst = output_path + "g0_00.pdb"
        

        #copy the original individual pdb file 
        shutil.copyfile(src, dst)

        #Create Mutated individuals 
        mut_rate = 0.3

        # # ## Create Population 0 
        gen = 0

        offspring_output_pdbfiles=[]
        argument=[]


        if checkpoint:
            try:
                # Load the checkpoint file
                logging.debug("Loading from checkpoint")

                # Load the checkpoint
                with open(f"{self.output}/checkpoint.pkl", 'rb') as cp_file:
                    cp = pickle.load(cp_file)

                pop = cp['population']
                ind0 = cp['ind0']
                ref_point = cp['ref_point']
                hof = cp['hof']
                ngen = cp['generation'] + 1
                record = cp['record']
                hv = cp['hv']
                logbook = cp['logbook']
                random.setstate(cp['rndstate'])
                population_output_pdbfiles = cp['population_output_pdbfiles']
            except FileNotFoundError as fne:
                logging.critical(f"Checkpoint file not found.\n{fne}")
                exit(1)
            
            except EOFError as ee:
                logging.critical(f"Checkpoint file corrupted.\n{ee}")
                exit(1)
        else:
            logging.info("Starting new run")

            #Delete old hallofame and population csv
            if(os.path.exists(SAVEHALLOFFAME_CSV)):os.remove(SAVEHALLOFFAME_CSV)
            if(os.path.exists(SAVEPOPGEN_CSV)):os.remove(SAVEPOPGEN_CSV)

            # HallofFame
            hof = tools.ParetoFront()
    
            #Create individual0
            logging.info(f"Creating Individual 0 from pdbfile: {self.pdbfile}")
            aa0 = self.my_protein_problem.create_individual0(self.pdbfile)
            logging.info(f"Individual-Original : {aa0}")
            ind0 = creator.Individual(aa0)  # Instantiate the Individual with fixed values
            ind0.father = "Original"
            ind0.pdb = dst
            
            logging.debug(f"ind0: {ind0}")
    
            fitness_indv0 = self.my_protein_problem.fitness(dst)
            
            # # # Set initial fitness value
            fitness_values = tuple(fitness_indv0[fitness_idx] for fitness_idx in self.fitness_idxs)
            ifitness_values = (*fitness_values, 0)
            logging.debug(f"Fitness Indv0: {ifitness_values}")

            ## Create parameters to run in parallel
            for i in range(1,popsize):
                if i < 10:
                    num_ind="0"+str(i)
                else:
                    num_ind=i
    
                output_file_path = self.output + "/g" + str(gen) + "/" + "g"+ str(gen) +"_" + str(num_ind) + ".pdb"
                offspring_output_pdbfiles.append(output_file_path)
                argument.append((pdbfile_path,output_file_path,mut_rate))
    
            logging.debug(f"Arguments:{argument}")
            
            #Mutate in parallel
            offspring_list = []
            offspring_list = self.my_protein_problem.mutate_population(argument)
    
            #Evaluate in parallel
            fitness = self.my_protein_problem.fitnessPop(offspring_output_pdbfiles)
        
            #Create Population 0
            ngen = 1
            i=0
            pop = []
            ##Add original individual
            ind0.fitness.values = (*fitness_values, 0)
            father_complete_sequence = self.get_complete_sequence(ind0.pdb, self.ligand_chain)
            ll_father = self.esm2.get_esm_ll(father_complete_sequence)

            pop.append(ind0)
            for indiv, pdb_file in zip(offspring_list, offspring_output_pdbfiles):
                ind = creator.Individual(indiv)
                ind.father = '0'
                ind.pdb = pdb_file
                ind.nmut = hamming_distance(ind0, indiv)
                ind.id = f'0-{ind.id}'

                
                ind_complete_sequence = self.get_complete_sequence(pdb_file, self.ligand_chain)
                ll = self.esm2.get_esm_ll(ind_complete_sequence)

                delta_ll = ll - ll_father
                fitness_values = tuple(fitness[i][fitness_idx] for fitness_idx in self.fitness_idxs)
                ind.fitness.values = (*fitness_values, delta_ll)
                logging.debug(f"Fitness: {ind.fitness.values}")
                pop.append(ind)
                i=i+1

            #Update HallOfFame
            hof.update(pop)

            # Compute Hypervolume
            front = sortNondominated(pop, k=len(pop), first_front_only=True)[0]
            ref_point = np.max([ind.fitness.values for ind in pop], axis=0) + 0.1
            hv_first = hypervolume(front, ref_point)
            hv = [hv_first]
    
            ## Add original individual pdb file to begining of the output file list            
            offspring_output_pdbfiles.insert(0,dst)
    
            #initial population statistics
            record = stats.compile(pop)
            logging.info("stats: %s", record)
            logbook.record(gen=0, **record)

    
            ###Save population to text file AA representation###
            gendir = self.output + "/g" + str(gen) + "/"
            with open(gendir+"/pop_g" + str(gen) + "_AA.txt", "w") as output_file:
                #printpopulation and fitness
                i=0  
                for ind in pop:
                    fit=ind.fitness.values
                    output_file.write((str(ind)) + " " + str(fit) + " "+ str(offspring_output_pdbfiles[i])+'\n')   
                    i=i+1    
            output_file.close() 
    
            population_output_pdbfiles = offspring_output_pdbfiles.copy()
            save_population_to_csv(pop, gen, SAVEPOPGEN_CSV)
        
        # ############################ 
        # ## Main evolutionary loop ##
        # ############################
        
        for gen in range(ngen,ngenerations):
            logging.info(f"-- Generation {gen} --")
            # Sort the population by fitness and select the elite individuals

            gendir = self.output + "/g" + str(gen)

            if not os.path.exists(gendir):
                try:
                   os.makedirs(gendir)
                except OSError:
                   logging.error("Creation of the directory %s failed" % gendir)

            ## -------- ##
            ## Mutation ##
            ## -------- ##
             
            #Empty offspring pdb files list
            offspring_output_pdbfiles.clear()
            
            argument.clear()
            i=0
            for mutant in pop:
                ######################
                #### Always Mutate ###
                ######################
                #Create Mutant
                indiv_to_mutate_pdb = population_output_pdbfiles[i]
                logging.info(f"pdbfile: {indiv_to_mutate_pdb}")
                logging.info(f"Individual-to-Mutate: {mutant}")
                if i < 10:
                    num_ind="0"+str(i)
                else:
                    num_ind=i

                output_file_path = self.output + "/tmp" + "/" + "temp_g"+ str(gen) +"_" + str(num_ind) + ".pdb"

                #Add to generation_pdb_files the new ones (pop + offspring)
                
                population_output_pdbfiles.append(output_file_path)
                offspring_output_pdbfiles.append(output_file_path)
                argument.append((indiv_to_mutate_pdb,output_file_path,mut_rate))
                i=i+1

            logging.info(f"Population_pdb_files: {population_output_pdbfiles}")
            logging.info("Before mutation in parallel")
            logging.info(f"Argument-gen{gen}: {argument}")
             
            #Mutate in parallel
            offspring = []
            offspring = self.my_protein_problem.mutate_population(argument)

            #Evaluate in parallel
            fitness = self.my_protein_problem.fitnessPop(offspring_output_pdbfiles)

            logging.info("Population before offspring")
            logging.info(f"Pop from Generation {gen} = {offspring}")
            logging.info(f"Fitness Pop = {fitness}")
            logging.info(f"Population output pdb files: {population_output_pdbfiles}")

            #Assign the mutants to offspring and replace as new population
            #Create new population 
            i=0
            pop_copy = pop.copy()
            ##Add individuals to population
            for indiv, pdb_file in zip(offspring, offspring_output_pdbfiles):
                ind = creator.Individual(indiv)
                ind.father = pop_copy[i].id  # Set the father of the individual
                ind.id = f'{gen}-{ind.id}'
                ind.nmut = hamming_distance(ind0, indiv)
                ind.pdb = pdb_file

                father_complete_sequence = self.get_complete_sequence(pop_copy[i].pdb, self.ligand_chain)
                ind_complete_sequence = self.get_complete_sequence(ind.pdb, self.ligand_chain)
                ll_father = self.esm2.get_esm_ll(father_complete_sequence)
                ll = self.esm2.get_esm_ll(ind_complete_sequence)

                delta_ll = ll - ll_father
                fitness_values = tuple(fitness[i][fitness_idx] for fitness_idx in self.fitness_idxs)
                ind.fitness.values = (*fitness_values, delta_ll)
                #Add the new individuals to population
                pop.append(ind)
                i=i+1

            logging.info(f"Pop = {pop} size: {len(pop)}")

            ##Select the new individual from new pop 
            offspring, selected_indices = self.unique_offspring(pop, toolbox.select, popsize)
            offspring = list(map(toolbox.clone, offspring))
            logging.info(f"Selected individuals from pop: {selected_indices}" )

            #Replace the selected individuals to new population
            pop.clear()
            pop = offspring 

            logging.info("New Population with elite and offspring ")
            logging.info(f"Pop = {pop} size: {len(pop)}")

            # Compute hypervolume
            front = sortNondominated(pop, k=len(pop), first_front_only=True)[0]

            hv_pop = hypervolume(front, ref_point)
            hv.append(hv_pop)

            ##Add to population elite individuals
            
            # Gather all the fitnesses in one list and printthe stats
            fits = [ind.fitness.values[0] for ind in pop]
            logging.info(f"Fitness : {fits}")
            
            new_generation_output_pdbfiles = [] 


            # Save the rest of pdb files from population    
            for i_, idx in enumerate(selected_indices):
                
                ##pdb files
                src1 = population_output_pdbfiles[idx]
                dst1 = self.output + "/g" + str(gen) + "/" + "g"+ str(gen) +"_" + str(i_) + ".pdb"
                offspring[i_].pdb = dst1
                #copy the original individual pdb file 
                shutil.copyfile(src1, dst1)

                logging.info("generation pdb files")
                logging.info(f"{src1} ---> {dst1}")                                

                new_generation_output_pdbfiles.append(dst1)

                ##score files- .sc(Energy Summary)
                src2 = str(src1) + ".sc"
                dst2 = str(dst1) + ".sc"
                shutil.copyfile(src2, dst2)

                logging.info("generation sc files")
                logging.info(f"{src2} ---> {dst2}")

            
            # Update HallOfFame
            hof.update(pop)

            logging.info(f"new_generation_output_pdbfiles={new_generation_output_pdbfiles}")
            ##Update population output files from new generation ###
            population_output_pdbfiles.clear()
            population_output_pdbfiles = new_generation_output_pdbfiles.copy()


            ###Save population to text file AA representation###
            with open(gendir+"/pop_g" + str(gen) + "_AA.txt", "w") as output_file:
                i=0    
                for ind in pop:
                    fit=ind.fitness.values 
                    output_file.write((str(ind)) + " " + str(fit) + " "+ str(self.output + "/g" + str(gen) + "/" + "g"+ str(gen) +"_" + str(i) + ".pdb")+'\n')   
                    i=i+1    
            output_file.close() 



            ###Compile Statistics
            record = stats.compile(pop)
            logging.info(f"stats: {record}")
            logbook.record(gen=gen, **record)

            save_population_to_csv(pop, gen, SAVEPOPGEN_CSV)

            # # Save the logbook to a file
            if gen % freq == 0 or gen == ngenerations - 1:
                logging.debug(f"Checkpoint reached at generation {gen}, saving logbook.")
                cp = dict(
                    population=pop,
                    ind0=ind0,
                    ref_point=ref_point,
                    hof=hof,
                    hv=hv,
                    generation=gen,
                    record=record,
                    logbook=logbook,
                    rndstate=random.getstate(),
                    population_output_pdbfiles=population_output_pdbfiles
                )
                with open(f"{self.output}/checkpoint.pkl", 'wb') as cp_file:
                    pickle.dump(cp, cp_file)

            ## Remove tmp files
            temp_files_path = os.path.join(self.output, "tmp")
            for file_temp in os.listdir(temp_files_path):
                file_temp_path = os.path.join(temp_files_path, file_temp)
                os.remove(file_temp_path)
            

        logging.info("-- End of Evolution --")


        # CSV to JSON
        csvToTree(SAVEPOPGEN_CSV, SAVEPOPGEN_CSV.replace('.csv', '.json'))
        save_HallofFame(hof, SAVEHALLOFFAME_CSV)

        # Group .sc files into .csv
        read_scfiles(self.output,SCFILE_CSV)
        
        logging.info("-- Saving Evolution Statistics--")
        logbook.header = "gen", "avg", "min", "max", "std"
        with open(SAVE_STATISTICS_CSV, "w") as stat_file:
            header = ["gen"] + [f"{stat}_obj{i}" for stat in ["avg", "min", "max", "std"] for i in range(nobj)]
            stat_file.write(",".join(header) + "\n")
        
            for entry in logbook:
                row = [str(entry["gen"])]
                for stat in ["avg", "min", "max", "std"]:
                    values = entry[stat] 
                    row += [str(v) for v in values]
                stat_file.write(",".join(row) + "\n")

        
        logging.info("--Saving Hypervolumes per generation--")
        with open(SAVEHYPERVOLUME_CSV, mode="w", newline='') as hv_file:
            writer = csv.writer(hv_file)
            writer.writerow(["Gen","HV"])
            for i,v in enumerate(hv):
                writer.writerow([i,v])