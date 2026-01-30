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

from prot_interface.prot_problemI import *
import prot_interface.prot_problemI as problem
import prot_interface.prot_settingsI as sets
import pyrosetta
import os
import shutil
import numpy
import random
from prot_interface.logging_config import setup_logging
from utils import csvToTree, hamming_distance, save_population_to_csv, read_scfiles
import logging
import pickle

# Initialize logging before anything else
setup_logging()
logger = logging.getLogger(__name__)

CSV_FILE = 'individuals.csv'

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

class deap_sga_protein:

    def __init__(self, scenario, algoritm_params, sim_params, fitness_idx, output, randomseed):
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
        self.fitness_idx = fitness_idx[0]
        self.output = output
        self.randomseed = randomseed
        self.scenario = scenario
        
        #Simulation Params
        self.partners = self.sim_params['partners']
        self.ligand_chain = self.sim_params['ligand_chain']
        self.pdbfile = self.sim_params['pdbfile']
        
        #Instantiate problem
        self.my_protein_problem = problem.prot_problem(self.scenario,self.partners,self.ligand_chain)


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



    def unique_offspring(self, population, selection_func, k, elite_idx):
        """
        Select k unique individuals using the given selection function, 
        returning both selected individuals and their original indices.

        :param population: List of individuals.
        :param selection_func: Selection function like toolbox.select.
        :param k: Number of individuals to select.
        :return: List of selected individuals (no duplicates), and their indices.
        """
        selected = set(elite_idx)
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
        SAVE_STATISTICS_CSV = os.path.join(PATH_STATISTICS, f"run{self.randomseed}_statistics.csv")

        if os.path.isfile(SAVEPOPGEN_CSV): os.remove(SAVEPOPGEN_CSV)
        if os.path.isfile(SCFILE_CSV): os.remove(SCFILE_CSV)
        if os.path.isfile(SAVE_STATISTICS_CSV): os.remove(SAVE_STATISTICS_CSV)

        os.makedirs(PATH_STATISTICS, exist_ok=True)
        
        popsize = self.algoritm_params['popsize']
        ngenerations = self.algoritm_params['gen']
        mutprob = self.algoritm_params['mutp']
        logging.info("###Algorithm Parameters###")
        logging.info(f"popsize = {popsize}")
        logging.info(f"ngenerations = {ngenerations}")
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
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))  # Minimization problem
        creator.create("Individual", CustomIndividual)

        # Initialize toolbox
        toolbox = base.Toolbox() 
        
        ##Selection operator implemented by framework##
        toolbox.register("select", tools.selTournament, tournsize=3)
        
        # Statistics
        stats = tools.Statistics(key=lambda ind: ind.fitness.values)
        stats.register("avg", numpy.mean)
        stats.register("std", numpy.std)
        stats.register("min", numpy.min)
        stats.register("max", numpy.max)

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

        ##Number of elite individuals per generation
        elite_size = int(0.1 * popsize)
        if (elite_size < 1):
            elite_size = 1      

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
                ngen = cp['generation'] + 1
                record = cp['record']
                logbook = cp['logbook']
                random.setstate(cp['rndstate'])
                population_output_pdbfiles = cp['population_output_pdbfiles']
            except FileNotFoundError as fne:
                logging.critical(f"Checkpoint file not found.\n{fne}")
                return 
            
            except EOFError as ee:
                logging.critical(f"Checkpoint file corrupted.\n{ee}")
                return 
        else:
            logging.info("Starting new run")

            #Delete old CSV_FILE
            if(os.path.exists(self.output + f"/{CSV_FILE}")):os.remove(self.output + f"/{CSV_FILE}")
    
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
            ind0.fitness.values = (fitness_indv0[self.fitness_idx],)  
            logging.debug(f"Fitness Indv0: {ind0.fitness.values}")

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
            pop.append(ind0)
            for indiv, pdb_file in zip(offspring_list, offspring_output_pdbfiles):
                ind = creator.Individual(indiv)
                ind.father = '0'
                ind.pdb = pdb_file
                ind.id = f'0-{ind.id}'
                ind.nmut = hamming_distance(ind0, indiv)
                ind.fitness.values = (fitness[i][self.fitness_idx],)  
                logging.debug(f"Fitness: {ind.fitness.values}")
                pop.append(ind)
                i=i+1
    
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
                               
            ##################################
            ## Select the elite individuals ##
            ##################################
            elite_inds = tools.selBest(pop, elite_size)
            elite_indexes = [pop.index(ind) for ind in elite_inds] 
            elite_inds = [toolbox.clone(ind) for ind in elite_inds]  # Clone to avoid overwriting
            logging.info(f"Elite Individuals {elite_inds}")
            logging.info(f"Elite Indexes {elite_indexes}")
            elite_pdb_files = [population_output_pdbfiles[i] for i in elite_indexes]
            logging.info(f"Elite individual file {elite_pdb_files}")
            ##Check Fitness
            for e_ind in elite_inds:
                logging.info(f"e_ind={e_ind} fitness:{e_ind.fitness.values}")

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
            for indiv in offspring:
                ind = creator.Individual(indiv)
                ind.father = pop_copy[i].id  # Set the father of the individual
                ind.id = f'{gen}-{ind.id}'
                ind.nmut = hamming_distance(ind0, indiv)
                ind.fitness.values = (fitness[i][self.fitness_idx],)  
                #Add the new individuals to population
                pop.append(ind)
                i=i+1

            logging.info(f"Pop = {pop} size: {len(pop)}")

            ##Select the new individual from new pop 
            offspring, selected_indices = self.unique_offspring(pop, toolbox.select, popsize-elite_size, elite_indexes)
            offspring = list(map(toolbox.clone, offspring))
            logging.info(f"Selected individuals from pop: {selected_indices}" )

            #Replace the selected individuals to new population
            pop.clear()
            pop = elite_inds + offspring 

            logging.info("New Population with elite and offspring ")
            logging.info(f"Pop = {pop} size: {len(pop)}")

            ##Add to population elite individuals
            
            # Gather all the fitnesses in one list and printthe stats
            fits = [ind.fitness.values[0] for ind in pop]
            logging.info(f"Fitness : {fits}")
                

            
            new_generation_output_pdbfiles = [] 
            
            #Save elite files to directory
             
            i=0
            for elite_file in elite_pdb_files:

                if i < 10:
                    num_ind="0"+str(i)
                else:
                    num_ind=i

                src1 = elite_file               
                ## 1-Change names
                dst1 = self.output + "/g" + str(gen) + "/" + "g"+ str(gen) +"_" + str(num_ind) + ".pdb"
                elite_inds[i].pdb = dst1

                ## 2-Copy to directory
                shutil.copyfile(src1, dst1)
                
                logging.info("Elite pdb files")
                logging.info(f"{src1} ---> {dst1}")                                
                
                #Append elite individuals
                new_generation_output_pdbfiles.append(dst1)

                ##score files- .sc(Energy Summary)
                src2 = str(src1) + ".sc"
                dst2 = str(dst1) + ".sc"
                shutil.copyfile(src2, dst2)                

                logging.info("Elite score files")
                logging.info(f"{src2} ---> {dst2}")                                

                i = i + 1 


            # Save the rest of pdb files from population    
            i=elite_size
            for i_, idx in enumerate(selected_indices):
                
                if i < 10:
                    num_ind="0"+str(i)
                else:
                    num_ind=i

                ##pdb files
                src1 = population_output_pdbfiles[idx]
                dst1 = self.output + "/g" + str(gen) + "/" + "g"+ str(gen) +"_" + str(num_ind) + ".pdb"
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

                i=i+1


            logging.info(f"new_generation_output_pdbfiles={new_generation_output_pdbfiles}")
            ##Update population output files from new generation ###
            population_output_pdbfiles.clear()
            population_output_pdbfiles = new_generation_output_pdbfiles.copy()    

            ###Save population to text file AA representation###
            with open(gendir+"/pop_g" + str(gen) + "_AA.txt", "w") as output_file:
                i=0    
                for ind in pop:
                    if i < 10:
                        num_ind="0"+str(i)
                    else:
                        num_ind=i
                    fit=ind.fitness.values
                    output_file.write((str(ind)) + " " + str(fit) + " "+ str(self.output + "/g" + str(gen) + "/" + "g"+ str(gen) +"_" + str(num_ind) + ".pdb")+'\n')
                    i=i+1
            output_file.close()



            ###Compile Statistics
            record = stats.compile(pop)
            logging.info("stats: ", record)
            logbook.record(gen=gen, **record)

            save_population_to_csv(pop, gen, SAVEPOPGEN_CSV)

            # # Save the logbook to a file
            if gen % freq == 0 or gen == ngenerations - 1:
                logging.debug(f"Checkpoint reached at generation {gen}, saving logbook.")
                cp = dict(
                    population=pop,
                    ind0=ind0,
                    generation=gen,
                    record=record,
                    logbook=logbook,
                    rndstate=random.getstate(),
                    population_output_pdbfiles=population_output_pdbfiles
                )
                with open(f"{self.output}/checkpoint.pkl", 'wb') as cp_file:
                    pickle.dump(cp, cp_file)

            # This mutation block is necesary to create the .pdb.txt files for the last generation
            if gen == (ngenerations - 1):
                argument.clear()
                i=0
                for mutant in pop:
                    indiv_to_mutate_pdb = population_output_pdbfiles[i]
                    output_file_path = self.output + "/tmp" + "/" + "temp_g"+ str(gen) +"_" + str(i) + ".pdb"
                    argument.append((indiv_to_mutate_pdb,output_file_path,mut_rate))                 
                    i=i+1
        
                offspring = self.my_protein_problem.mutate_population(argument)
    
            ## Remove tmp files
            temp_files_path = os.path.join(self.output, "tmp")
            for file_temp in os.listdir(temp_files_path):
                file_temp_path = os.path.join(temp_files_path, file_temp)
                os.remove(file_temp_path)
    

        logging.info("-- End of Evolution --")

        # CSV to JSON
        csvToTree(SAVEPOPGEN_CSV, SAVEPOPGEN_CSV.replace('.csv', '.json'))

        # Group .sc files into .csv
        read_scfiles(self.output,SCFILE_CSV)
        
        logging.info("-- Saving Evolution Statistics--")
        logbook.header = "gen", "avg", "min", "max", "std"
        gen, avg, min, max, std = logbook.select("gen", "avg", "min", "max", "std") 
        
        with open(SAVE_STATISTICS_CSV, "w") as stat_file:                    
            genrt, avg, min, max, std = logbook.select("gen", "avg", "min", "max", "std") 
            stat_file.write("gen,avg,min,max,std\n")
            for i in range(0, int(ngenerations)): 
                stat_file.write(str(genrt[i])+","+str(avg[i])+","+str(min[i])+","+str(max[i])+","+str(std[i])+"\n")   
        stat_file.close() 

