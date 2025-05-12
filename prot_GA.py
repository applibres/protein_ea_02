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


class deap_sga_protein:

    def __init__(self, scenario, algoritm_params, sim_params, output, randomseed):
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
        pyrosetta.init()




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
            print("Warning: Could not find enough unique individuals, filled with randoms.")
            while len(unique_inds) < k:
                rand_idx = random.randint(0, len(population) - 1)
                if rand_idx not in selected:
                    unique_inds.append(population[rand_idx])
                    indices.append(rand_idx)
                    selected.add(rand_idx)

        return unique_inds, indices

    def convert_to_characters(self,individual):
         aa_list=[]
         for i in individual:
            aa = chr(i+64)
            aa_list.append(aa)
         return aa_list   

    # Custom simple evolutionary algorithm
    def run(self):
        #Capture parameters

        #Algorithm Params

        FITNESS_INDEX = 2 #index to capture the fitness value
        popsize = self.algoritm_params['popsize']
        ngenerations = self.algoritm_params['gen']
        nobj = self.algoritm_params['obj']
        mutprob = self.algoritm_params['mutp']
        print("###Algorithm Parameters###")
        print (f"popsize = {popsize}")
        print (f"ngenerations = {ngenerations}")
        print (f"nobj = {nobj}")
        print (f"mut prob = {mutprob}")
    
        
        #Create tmp directory
        output_path = self.output+"/tmp/"
        print (f"output_path = {output_path}")
        

        #Create directory if not exist
        if os.path.isdir(output_path):
            print (f"{output_path} directory already exist")
        else:
            os.makedirs(output_path)
            print (f"{output_path} directory created")
    


        output_path = self.output+"/g0/"
        print (f"output_path = {output_path}")
        

        #Create directory if not exist
        if os.path.isdir(output_path):
            print (f"{output_path} directory already exist")
        else:
            os.makedirs(output_path)
            print (f"{output_path} directory created")





        ##Declare FitnessMinimization and Individual
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))  # Minimization problem
        creator.create("Individual", list, fitness=creator.FitnessMin)

        # Initialize toolbox
        toolbox = base.Toolbox() 
        
        ##Selection operator implemented by framework##
        toolbox.register("select", tools.selTournament, tournsize=3) 
        #toolbox.register("select", tools.selBest)       
        #toolbox.register("select", tools.selRoulette)
        #toolbox.register("select", tools.selRandom)
        
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

        #Create individual0
        print(f"Creating Individual 0 from pdbfile: {self.pdbfile}")
        indiv0, aa0 = self.my_protein_problem.create_individual0(self.pdbfile)
        print (f"Individual-Original : {aa0}")
        #print (f"Indiv size : {len(indiv0)}")
        ind0 = creator.Individual(indiv0)  # Instantiate the Individual with fixed values
        
        print (f"ind0: {ind0}")
        #print (f"Indiv size : {len(ind0)}")
        
        # # #Get the fitness value
        pdbfile_path = sets.CONFIG_PATH + self.scenario + "/" +self.pdbfile 
        
        src = pdbfile_path
        #print("src: ",src)         
        dst = output_path + "g0_00.pdb"
        #copy the original individual pdb file 
        shutil.copyfile(src, dst)

        fitness = self.my_protein_problem.fitness(dst)
        
        # # # Set initial fitness value
        ind0.fitness.values = (fitness[FITNESS_INDEX],)  
        #print (f"Fitness: {fitness}")
        print (f"Fitness Indv0: {ind0.fitness.values}")

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

        ## Create parameters to run in parallel
        for i in range(1,popsize):
            #population_pdbfiles.append(pdbfile_path)
            if i < 10:
                num_ind="0"+str(i)
            else:
                num_ind=i

            output_file_path = self.output + "/g" + str(gen) + "/" + "g"+ str(gen) +"_" + str(num_ind) + ".pdb"
            offspring_output_pdbfiles.append(output_file_path)
            argument.append((pdbfile_path,output_file_path,mut_rate))

        #print (f"Population pdb files: {population_pdbfiles}")
        #print (f"Population output pdb files: {population_output_pdbfiles}")

        print(f"Arguments:{argument}")
        
        #Mutate in parallel
        offspring_list = []
        offspring_list = self.my_protein_problem.mutate_population(argument)

        #Evaluate in parallel
        fitness = self.my_protein_problem.fitnessPop(offspring_output_pdbfiles)

        #print(f"Pop = {population}")
        #print(f"Fitness Pop = {fitness}")
        #print (f"Population output pdb files: {population_output_pdbfiles}")
        
        #Create Population 0
        i=0
        pop = []
        ##Add original individual
        pop.append(ind0)
        for indiv in offspring_list:
            ind = creator.Individual(indiv) 
            ind.fitness.values = (fitness[i][FITNESS_INDEX],)  
            print (f"Fitness: {ind.fitness.values}")
            pop.append(ind)
            i=i+1

        ## Add original individual pdb file to begining of the output file list            
        offspring_output_pdbfiles.insert(0,dst)

        #print (f"Population 0: {pop}")

        #initial population statistics
        record = stats.compile(pop)
        print("stats: ", record)
        logbook.record(gen=0, **record)

 
        ###Save population to text file integer representation###
        gendir = self.output + "/g" + str(gen) + "/"
        with open(gendir+"/pop_g" + str(gen) + ".txt", "w") as output_file:
            #print population and fitness
            #cp_file.write('\n'.join(map(str, population)))
            #print hall of fame
            #cp_file.write(halloffame)
            i=0    
            for ind in pop:
                fit=ind.fitness.values
                #output_file.write(str(ind) + " " + str(fit)+'\n')   
                output_file.write((str(ind)) + " " + str(fit) + " "+ str(offspring_output_pdbfiles[i])+'\n')   
                i=i+1    
        output_file.close() 


        ###Save population to text file AA representation###
        gendir = self.output + "/g" + str(gen) + "/"
        with open(gendir+"/pop_g" + str(gen) + "_AA.txt", "w") as output_file:
            #print population and fitness
            #cp_file.write('\n'.join(map(str, population)))
            #print hall of fame
            #cp_file.write(halloffame)
            i=0 
            ind_AA=[]   
            for ind in pop:
                fit=ind.fitness.values
                #output_file.write(str(ind) + " " + str(fit)+'\n') 
                ind_AA=self.convert_to_characters(ind)  
                output_file.write((str(ind_AA)) + " " + str(fit) + " "+ str(offspring_output_pdbfiles[i])+'\n')   
                i=i+1    
        output_file.close() 


        #offspring = toolbox.select(pop, popsize)
        #offspring = list(map(toolbox.clone, offspring))
        #print (f"offspring: {offspring}")                 
 
        population_output_pdbfiles = offspring_output_pdbfiles.copy()
        
        # ############################ 
        # ## Main evolutionary loop ##
        # ############################
        
        for gen in range(1,ngenerations):
            print(f"-- Generation {gen} --")
            # Sort the population by fitness and select the elite individuals
            #elite_individuals = sorted(population, key=lambda ind: ind.fitness.values, reverse=False)[:elite_size]

            gendir = self.output + "/g" + str(gen)

            if not os.path.exists(gendir):
                try:
                   os.makedirs(gendir)
                except OSError:
                   print ("Creation of the directory %s failed" % gendir)
                #else:
                #    print ("Successfully created the directory %s " % path)
                               
            ##################################
            ## Select the elite individuals ##
            ##################################
            elite_inds = tools.selBest(pop, elite_size)
            elite_inds = [toolbox.clone(ind) for ind in elite_inds]  # Clone to avoid overwriting
            elite_indexes = [pop.index(ind) for ind in elite_inds] 
            print (f"Elite Individuals {elite_inds}")
            print (f"Elite Indexes {elite_indexes}")
            elite_pdb_files = [population_output_pdbfiles[i] for i in elite_indexes]
            print (f"Elite individual file {elite_pdb_files}")
            ##Check Fitness
            for e_ind in elite_inds:
                print(f"e_ind={e_ind} fitness:{e_ind.fitness.values}")

            #####

            
            #offspring_pdb_files = population_output_pdbfiles.copy()           
            #generation_pdb_files = population_output_pdbfiles.copy()
            #generation_pdb_files = []
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
                print(f"pdbfile: {indiv_to_mutate_pdb}")
                print (f"Individual-to-Mutate: {mutant}")
                if i < 10:
                    num_ind="0"+str(i)
                else:
                    num_ind=i

                #output_file_path = self.output + "/g" + str(gen) + "/" + "temp_g"+ str(gen) +"_" + str(num_ind) + ".pdb"
                output_file_path = self.output + "/tmp" + "/" + "temp_g"+ str(gen) +"_" + str(num_ind) + ".pdb"

                 #Add to generation_pdb_files the new ones (pop + offspring)
                 #generation_pdb_files.append(output_file_path)
                
                population_output_pdbfiles.append(output_file_path)
                offspring_output_pdbfiles.append(output_file_path)
                argument.append((indiv_to_mutate_pdb,output_file_path,mut_rate))                 
                i=i+1

                #print(f"Generation_pdb_files: {generation_pdb_files}")
            print(f"Population_pdb_files: {population_output_pdbfiles}")
            print("Before mutation in parallel")
            print(f"Argument: {argument}")
             
            #Mutate in parallel
            offspring = []
            offspring = self.my_protein_problem.mutate_population(argument)

            
             #print("Before Evaluation in paralell")
             #print(f"Files to get fitness: {population_output_pdbfiles}") 

            #Evaluate in parallel
            fitness = self.my_protein_problem.fitnessPop(offspring_output_pdbfiles)

            print("Population before offspring")
            print(f"Pop from Generation {gen} = {offspring}")
            print(f"Fitness Pop = {fitness}")
            print (f"Population output pdb files: {population_output_pdbfiles}")

            #Assign the mutants to offspring and replace as new population
            #Create new population 
            i=0
            #pop = []
            
            ##Add individuals to population
            for indiv in offspring:
                ind = creator.Individual(indiv) 
                ind.fitness.values = (fitness[i][FITNESS_INDEX],)  
                #print (f"Fitness: {ind.fitness.values}")
                #Add the new individuals to population
                pop.append(ind)
                i=i+1

            print(f"Pop = {pop} size: {len(pop)}")

            ##Select the new individual from new pop 
            offspring, selected_indices = self.unique_offspring(pop, toolbox.select, popsize-elite_size)
            offspring = list(map(toolbox.clone, offspring))
            #print (f"offspring: {offspring} size: {len(offspring)} ") 
            print (f"Selected individuals from pop: {selected_indices}" )

            #Replace the selected individuals to new population
            pop.clear()
            pop = elite_inds + offspring  

            print("New Population with elite and offspring ")
            print(f"Pop = {pop} size: {len(pop)}")

        #     ##Add to population elite individuals
        #     #pop.append(elite_inds)
            
            # Gather all the fitnesses in one list and print the stats
            fits = [ind.fitness.values[0] for ind in pop]
            print(f"Fitness : {fits}")

            
            new_generation_output_pdbfiles = [] 
            
            #Save elite files to directory
             
            i=0
            for elite_file in elite_pdb_files: 
                
                src1 = elite_file               
                ## 1-Change names
                dst1 = self.output + "/g" + str(gen) + "/" + "g"+ str(gen) +"_" + str(i) + ".pdb"

                ## 2-Copy to directory
                shutil.copyfile(src1, dst1)
                
                print("Elite pdb files")
                print(f"{src1} ---> {dst1}")                                
                
                #Append elite individuals
                new_generation_output_pdbfiles.append(dst1)

                ##score files- .sc(Energy Summary)
                src2 = str(src1) + ".sc"
                dst2 = str(dst1) + ".sc"
                shutil.copyfile(src2, dst2)                

                print("Elite score files")
                print(f"{src2} ---> {dst2}")                                

                i = i + 1 


            # Save the rest of pdb files from population    
            i=elite_size
            for idx in selected_indices:
                
                ##pdb files
                src1 = population_output_pdbfiles[idx]
                #print("src: ",src)         
                dst1 = self.output + "/g" + str(gen) + "/" + "g"+ str(gen) +"_" + str(i) + ".pdb"
                #dst = output_path + "g0_00.pdb"
                #copy the original individual pdb file 
                shutil.copyfile(src1, dst1)

                print("generation pdb files")
                print(f"{src1} ---> {dst1}")                                

                new_generation_output_pdbfiles.append(dst1)

                ##score files- .sc(Energy Summary)
                src2 = str(src1) + ".sc"
                dst2 = str(dst1) + ".sc"
                shutil.copyfile(src2, dst2)

                print("generation sc files")
                print(f"{src2} ---> {dst2}")                                

                i=i+1


            print(f"new_generation_output_pdbfiles={new_generation_output_pdbfiles}")
            ##Update population output files from new generation ###
            population_output_pdbfiles.clear()
            population_output_pdbfiles = new_generation_output_pdbfiles.copy()    
                
            ###Save population to text file integer representation###
            with open(gendir+"/pop_g" + str(gen) + ".txt", "w") as output_file:
                i=0    
                for ind in pop:
                    fit=ind.fitness.values
                    #output_file.write(str(ind) + " " + str(fit)+'\n')   
                    #output_file.write((str(ind)) + " " + str(fit) + " "+ str(population_output_pdbfiles[i])+'\n')   
                    output_file.write((str(ind)) + " " + str(fit) + " "+ str(self.output + "/g" + str(gen) + "/" + "g"+ str(gen) +"_" + str(i) + ".pdb")+'\n')   

                    i=i+1    
            output_file.close() 

            ###Save population to text file AA representation###
            ind_AA=[]
            with open(gendir+"/pop_g" + str(gen) + "_AA.txt", "w") as output_file:
                i=0    
                for ind in pop:
                    fit=ind.fitness.values
                    #output_file.write(str(ind) + " " + str(fit)+'\n')   
                    ind_AA=self.convert_to_characters(ind)  
                    output_file.write((str(ind_AA)) + " " + str(fit) + " "+ str(self.output + "/g" + str(gen) + "/" + "g"+ str(gen) +"_" + str(i) + ".pdb")+'\n')   
                    i=i+1    
            output_file.close() 



            ###Compile Statistics
            record = stats.compile(pop)
            print("stats: ", record)
            logbook.record(gen=gen, **record)

        print("-- End of Evolution --")

        ## Remove tmp files
        #remove_command = "rm -r "+ self.output+"/tmp/*"
        #os.system(remove_command)
        
        print("-- Saving Evolution Statistics--")
        logbook.header = "gen", "avg", "min", "max", "std"
        gen, avg, min, max, std = logbook.select("gen", "avg", "min", "max", "std") 
        
        with open(self.output+"/run"+str(self.randomseed)+"_statistics.csv", "w") as stat_file:                    
            genrt, avg, min, max, std = logbook.select("gen", "avg", "min", "max", "std") 
            stat_file.write("gen,avg,min,max,std\n")
            for i in range(0, int(ngenerations)): 
                stat_file.write(str(genrt[i])+","+str(avg[i])+","+str(min[i])+","+str(max[i])+","+str(std[i])+"\n")   
        stat_file.close() 


