#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 30/01/2025

@author: Rolando Armas
Yachay Tech University
Phage Therapy Group
"""

''' Call to protein evolutionary algorithm '''

import sys
import timeit
import prot_interface.prot_parserI as parser
#import sga_mut_protein as sga
#import sga_mut_protein02 as sga
import prot_GA as sga

# Main function to run the genetic algorithm
def main():

    # command line arguments
    if len(sys.argv) != 5:
        print('usage: ', sys.argv[0], '<scenario> <sea/moea> <sim params> <algo params>')
        sys.exit(-1)

    print("Init Main")
    
    SCENARIO       = sys.argv[1] # Scenario Name
    # Given the scenario, set the pdb file and mutation limits
    
    ALGO_NAME      = sys.argv[2] # Single Objective (sea) / Multi Objective (moea)   
    SIM_PARAM_STR  = sys.argv[3]
    ALGO_PARAM_STR = sys.argv[4]
    #print("SIM_PARAM_STR ",SIM_PARAM_STR)
    #print("ALGO_PARAM_STR ",ALGO_PARAM_STR)

    ALGO_PARAMS    = parser.parse_params(ALGO_PARAM_STR)
    SIM_PARAMS     = parser.parse_params(SIM_PARAM_STR)

    print("SCENARIO:",SCENARIO)
    print("ALGO_NAME:",ALGO_NAME)
    print("ALGO_PARAMS:",ALGO_PARAMS)
    print("SIM_PARAMS:",SIM_PARAMS)
    
    output="../output"

    #mutation limits
    #limit_inf=1
    #limit_sup=420 
    #elite_size=2
	
    tic=timeit.default_timer()

    randomseed=15
    output=output+"/run"+str(randomseed)

    # Run the genetic algorithm
    sga.deap_sga_protein(SCENARIO, ALGO_PARAMS, SIM_PARAMS, output, randomseed).run()

    toc=timeit.default_timer()

    print("Execution Time=",toc-tic)


if __name__ == "__main__":
    main()
