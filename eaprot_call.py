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
import prot_GA as sga
import prot_mob_GA as sga_mob
from prot_interface.logging_config import setup_logging
import logging
import os

# Initialize logging before anything else
setup_logging()
logger = logging.getLogger(__name__)

# Main function to run the genetic algorithm
def main():

    # command line arguments
    if len(sys.argv) != 8:
        print(len(sys.argv))
        logging.critical(f'usage: {sys.argv[0]} <scenario> <sea/moea> <sim params> <algo params> <fitness_idxs> <restore> <checkpoint>')
        sys.exit(-1)

    logging.info("Init Main")
    
    SCENARIO       = sys.argv[1] # Scenario Name
    # Given the scenario, set the pdb file and mutation limits
    
    ALGO_NAME      = sys.argv[2] # Single Objective (sea) / Multi Objective (moea)   
    SIM_PARAM_STR  = sys.argv[3]
    ALGO_PARAM_STR = sys.argv[4]
    FITNESS_IDXS   = sys.argv[5]
    CHECKPOINT_STR = str(sys.argv[6])
    CHECKPOINT     = True if CHECKPOINT_STR == "True" else False
    FREQ           = int(sys.argv[7])

    ALGO_PARAMS    = parser.parse_params(ALGO_PARAM_STR)
    SIM_PARAMS     = parser.parse_params(SIM_PARAM_STR)
    FITNESS_IDXS   = parser.parse_list(FITNESS_IDXS)

    logging.info("SCENARIO: %s",SCENARIO)
    logging.info("ALGO_NAME: %s",ALGO_NAME)
    logging.info("ALGO_PARAMS: %s",ALGO_PARAMS)
    logging.info("SIM_PARAMS: %s",SIM_PARAMS)
    logging.info("FITNESS_IDXS: %s",FITNESS_IDXS)
    
    output="../output"
	
    tic=timeit.default_timer()

    randomseed=15
    output=output+"/run"+str(randomseed)

    # Create output directory
    os.makedirs(output, exist_ok=True)
    # Create checkpoint.pkl if it isn't exist
    checkpoint_path = f'{output}/checkpoint.pkl'
    if not os.path.exists(checkpoint_path):
        with open(checkpoint_path, 'wb') as f:
            pass  
    

    # Run the genetic algorithm
    if(len(FITNESS_IDXS)>1):
        sga_mob.deap_mob_sga_protein(SCENARIO, ALGO_PARAMS, SIM_PARAMS, FITNESS_IDXS, output, randomseed).run(checkpoint=CHECKPOINT, freq=FREQ)
    else:
        sga.deap_sga_protein(SCENARIO, ALGO_PARAMS, SIM_PARAMS, FITNESS_IDXS, output, randomseed).run(checkpoint=CHECKPOINT, freq=FREQ)

    toc=timeit.default_timer()

    logging.info("Execution Time = %.2f seconds", toc - tic)


if __name__ == "__main__":
    main()
