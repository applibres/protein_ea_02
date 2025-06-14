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
from prot_interface.logging_config import setup_logging
import logging

# Initialize logging before anything else
setup_logging()
logger = logging.getLogger(__name__)

# Main function to run the genetic algorithm
def main():

    # command line arguments
    if len(sys.argv) != 5:
        logging.critical('usage: ', sys.argv[0], '<scenario> <sea/moea> <sim params> <algo params>')
        sys.exit(-1)

    logging.info("Init Main")
    
    SCENARIO       = sys.argv[1] # Scenario Name
    # Given the scenario, set the pdb file and mutation limits
    
    ALGO_NAME      = sys.argv[2] # Single Objective (sea) / Multi Objective (moea)   
    SIM_PARAM_STR  = sys.argv[3]
    ALGO_PARAM_STR = sys.argv[4]

    ALGO_PARAMS    = parser.parse_params(ALGO_PARAM_STR)
    SIM_PARAMS     = parser.parse_params(SIM_PARAM_STR)

    logging.info("SCENARIO: %s",SCENARIO)
    logging.info("ALGO_NAME: %s",ALGO_NAME)
    logging.info("ALGO_PARAMS: %s",ALGO_PARAMS)
    logging.info("SIM_PARAMS: %s",SIM_PARAMS)
    
    output="../output"
	
    tic=timeit.default_timer()

    randomseed=15
    output=output+"/run"+str(randomseed)

    # Run the genetic algorithm
    sga.deap_sga_protein(SCENARIO, ALGO_PARAMS, SIM_PARAMS, output, randomseed).run()

    toc=timeit.default_timer()

    logging.info("Execution Time = %.2f seconds", toc - tic)


if __name__ == "__main__":
    main()
