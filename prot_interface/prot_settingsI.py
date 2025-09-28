#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 19/03/2025

@author: Rolando Armas
Yachay Tech University
Phage Therapy Group
"""

# Config Parameters
#######################################



#LIGAND_CHAIN = 'C'  # Specify the ligand chain (e.g., 'C')
#PARTNERS = 'A_C' # Specify partners (e.g., 'A_C')
START_PATTERN = "##### PAIRWISE SHORT-RANGE ENERGIES #####"

###Params rosetta-commons docker image
ROSETTA_BIN = "/usr/local/bin/"
CONFIG_PATH = "/app/scenarios/"
INTERF_EN = "interface_energy.default.linuxgccrelease -s"

###Params Cedia Cluster###
#ROSETTA_BIN = "/usr/local/bin/"
#CONFIG_PATH = "/root/Bio/scenarios/"
#INTERF_EN = "interface_energy.cxx11threadserialization.linuxgccrelease -s"


FACE1_FILE_NAME = "faceA.txt"
FACE2_FILE_NAME = "faceC.txt"
MSA_MATRIX = "MSA_matrix_old.tsv"

###Params Cedia Cluster###
INTERF_AN = "InterfaceAnalyzer.default.linuxgccrelease -s"

#INTERF_AN = "InterfaceAnalyzer.static.macosclangrelease -s"
SCORE_INDEXES = [4,7,5,11,23] #score positions for dG_separated, dSASA_int, dG_separated/dSASAx100, hbonds_int
FLAGS = "-compute_packstat true -tracer_data_print false -pack_input true " \
        "-pack_separated true -add_regular_scores_to_scorefile true " \
        "-atomic_burial_cutoff 0.01 -sasa_calculator_probe_radius 1.4 " \
        "-pose_metrics::interface_cutoff 8.0 -use_input_sc " \
        "-out:file:score_only "

