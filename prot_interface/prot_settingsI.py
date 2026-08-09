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
MSA_MATRIX = "MSA_matrix.tsv"

###Params Cedia Cluster###
INTERF_AN = "InterfaceAnalyzer.default.linuxgccrelease -s"

#INTERF_AN = "InterfaceAnalyzer.static.macosclangrelease -s"
#SCORE_INDEXES = [4,7,5,11,23] #score positions for dG_separated, dSASA_int, dG_separated/dSASAx100, hbonds_int
SCORE_INDEXES = [
    29, # packstat -> MAXIMIZAR (Calidad del empaquetamiento; valor ideal > 0.65)
    34, # sc_value -> MAXIMIZAR (Complementariedad de formas; valor ideal > 0.60) 
    0,  # total_score -> MINIMIZAR (Estabilidad global del complejo) [3, 4]
    9,  # delta_unsatHbonds -> MINIMIZAR (Penalización por polares no satisfechos; ideal cercano a 0) [5]
    16, # fa_rep -> MINIMIZAR (Repulsión estérica/choques; debe mantenerse bajo para ser físicamente posible) 
    30, # per_residue_energy_int -> MINIMIZAR (Energía promedio por residuo en la interfaz) 
    7,  # dSASA_int -> MAXIMIZAR (Área enterrada; valor ideal entre 1200-2000 A^2)
    5,  # dG_separated/dSASAx100 -> MINIMIZAR
    23] # hbonds_int
SCORE_OBJECTIVE = [1,1,-1,-1,-1,-1,1,-1,1]
FLAGS = "-compute_packstat true -tracer_data_print false -pack_input false " \
        "-pack_separated true -add_regular_scores_to_scorefile true " \
        "-atomic_burial_cutoff 0.01 -sasa_calculator_probe_radius 1.4 " \
        "-pose_metrics::interface_cutoff 8.0 -use_input_sc " \
        "-out:file:score_only "

