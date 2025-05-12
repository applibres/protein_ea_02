# protein_ea_02

Evolutionary Algorithms for Protein Design 

Phage Therapy Group Yachay Tech 

## Introduction:
This project aims to implement evolutionary algorithms to study protein design.    

## Installation:
You need to install the software listed below: 

1) PyRosetta: PyRosetta is an interactive Python-based interface to the powerful Rosetta molecular modeling suite. It enables users to design their own custom molecular modeling algorithms using Rosetta sampling methods and energy functions. To install PyRosseta, you can follow the instruction given by the official web site. An alternative is to use rosetta image to run as a container with docker. You can use the next link to https://hub.docker.com/r/rosettacommons/rosetta or if you have docker, you can use the next docker command: pull rosettacommons/rosetta  

 

2) Distributed Evolutionary Algorithms in Python (DEAP): DEAP is a novel evolutionary computation framework for rapid prototyping and testing of ideas. It seeks to make algorithms explicit and data structures transparent. You can install deap directly using the next command: pip install deap. If you are running pyrosetta as a container, please don't forget to enter the container first and after that install deap.

3) Install the next python libraries: numpy, random, os, shutil 

4) Change the respective parameters values inside file prot_interface/prot_settingsI.py. The following parameters need to set according to your local repository and  rosetta-commons instalation:

ROSETTA_BIN = Main directory path wher is rosetta-commons binary files. ("/rosetta.binary.m1.release-371/main/source/bin/"

CONFIG_PATH = Main directory where are located the scenarios to test ("/scenarios/")

INTERF_EN = Name of interface_energy binary file ("interface_energy.static.macosclangrelease -s")

FACE1_FILE_NAME = Face A name file ("faceA.txt")

FACE2_FILE_NAME = Face C name file ("faceC.txt")

MSA_MATRIX = Name of MSA Matrix ("MSA_matrix.tsv") 

INTERF_AN = Name of InterfaceAnalyzer binary file ("InterfaceAnalyzer.static.macosclangrelease -s")

SCORE_INDEXES = [4,7,5,23] #score positions for dG_separated, dSASA_int, dG_separated/dSASAx100, hbonds_int

FLAGS = "-compute_packstat true -tracer_data_print false -pack_input true " \
        "-pack_separated true -add_regular_scores_to_scorefile true " \
        "-atomic_burial_cutoff 0.01 -sasa_calculator_probe_radius 1.4 " \
        "-pose_metrics::interface_cutoff 8.0 -use_input_sc " \
        "-out:file:score_only "


## Basic Use
You can download the project’s content as a zip file. Uncompressed it in the corresponding location that you decide to copy the project’s content. The directory contains the next structure: 

 
prot_interface: Folder with python scripts and classes which codes the interface with pyRosetta 

scenarios: Folder which contains the protein file in pdb format 

protmut_run.sh: Main shell script that calls to routines to evolutionary process 

eaprot_call.py: Main python script that parses the main evolutionary algorithm parameters and call to evolutionary python routine 

prot_GA.py: Evolutionary python algorithm 


To run the evolutionary algorithm, you can go inside of the directory were is located the project code and you can call to shell script file as is shown in the next line: 

## ./protmut_run.sh test04 sea pdbfile=protein01.pdb,partners=A_C,ligand_chain=C gen=4,popsize=8,obj=1,mutp=1 > ../output/logtest.txt


At following the parameters description: 

test04: the name of the scenario. Inside this folder there are the next files: faceA.txt, faceC.txt, MSA_Matrix.tsv, protein01.pdb, protein01.pdb.sc, protein01.pdb.txt   

sea: the name of the evolutionary algorithm (sea: single objective evolutionary algorithm) 

pdbfile: the name of the pdb file 

partners: The chains that conform the protein structure 

ligand_chain: The chain that will be used for evolution process

gen= number of generations  

popsize= population size  

obj= number of evaluation functions, for single objective evolutionary algorithm, obj=1 

mutp= mutation probability (not used)

