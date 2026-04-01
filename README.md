# protein_ea_02

Evolutionary Algorithms for Protein Design

Phage Therapy Group Yachay Tech

## Introduction:
This project aims to implement evolutionary algorithms to study protein design.

## Installation:
You need to install the software listed below:

1) Create docker image from Dockerfile using

'''bash

docker build -t protein_ea_02 .

'''

2) Run the docker image

'''bash

docker run -it --rm protein_ea_02

'''

MultiObjective mode makes use of ESM2 model, so if your device has GPU, you can use the next command:

'''bash

docker run -it --rm --gpus all protein_ea_02

'''

3) Change the respective parameters values inside file prot_interface/prot_settingsI.py. The following parameters need to set according to your local repository and  rosetta-commons instalation:

ROSETTA_BIN = Main directory path wher is rosetta-commons binary files. ("/rosetta.binary.m1.release-371/main/source/bin/")

CONFIG_PATH = Main directory where are located the scenarios to test ("/scenarios/")

INTERF_EN = Name of interface_energy binary file ("interface_energy.static.macosclangrelease -s")

FACE1_FILE_NAME = Face A name file ("faceA.txt")

FACE2_FILE_NAME = Face C name file ("faceC.txt")

MSA_MATRIX = Name of MSA Matrix ("MSA_matrix.tsv")

INTERF_AN = Name of InterfaceAnalyzer binary file ("InterfaceAnalyzer.static.macosclangrelease -s")

SCORE_INDEXES = [4,7,5,11,23] #score positions for dG_separated, dSASA_int, dG_separated/dSASAx100, hbonds_int

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

prot_GA.py: Evolutionary python algorithm for single objective mode

prot_mob_GA.py: Evolutionary python algorithm for multi objective mode


To run the evolutionary algorithm, you can go inside of the directory were is located the project code and you can call to shell script file as is shown in the next line:

## ./protmut_run.sh test04 sea pdbfile=protein01.pdb,partners=A_C,ligand_chain=C gen=5,popsize=5,mutp=1 fitness_idxs=2,3 fitness_weights=-1,1 checkpoint=False checks=2 mobj=True randomseed=15


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

fitness_idxs= score positions for dG_separated, dSASA_int, dG_separated/dSASAx100, hbonds_int

fitness_weights= weights for each fitness function (-1 for minimize or 1 for maximize)

checkpoint= Restore from a previous checkpoint

checks= checkpoint frequency

mobj= algorith mode. If mobj=True uses fintess_idxs scores and ll score from esm2

random_seed= random seed