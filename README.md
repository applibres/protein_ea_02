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

 

## Basic Use
You can download the project’s content as a zip file. Uncompressed it in the corresponding location that you decide to copy the project’s content. The directory contains the next structure: 

 
prot_interface: Folder with python scripts and classes which codes the interface with pyRosetta 

scenarios: Folder which contains the protein file in pdb format 

protmut_run.sh: Main shell script that calls to routines to evolutionary process 

eaprot_call.py: Main python script that parses the main evolutionary algorithm parameters and call to evolutionary python routine 

prot_GA.py: Evolutionary python algorithm 


To run the evolutionary algorithm, you can go inside of the directory were is located the project code and you can call to shell script file as is shown in the next line: 

./protmut_run.sh test04 sea pdbfile=protein01.pdb,partners=A_C,ligand_chain=C gen=4,popsize=8,obj=1,mutp=1 > ../output/logtest.txt


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

