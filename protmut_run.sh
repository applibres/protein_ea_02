#!/bin/bash

CLI="eaprot_call.py"

# command line arguments
#$./protmut_run.sh test04 sea pdbfile=protein01.pdb,partners=A_C,ligand_chain=C gen=5,popsize=5,obj=1,mutp=1 False 2

if [ $# -lt 6 ]; then
    echo "usage: ./protmut_run.sh <scenario> <algo> <sim_params> <algo_params> <checkpoint> <freq>"
    echo "for example:"
    echo "$./protmut_run.sh test04 sea pdbfile=protein01.pdb,partners=A_C,ligand_chain=C gen=5,popsize=5,obj=1 mutp=1 False,2"
    exit
fi

# Create checkpoint.pkl if they do not exist
if [ ! -f checkpoint.pkl ]; then
    touch checkpoint.pkl
fi

SCENARIO=$1
ALGO=$2
SIMU_PARAMS=$3
ALGO_PARAMS=$4
CHECKPOINT=$5
FREQ=$6

# run
#python eaprot_call.py test1 sea pdbfile=../scenarios/test01/protein01.pdb,partners=A_B,mut_init=235,mut_end=400 gen=5,popsize=5,obj=1,mutp=1 False 2

python $CLI $SCENARIO $ALGO $SIMU_PARAMS $ALGO_PARAMS $CHECKPOINT $FREQ
