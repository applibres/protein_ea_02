#!/bin/bash

CLI="eaprot_call.py"

# command line arguments
#$./protmut_run.sh test04 sea pdbfile=protein01.pdb,partners=A_C,ligand_chain=C gen=5,popsize=5,mutp=1 fitness_idxs=2,3 checkpoint=False checks=2 mobj=True randomseed=15

if [ $# -lt 8 ]; then
    echo "usage: ./protmut_run.sh <scenario> <algo> <sim_params> <algo_params> <checkpoint> <fitness_idxs> <freq> <multiobj> <randomseed>"
    echo "for example:"
    echo "$./protmut_run.sh test04 sea pdbfile=protein01.pdb,partners=A_C,ligand_chain=C gen=5,popsize=5,mutp=1 fitness_idxs=2,3 checkpoint=False checks=2 mobj=True randomseed=15"
    exit 1
fi

SCENARIO=$1
ALGO=$2
SIMU_PARAMS=$3
ALGO_PARAMS=$4
FITNSESSS_IDXS=$5
CHECKPOINT=$6
FREQ=$7
MOBJ=$8
RANDOMSEED=$9


if [ -f loggest.txt ] && [ "$CHECKPOINT" == "False" ]; then
    > loggest.txt
fi

# run
#python eaprot_call.py test1 sea pdbfile=../scenarios/test01/protein01.pdb,partners=A_B,mut_init=235,mut_end=400 gen=5,popsize=5,mutp=1 fitness_idxs=2,3 False 2 True 15

python $CLI $SCENARIO $ALGO $SIMU_PARAMS $ALGO_PARAMS $FITNSESSS_IDXS $CHECKPOINT $FREQ $MOBJ $RANDOMSEED