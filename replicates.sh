#!/bin/bash

INIT=$1
END=$2

if [ $# -lt 2 ]; then
    echo "usage: ./exec_test.sh <INIT> <END>"
    echo "example:"
    echo "./exec_test.sh 10 30"
    exit 1
fi

for (( K=INIT; K<=END; K+=1 )); do 
  ./protmut_run.sh test04 sea pdbfile=protein01.pdb,partners=A_C,ligand_chain=C gen=20,popsize=25,mutp=0.3 fitness_idsx=2 False 2 False $K
  if [ $? -ne 0 ]; then
    echo "-----------------------------Error in replicate $K---------------------------" 
    exit 1
  fi
  echo "====================Successfull execution for $K==================================="
  sleep 30
done

echo "++++++++++++++++++++++Succesfull test++++++++++++++++++++++++++++"
