#!/bin/bash

# 1. Guardar el tiempo de inicio
START_TIME=$SECONDS

INIT=$1
END=$2

if [ $# -lt 2 ]; then
    echo "usage: ./replicates.sh <INIT> <END>"
    echo "example:"
    echo "./replicates.sh 10 30"
    exit 1
fi

for (( K=INIT; K<=END; K+=1 )); do 
  ./protmut_run.sh test06 sea pdbfile=9Q1V_prepared_clean_relaxed.pdb,partners=A_B,ligand_chain=B gen=20,popsize=25,mutp=0.3 fitness_idsx=5,7 checkpoint=False checks=2 mobj=True randomseed=$K
  if [ $? -ne 0 ]; then
    echo "-----------------------------Error in replicate $K---------------------------" 
    exit 1
  fi
  echo "====================Successfull execution for $K==================================="
  sleep 30
done

# 2. Calcular la diferencia
ELAPSED=$(( SECONDS - START_TIME ))

# 3. Formatear y mostrar (opcionalmente en min:seg)
echo "++++++++++++++++++++++Succesfull test++++++++++++++++++++++++++++"
echo "Tiempo total de ejecución: $((ELAPSED / 60))m $((ELAPSED % 60))s"
