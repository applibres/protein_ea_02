#!/bin/bash

# 1. Guardar el tiempo de inicio
START_TIME=$SECONDS

INIT=$1
END=$2

if [ $# -lt 2 ]; then
    echo "usage: ./exec_test.sh <INIT> <END>"
    echo "example:"
    echo "./exec_test.sh 10 30"
    exit 1
fi

for (( K=INIT; K<=END; K+=1 )); do 
  ./protmut_run.sh test05 sea pdbfile=6M0J_nowaters.pdb,partners=A_E,ligand_chain=E gen=20,popsize=25,mutp=0.3 fitness_idsx=2 False 2 True $K
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