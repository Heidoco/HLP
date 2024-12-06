#!/bin/bash

# Definir o valor de alpha
alpha=0.75

# Lista de revenues
revenues=(20 30 50)

# Lista de tamanhos de datasets
datasets=(25 40 50 75 100 125 150 200)

# Número de execuções para cada combinação
num_execucoes=10

# Loop através dos tamanhos de datasets
for size in "${datasets[@]}"; do
    for suffix in T L; do
        dataset="${size}${suffix}.txt"
        for revenue in "${revenues[@]}"; do
            for ((i=1; i<=num_execucoes; i++)); do
                echo "Execução 1 thread $i para dataset $dataset com revenue $revenue"
                    # Executar a.out e capturar tempo em tempoh.txt
                    time ./multi "$dataset" "$alpha" "$revenue";
                done
            done
    done
done
