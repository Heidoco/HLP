#!/bin/bash

# Definir o valor de alpha
alpha=0.75

# Lista de revenues
revenues=(20 30 50)

# Lista de tamanhos de datasets
datasets=(100 125)

# Número de execuções para cada combinação
num_execucoes=10

# Loop através dos tamanhos de datasets
for size in "${datasets[@]}"; do
    echo "Execução"
    # Loop para os sufixos L e T
    for suffix in L T; do
        # Nome do arquivo de dataset
        dataset="${size}${suffix}.txt"
        echo "Execuçãoo"
        # Verificar se o arquivo existe antes de executar
        if [ -f "$dataset" ]; then
            # Loop através dos revenues
            for revenue in "${revenues[@]}"; do
                # Executar o programa 10 vezes para cada combinação
                for ((i=1; i<=num_execucoes; i++)); do
                    echo "Execução $i para dataset $dataset com revenue $revenue"
                    ./a.out "$dataset" "$alpha" "$revenue"
                done
            done
        else
            echo "Arquivo $dataset não encontrado."
        fi
    done
done
