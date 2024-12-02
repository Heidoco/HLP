#!/bin/bash

# Definir o valor de alpha
alpha=0.75

# Lista de revenues
revenues=(20 30 50)

# Lista de tamanhos de datasets
datasets=(100 125)

# Número de execuções para cada combinação
num_execucoes=1
echo "75L.txt" "0.75" "50"
./form3condAP "75L.txt" "0.75" "50"

# Loop através dos tamanhos de datasets
for size in "${datasets[@]}"; do
    # Loop para os sufixos L e T
    for suffix in T L; do
        # Nome do arquivo de dataset
        dataset="${size}${suffix}.txt"
        
        # Verificar se o arquivo existe antes de executar
        if [ -f "$dataset" ]; then
            # Loop através dos revenues
            for revenue in "${revenues[@]}"; do
                # Executar o programa 10 vezes para cada combinação
                for ((i=1; i<=num_execucoes; i++)); do
                    echo "Execução $i para dataset $dataset com revenue $revenue"
                    ./form3condAP "$dataset" "$alpha" "$revenue"
                done
            done
        else
            echo "Arquivo $dataset não encontrado."
        fi
    done
done


