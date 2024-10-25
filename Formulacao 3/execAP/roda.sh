#!/bin/bash

# Array de valores para revenue
revenues=(20 30 50)

# Caminho do arquivo de entrada
input_files=("10L.txt" "10T.txt" "20L.txt" "20T.txt" "25L.txt" "25T.txt" "40L.txt" "40T.txt")

# Loop aninhado para gerar todas as combinações
for input_file in "${input_files[@]}"; do
    for revenue in "${revenues[@]}"; do
        # Executa o programa com os parâmetros atuais
        ./form3condAP $input_file 0.75 $revenue
    done
done

