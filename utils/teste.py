import json
import re

def parse_file(filename):
    with open(filename, 'r', encoding='utf-8') as file:
        content = file.read()

    # Dividir o conteúdo em blocos usando a linha de separação
    blocks = content.strip().split('------------------------------------------------------')

    data_list = []

    for block in blocks:
        block = block.strip()
        if not block:
            continue  # Pular blocos vazios

        lines = block.split('\n')

        data = {}

        # Extrair o nome do arquivo (primeira linha)
        data['arquivo'] = lines[0].strip()

        for line in lines[1:]:
            line = line.strip()
            if not line:
                continue

            # Usar expressão regular para capturar chave e valor
            match = re.match(r'([\w\s]+):\s*(.*)', line)
            if match:
                key = match.group(1).strip().lower().replace(' ', '_')
                value = match.group(2).strip()

                if key in ['alpha', 'receita', 'valor_inicial', 'valor_solucao', 'tempo']:
                    # Converter valores numéricos para o tipo correto
                    data[key] = float(value.replace(',', '.'))
                elif key == 'hubs_instalados':
                    # Verificar se há múltiplos valores
                    values = re.split(r'[\t\s]+', value)
                    data[key] = [int(v) for v in values if v]
                elif key == 'conexões_diretas_instaladas':
                    # Extrair as conexões diretas
                    connections = re.findall(r'\((\d+),(\d+)\)', value)
                    data[key] = [(int(a), int(b)) for a, b in connections]
            elif line.startswith('Conexões diretas instaladas:'):
                # Capturar conexões diretas que estão na mesma linha
                connections_line = line[len('Conexões diretas instaladas:'):].strip()
                connections = re.findall(r'\((\d+),(\d+)\)', connections_line)
                data['conexões_diretas_instaladas'] = [(int(a), int(b)) for a, b in connections]
            elif line.startswith('arcos instalados:'):
                # Pode ser ignorado se não houver dados após essa linha
                continue

        data_list.append(data)

    return data_list

def save_to_json(data_list, output_filename):
    with open(output_filename, 'w', encoding='utf-8') as outfile:
        json.dump(data_list, outfile, ensure_ascii=False, indent=4)

# Uso
data_list = parse_file('Resultados-configuracoes.txt')
save_to_json(data_list, 'dados.json')
