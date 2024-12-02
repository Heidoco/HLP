#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <list>
#include <math.h>
#include <time.h>
#include <stdexcept>
#include <algorithm>  
#include <random>     

using namespace std;

#define EPSILON 1e-07 /// Parâmetro usado para representar um número pequeno
#define M 10e10       /// Parâmetro usado para representar um número suficientemente grande

/// ==============================================
/// Recuperar índices das variáveis
/// ==============================================

inline int km(int n, int k, int m)
{
    if (k < m)
    {
        return (n - 1) * k + m - 1;
    }
    else if (k >= m)
    {
        return (n - 1) * k + m;
    }
}

class Heuristicas
{
public:
    // Parâmetros das instâncias
    int n;                    // Quantidade de nós
    float alpha;              // Fator de desconto no custo de transporte em um link entre hubs
    float alpha1;             // Coletar valor de alpha via linha de comando
    float receita;            // Coletar valor da receita via linha de comando
    vector<double> codx;      // Coordenada x dos nós (AP)
    vector<double> cody;      // Coordenada y dos nós (AP)
    vector<vector<double>> w; // declara um vetor de vetores para representar uma matriz com n linhas e n colunas - Quantidade de demanda a ser enviada entre os nós i e j
    vector<vector<double>> r; // Receita obtida por enviar uma unidade de demanda entre os nós i e j
    vector<vector<double>> c; // Custos por enviar uma unidade de demanda entre os nós i e j
    vector<double> s;         // Custos fixos de instalação de um hub
    vector<vector<double>> g; // Custos de operação nos links inter-hubs
    vector<vector<double>> q; // Custo operação conexao direta
    float soma;               // Soma para obter a média geral dos custos fixos de instalação hub - Dados AP
                              // Variáveis das instâncias
    vector<double> h;         // Vetor binário para os hubs ativos (tamanho n)
    vector<double> h_copia;
    vector<double> h_linha; // Vetor retornado pelas vizinhanças
    vector<double> h_star;

    vector<double> z; // Vetor binário para indicar os arcos entre hubs ativos (tamanho n*(n-1))
    vector<double> z_copia;
    vector<double> z_linha; // Vetor retornado pelas vizinhancas
    vector<double> z_star;

    vector<double> e; // Vetor binário para indicar as conexoes diretas entre nós (tamanho n*(n-1))
    vector<double> e_copia;
    vector<double> e_linha;
    vector<double> e_star;
    vector<double> e_zero;
    vector<double> e_registro;

    double valor_solucao = 0.0; // Valor da solução do problema
    double valor_solucao_linha; // Valor da função de avaliação retornado pelas vizinhanças
    double valor_solucao_star;
    double inicio_CPU;
    double fim_CPU;

public:
    Heuristicas() {} /// Construtor da classe
    void Read_main_arg(int argc, char *argv[]);
    void Read_data(char *arq);
    vector<double> sol_um_hub(); // Retorna um vetor de duas posições com o hub e o valor da solução com um único hub (se houver solução)
    double Shortest_path_algorithm();

    void vizinhanca1();           // insere um hub
    void vizinhanca2();           // remove um hub e os arcos incidentes sobre ele
    void vizinhanca3();           // adiciona um arco hub entre hubs ativos
    void vizinhanca4();           // remove um arco hub
    void vizinhanca5();           // troca um hub ativo por um não ativo e remove os arcos hubs incidentes sobre o hub desativado
    void vizinhanca6();           // insere um hub e o conecta aos demais hubs instalados

    void VND();
    void Perturbacao(int nivel);
    void perturbacao_v5();
    void SmartILS(int iter_max);
};

int main(int argc, char *argv[])
{
    // Cria uma instância da classe Heuristicas
    Heuristicas *ht = new Heuristicas();

    //========= Entrada dos dados =========//

    // Lê argumentos da linha de comando
    ht->Read_main_arg(argc, argv);
    // Lê os dados do arquivo especificado
    ht->Read_data(argv[1]);

    //========= Solução Inicial =========//

    // Inicia a contagem do tempo
    ht->inicio_CPU = clock();

    // Vetor para armazenar o resultado da solução inicial com um hub
    vector<double> vetor_um_hub;
    vetor_um_hub = ht->sol_um_hub();
    // Descomente as linhas abaixo para imprimir o hub ativo e o valor da solução
    // cout << "hub ativo: " << vetor_um_hub[0] << endl;
    // cout << "valor da solução com um hub ativo: " << vetor_um_hub[1] << endl;

    // Ativa o hub selecionado e define o valor da solução inicial
    ht->h[static_cast<int>(vetor_um_hub[0])] = 1;
    ht->valor_solucao = vetor_um_hub[1];

    // Armazena o valor inicial da solução
    double valor_inicial = ht->valor_solucao;

    //========= Heurística E-ILS-RVND =========//

    // Inicializa a semente para geração de números aleatórios
    srand(time(NULL));
    // Aplica o método SmartILS com número máximo de iterações igual a 4
    ht->SmartILS(4);

    // Finaliza a contagem do tempo
    ht->fim_CPU = clock();
    // Descomente a linha abaixo para imprimir o tempo de execução
    // printf("Tempo de execução = %10.2f segundos\n", (double)(ht->fim_CPU - ht->inicio_CPU) / CLOCKS_PER_SEC);
    // Descomente a linha abaixo para imprimir o valor da solução após o SmartILS
    // cout << "solução SmartILS: " << fixed << setprecision(2) << ht->valor_solucao_star << endl;

    //=========== Arquivo de Saída ===========//

    // Abre o arquivo para escrever os resultados da heurística
    ofstream arq_saida("Resultados-heuristica-vizinhancasall.txt", std::ofstream::app);
    if (!arq_saida)
    {
        cerr << "Erro ao abrir o arquivo de saída.\n";
        exit(0);
    }

    // Escreve os resultados no arquivo
    arq_saida << argv[1] << "\t" << ht->alpha << "\t" << ht->r[0][0] << "\t"
              << fixed << setprecision(2) << valor_inicial << "\t"
              << fixed << setprecision(2) << ht->valor_solucao_star << "\t"
              << fixed << setprecision(2) << (double)(ht->fim_CPU - ht->inicio_CPU) / CLOCKS_PER_SEC << endl;

    // Abre o arquivo para escrever a configuração da rede (opcional)
    ofstream arq_saida2("Resultados-configuracoes.txt", std::ofstream::app); // Ativar se quiser imprimir a configuração da rede
    if (!arq_saida2)
    {
        cerr << "Erro ao abrir o arquivo de configuração.\n";
        exit(0);
    }

    // Escreve detalhes da configuração no arquivo
    arq_saida2 << "\n\n" << argv[1]
               << "\nalpha: " << ht->alpha
               << "\nreceita: " << ht->r[0][0]
               << "\nvalor_inicial: " << fixed << setprecision(2) << valor_inicial
               << "\nvalor_solucao: " << fixed << setprecision(2) << ht->valor_solucao_star
               << "\ntempo: " << fixed << setprecision(2) << (double)(ht->fim_CPU - ht->inicio_CPU) / CLOCKS_PER_SEC << endl;

    // Escreve os hubs instalados
    arq_saida2 << "hubs instalados: ";
    for (int k = 0; k < ht->n; k++)
    {
        if (ht->h_star[k] > (1 - EPSILON))
        {
            arq_saida2 << k << "\t";
        }
    }

    // Escreve os arcos instalados entre hubs
    arq_saida2 << "\n" << "arcos instalados:";
    for (int k = 0; k < ht->n; k++)
    {
        for (int m = 0; m < ht->n; m++)
        {
            if (m != k)
            {
                if (ht->z_star[km(ht->n, k, m)] > (1 - EPSILON))
                {
                    arq_saida2 << "(" << k << "," << m << ")" << "\t";
                }
            }
        }
    }

    // Escreve as conexões diretas instaladas
    arq_saida2 << "\n" << "Conexões diretas instaladas: ";
    for (int i = 0; i < ht->n; i++)
    {
        for (int j = 0; j < ht->n; j++)
        {
            if (i != j)
            {
                if (ht->e_star[km(ht->n, i, j)] > (1 - EPSILON))
                {
                    arq_saida2 << "(" << i << "," << j << ")" << "\t";
                }
            }
        }
    }
    arq_saida2 << "\n" << "------------------------------------------------------";

    // Libera a memória alocada
    delete ht;

    return 0;
}


void Heuristicas::Read_main_arg(int argc, char *argv[])
{
    alpha1 = (argc > 2) ? atof(argv[2]) : 0.75; // valor de alpha via linha de comando
    receita = (argc > 3) ? atof(argv[3]) : 20;  // valor da receita via linha de comando
}

void Heuristicas::Read_data(char name[])
{
    // Abre o arquivo de entrada com os dados da instância
    ifstream arq(name);
    if (!arq)
    {
        cerr << "Erro ao abrir o arquivo de dados.\n";
        exit(0);
    }

    // Lê o número de nós (n) da instância
    arq >> n;

    // Inicializa os vetores e matrizes de acordo com o número de nós
    codx = vector<double>(n);                         // Coordenadas x dos nós
    cody = vector<double>(n);                         // Coordenadas y dos nós
    w = vector<vector<double>>(n, vector<double>(n)); // Matriz de demanda entre os nós i e j
    r = vector<vector<double>>(n, vector<double>(n)); // Matriz de receita entre os nós i e j
    c = vector<vector<double>>(n, vector<double>(n)); // Matriz de custos entre os nós i e j
    s = vector<double>(n);                            // Custos fixos de instalação dos hubs
    g = vector<vector<double>>(n, vector<double>(n)); // Custos de operação nos links inter-hubs
    q = vector<vector<double>>(n, vector<double>(n, 0.0)); // Custos de operação das conexões diretas
    h = vector<double>(n, 0);                         // Vetor binário indicando hubs ativos (inicialmente zeros)
    z = vector<double>(n * (n - 1), 0);               // Vetor binário para arcos entre hubs ativos
    e = vector<double>(n * (n - 1), 0);               // Vetor binário para conexões diretas entre nós
    e_registro = vector<double>(n * (n - 1), 0.0);    // Registro das conexões diretas
    e_linha = vector<double>(n * (n - 1), 0);         // Vetor auxiliar para conexões diretas
    e_zero = vector<double>(n * (n - 1), 0);          // Vetor de zeros para inicialização

    soma = 0; // Variável para armazenar a soma dos custos fixos de instalação dos hubs

    // Define o fator de desconto alpha (recebido via linha de comando)
    alpha = alpha1;

    // Inicializa a matriz de receita r[i][j] com o valor de 'receita' (recebido via linha de comando)
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            r[i][j] = receita; // Receita pode ser 20, 30 ou 50 (AP)
        }
    }

    // Lê as coordenadas x e y dos nós a partir do arquivo
    for (int i = 0; i < n; i++)
    {
        arq >> codx[i];
        arq >> cody[i];
    }

    // Calcula a matriz de custos c[i][j] entre os nós i e j
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            // Calcula a distância euclidiana multiplicada por um fator de escala (0.001)
            double dx = codx[i] - codx[j];
            double dy = cody[i] - cody[j];
            c[i][j] = 0.001 * sqrt(dx * dx + dy * dy);
        }
    }

    // Lê a matriz de demanda w[i][j] a partir do arquivo
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            arq >> w[i][j];
        }
    }

    // Lê os custos fixos de instalação dos hubs s[i] e calcula a soma total
    for (int i = 0; i < n; i++)
    {
        arq >> s[i];
        s[i] = s[i] / 10; // Ajuste do custo (dividido por 10 conforme AP)
        soma += s[i];     // Acumula para calcular a média
    }

    // Calcula os custos de operação g[i][j] nos links inter-hubs e os custos q[i][j] das conexões diretas
    double media_custos_instalacao = soma / n;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            // Custo de operação nos links inter-hubs
            g[i][j] = 0.1 * media_custos_instalacao;
            // Custo de operação das conexões diretas (20% de g[i][j])
            q[i][j] = 0.20 * g[i][j];
        }
    }

    // Fecha o arquivo após a leitura
    arq.close();
}


vector<double> Heuristicas::sol_um_hub()
{
    // Solução inicial considerando um único hub

    double phi = 0.0;        // Valor máximo da função objetivo encontrado
    double phi_aux = 0.0;    // Valor auxiliar da função objetivo para comparações
    int hub_ativo = -1;      // Índice do hub ativo (inicialmente nenhum)
    vector<vector<double>> l(n, vector<double>(n));    // Matriz de lucro usando o hub
    vector<vector<double>> le(n, vector<double>(n));   // Matriz de lucro usando conexão direta
    vector<double> vet_sol_um_hub(2);                  // Vetor para armazenar o hub ativo e o valor máximo de phi

    // Percorre todos os possíveis hubs para encontrar aquele que maximiza o lucro
    for (int k = 0; k < n; k++)
    {
        phi_aux = -s[k];     // Subtrai o custo fixo de instalação do hub k
        e_linha = e_zero;    // Reinicia o vetor de conexões diretas

        // Calcula o lucro para cada par de nós (i, j)
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                // Lucro ao passar pelo hub k
                l[i][j] = r[i][j] * w[i][j] - ( (c[i][k] + c[k][j]) * w[i][j] );

                // Lucro por conexão direta entre i e j
                le[i][j] = r[i][j] * w[i][j] - c[i][j] * w[i][j] - q[i][j];

                // Verifica qual opção é mais lucrativa e se há lucro positivo
                if (l[i][j] > EPSILON && l[i][j] > le[i][j])
                {
                    phi_aux += l[i][j];  // Adiciona o lucro via hub k
                }
                else if (le[i][j] > EPSILON && le[i][j] > l[i][j])
                {
                    phi_aux += le[i][j];                  // Adiciona o lucro via conexão direta
                    e_linha[km(n, i, j)] = 1;             // Marca a conexão direta como utilizada
                }
            }
        }

        // Verifica se o lucro total com o hub k é maior que o melhor encontrado até agora
        if (phi_aux >= phi)
        {
            e = e_linha;     // Atualiza o vetor de conexões diretas com a melhor configuração
            phi = phi_aux;   // Atualiza o valor máximo de phi
            hub_ativo = k;   // Registra o hub ativo correspondente
        }
    }

    // Armazena o hub ativo e o valor máximo de phi no vetor de retorno
    vet_sol_um_hub[0] = hub_ativo;
    vet_sol_um_hub[1] = phi;

    // Ativa o hub selecionado no vetor h
    if (hub_ativo >= EPSILON)
        h[hub_ativo] = 1;

    return vet_sol_um_hub;
}


double Heuristicas::Shortest_path_algorithm()
{
    // ==============================================
    // Algoritmo para calcular o caminho mais curto e o valor da solução
    // ==============================================

    // Matrizes para armazenar os caminhos mais curtos
    vector<vector<double>> hub_net_sp_matrix(n, vector<double>(n, M)); // Caminhos mais curtos na rede de hubs
    vector<vector<double>> sp_matrix(n, vector<double>(n, M));         // Caminhos mais curtos para todos os nós

    double sol_value = 0.0;    // Valor da solução (lucro total)
    list<int> hub_list;        // Lista dos hubs ativos
    vector<vector<double>> l(n, vector<double>(n));   // Lucro via hubs
    vector<vector<double>> le(n, vector<double>(n));  // Lucro via conexões diretas

    // ==============================================
    // 1. Monta a lista de hubs ativos e subtrai os custos de instalação
    // ==============================================

    for (int k = 0; k < n; k++)
    {
        if (h[k] >= (1 - EPSILON))
        {
            hub_list.push_back(k);    // Adiciona o hub ativo à lista
            sol_value -= s[k];        // Subtrai o custo fixo de instalação do hub k
        }
    }

    // ==============================================
    // 2. Calcula os caminhos mais curtos entre os hubs na rede de hubs
    // ==============================================

    // Inicializa os custos diretos entre hubs ativos
    for (auto k = hub_list.begin(); k != hub_list.end(); k++)
    {
        for (auto m = hub_list.begin(); m != hub_list.end(); m++)
        {
            if (*k != *m)
            {
                if (z[km(n, *k, *m)] >= (1 - EPSILON))
                {
                    // Custo direto entre hubs conectados
                    hub_net_sp_matrix[*k][*m] = alpha * c[*k][*m];
                    sol_value -= g[*k][*m];    // Subtrai o custo de operação do link inter-hub
                }
                else
                {
                    // Não há conexão direta; mantém o custo alto (M)
                    hub_net_sp_matrix[*k][*m] = M;
                }
            }
        }
    }

    // Aplica o algoritmo de Floyd-Warshall para encontrar os caminhos mais curtos entre os hubs
    for (auto k = hub_list.begin(); k != hub_list.end(); k++)
    {
        for (auto m = hub_list.begin(); m != hub_list.end(); m++)
        {
            for (auto l = hub_list.begin(); l != hub_list.end(); l++)
            {
                if ((*m != *l) && (*m != *k) && (*l != *k))
                {
                    double cost_aux = hub_net_sp_matrix[*m][*k] + hub_net_sp_matrix[*k][*l];
                    if (cost_aux < hub_net_sp_matrix[*m][*l])
                    {
                        hub_net_sp_matrix[*m][*l] = cost_aux;  // Atualiza o caminho mais curto
                    }
                }
            }
        }
    }

    // ==============================================
    // 3. Calcula os caminhos mais curtos para todos os pares de nós (i, j)
    // ==============================================

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            sp_matrix[i][j] = M;  // Inicializa com um valor grande

            // Considera todos os hubs ativos como possíveis intermediários
            for (auto k = hub_list.begin(); k != hub_list.end(); k++)
            {
                // Custo passando por um único hub (*k)
                sp_matrix[i][j] = min(c[i][*k] + c[*k][j], sp_matrix[i][j]);

                for (auto m = hub_list.begin(); m != hub_list.end(); m++)
                {
                    if (*k != *m)
                    {
                        // Custo passando por dois hubs (*k e *m)
                        double cost_aux = c[i][*k] + hub_net_sp_matrix[*k][*m] + c[*m][j];
                        if (sp_matrix[i][j] > cost_aux)
                        {
                            sp_matrix[i][j] = cost_aux;  // Atualiza o caminho mais curto
                        }
                    }
                }
            }

            // ==============================================
            // 4. Calcula o lucro para cada par de nós (i, j)
            // ==============================================

            // Lucro usando a rede de hubs
            l[i][j] = (r[i][j] - sp_matrix[i][j]) * w[i][j];

            // Lucro usando conexão direta
            le[i][j] = (r[i][j] - c[i][j]) * w[i][j] - q[i][j];

            // Seleciona a opção com maior lucro positivo
            if (l[i][j] > EPSILON && l[i][j] > le[i][j])
            {
                sol_value += l[i][j];  // Adiciona o lucro via hubs
            }
            else if (le[i][j] > EPSILON && le[i][j] > l[i][j])
            {
                sol_value += le[i][j];                 // Adiciona o lucro via conexão direta
                e_linha[km(n, i, j)] = 1;              // Marca a conexão direta como utilizada
            }
        }
    }

    // ==============================================
    // 5. Limpeza e retorno do valor da solução
    // ==============================================

    hub_net_sp_matrix.clear();  // Limpa a matriz de caminhos entre hubs
    hub_list.clear();           // Limpa a lista de hubs ativos

    return sol_value;           // Retorna o valor total da solução
}



// Estruturas de vizinhança

// Vizinhança 1: Inserção de um hub
void Heuristicas::vizinhanca1()
{
    

    vector<int> hubs_inativos;  // Vetor para armazenar os índices dos hubs inativos
    double f_avaliacao = 0.0;   // Variável para armazenar o valor da função objetivo na iteração atual
    int u;                      // Tamanho da vizinhança que será explorada
    int posicao;                // Posição aleatória selecionada no vetor de hubs inativos

    // Salva o estado atual dos hubs e arcos
    h_copia = h;                // Copia do vetor de hubs ativos
    h_linha = h;                // Inicializa o vetor h_linha com o estado atual
    valor_solucao_linha = valor_solucao; // Valor atual da solução
    z_linha = z;                // Copia do vetor de arcos entre hubs

    // Coleta os hubs inativos
    for (int k = 0; k < n; k++)
    {
        if (h[k] < EPSILON)
        {
            hubs_inativos.push_back(k); // Armazena o índice do hub inativo
        }
    }

    // Verifica se há hubs inativos para explorar
    if (!hubs_inativos.empty())
    {
        u = hubs_inativos.size(); // Define o tamanho da vizinhança (pode ser ajustado)

        // Explora a vizinhança adicionando hubs inativos
        for (int j = 0; j < u; j++)
        {
            // Seleciona aleatoriamente um hub inativo
            posicao = rand() % hubs_inativos.size();
            int hub_selecionado = hubs_inativos[posicao];

            // Ativa o hub selecionado
            h[hub_selecionado] = 1;

            // Avalia a nova solução
            f_avaliacao = Shortest_path_algorithm();

            // Verifica se a nova solução é melhor que a melhor encontrada até agora
            if (f_avaliacao > valor_solucao_linha)
            {
                h_linha = h;                   // Atualiza o melhor vetor de hubs
                valor_solucao_linha = f_avaliacao; // Atualiza o melhor valor da solução
            }

            // Remove o hub selecionado da lista de inativos
            hubs_inativos.erase(hubs_inativos.begin() + posicao);

            // Restaura o estado original dos hubs para a próxima iteração
            h = h_copia;
        }
    }

    // Limpa o vetor de hubs inativos
    hubs_inativos.clear();
}


void Heuristicas::vizinhanca6()
{

    vector<double> hubs_inativos;
    vector<double> hubs_ativos;
    double f_avaliacao = 0.0;
    int u; // tamanho da vizinhanca que será explorada
    int posicao;

    h_copia = h;
    h_linha = h;
    valor_solucao_linha = valor_solucao;
    z_copia = z;
    z_linha = z;

    // cout << "sol_linha inicial: " << valor_solucao_linha << endl;

    for (int k = 0; k < n; k++)
    { // coleta os hubs inativos
        if (h[k] < EPSILON)
        {
            hubs_inativos.push_back(k);
            // cout << "posicao inativos: " << k << endl;
        }
        else
        {
            hubs_ativos.push_back(k);
        }
    }

    if (!hubs_inativos.empty())
    {
        u = 1 * (hubs_inativos.size());
        // cout << "valor de m: " << m << endl;
        for (int j = 0; j < u; j++)
        {
            if (j == 0)
            {
                posicao = rand() % hubs_inativos.size();
                // cout << "posicao 0: " << hubs_inativos[posicao] << endl;
                h[hubs_inativos[posicao]] = 1;
                for (int k = 0; k < hubs_ativos.size(); k++)
                {
                    z[km(n, hubs_inativos[posicao], hubs_ativos[k])] = 1;
                    z[km(n, hubs_ativos[k], hubs_inativos[posicao])] = 1;
                }
                valor_solucao_linha = Shortest_path_algorithm();
                h_linha = h;
                z_linha = z;
                hubs_inativos.erase(hubs_inativos.begin() + posicao);
                h = h_copia;
                z = z_copia;
            }
            else
            {
                posicao = rand() % hubs_inativos.size();
                //  cout << "posicao depois: " << hubs_inativos[posicao] << endl;
                h[hubs_inativos[posicao]] = 1;
                for (int k = 0; k < hubs_ativos.size(); k++)
                {
                    z[km(n, hubs_inativos[posicao], hubs_ativos[k])] = 1;
                    z[km(n, hubs_ativos[k], hubs_inativos[posicao])] = 1;
                }
                f_avaliacao = Shortest_path_algorithm();
                // cout << "f_avaliacao: " << f_avaliacao << endl;
                // cout << "sol_linha: " << valor_solucao_linha << endl;
                if (f_avaliacao > valor_solucao_linha)
                {
                    h_linha = h;
                    z_linha = z;
                    valor_solucao_linha = f_avaliacao;
                }
                hubs_inativos.erase(hubs_inativos.begin() + posicao);
                h = h_copia;
                z = z_copia;
            }
        }
    }
    hubs_inativos.clear();
    hubs_ativos.clear();
}

void Heuristicas::vizinhanca2()
{
    // ==============================================
    // Vizinhança 2: Remoção de um hub ativo
    // ==============================================

    vector<int> hubs_ativos;    // Vetor para armazenar os índices dos hubs ativos
    vector<int> arc_hub1;       // Vetor para armazenar o primeiro nó dos arcos entre hubs ativos
    vector<int> arc_hub2;       // Vetor para armazenar o segundo nó dos arcos entre hubs ativos
    double f_avaliacao = 0.0;   // Variável para armazenar o valor da função objetivo na iteração atual
    int u;                      // Tamanho da vizinhança que será explorada
    int posicao;                // Posição aleatória selecionada no vetor de hubs ativos

    // Salva o estado atual dos hubs e arcos
    h_copia = h;                // Cópia do vetor de hubs ativos
    h_linha = h;                // Inicializa o vetor h_linha com o estado atual
    z_copia = z;                // Cópia do vetor de arcos entre hubs
    z_linha = z;                // Inicializa z_linha com o estado atual
    valor_solucao_linha = valor_solucao; // Valor atual da solução

    // ==============================================
    // Coleta dos hubs ativos
    // ==============================================
    for (int k = 0; k < n; k++)
    {
        if (h[k] > (1 - EPSILON))
        {
            hubs_ativos.push_back(k); // Armazena o índice do hub ativo
        }
    }

    // ==============================================
    // Coleta dos arcos entre hubs ativos
    // ==============================================
    for (int k = 0; k < n; k++)
    {
        for (int m = 0; m < n; m++)
        {
            if (m != k)
            {
                if (z[km(n, k, m)] > (1 - EPSILON))
                {
                    arc_hub1.push_back(k);
                    arc_hub2.push_back(m);
                }
            }
        }
    }

    // ==============================================
    // Exploração da vizinhança: remoção de um hub ativo
    // ==============================================
    if (hubs_ativos.size() > 1)
    {
        u = hubs_ativos.size(); // Define o tamanho da vizinhança como o número de hubs ativos

        for (int j = 0; j < u; j++)
        {
            // Seleciona aleatoriamente um hub ativo para desativar
            posicao = rand() % hubs_ativos.size();
            int hub_selecionado = hubs_ativos[posicao];

            // Desativa o hub selecionado
            h[hub_selecionado] = 0;

            // Remove os arcos incidentes ao hub desativado
            for (size_t i = 0; i < arc_hub1.size(); i++)
            {
                if (arc_hub1[i] == hub_selecionado || arc_hub2[i] == hub_selecionado)
                {
                    z[km(n, arc_hub1[i], arc_hub2[i])] = 0;
                }
            }

            // Avalia a nova solução
            f_avaliacao = Shortest_path_algorithm();

            // Verifica se a nova solução é melhor que a melhor encontrada até agora
            if (f_avaliacao > valor_solucao_linha)
            {
                h_linha = h;                   // Atualiza o melhor vetor de hubs
                z_linha = z;                   // Atualiza o melhor vetor de arcos
                valor_solucao_linha = f_avaliacao; // Atualiza o melhor valor da solução
            }

            // Remove o hub selecionado da lista de hubs ativos
            hubs_ativos.erase(hubs_ativos.begin() + posicao);

            // Restaura o estado original dos hubs e arcos para a próxima iteração
            h = h_copia;
            z = z_copia;
        }
    }

    // Limpa os vetores auxiliares
    hubs_ativos.clear();
    arc_hub1.clear();
    arc_hub2.clear();
}


/**
 * @brief Vizinhança 3: Adiciona um arco entre hubs ativos
 *
 * Esta função implementa a terceira estrutura de vizinhança do algoritmo heurístico,
 * onde se adiciona um arco entre dois hubs ativos que ainda não possuem conexão direta.
 * O objetivo é explorar soluções vizinhas à atual, adicionando arcos entre hubs ativos
 * e avaliando se a nova configuração melhora o valor da função objetivo.
 */
void Heuristicas::vizinhanca3()
{
    vector<int> hubs_ativos;  // Vetor para armazenar os índices dos hubs ativos
    vector<int> arc_hub1;     // Vetor para armazenar o primeiro nó dos arcos possíveis entre hubs
    vector<int> arc_hub2;     // Vetor para armazenar o segundo nó dos arcos possíveis entre hubs
    double f_avaliacao = 0.0; // Valor da função objetivo na iteração atual
    int u;                    // Tamanho da vizinhança que será explorada
    int posicao;              // Posição aleatória selecionada no vetor de arcos possíveis

    // Salva o estado atual dos hubs e arcos
    h_linha = h;
    z_copia = z;
    z_linha = z;
    valor_solucao_linha = valor_solucao;

    // Coleta os hubs ativos
    for (int k = 0; k < n; k++)
    {
        if (h[k] > (1 - EPSILON))
        {
            hubs_ativos.push_back(k);
        }
    }

    // Verifica se há pelo menos dois hubs ativos para conectar
    if (hubs_ativos.size() > 1)
    {
        // Monta o conjunto de arcos possíveis entre hubs ativos que ainda não estão conectados
        for (size_t i = 0; i < hubs_ativos.size(); i++)
        {
            for (size_t j = 0; j < hubs_ativos.size(); j++)
            {
                if ((i != j) && (z[km(n, hubs_ativos[i], hubs_ativos[j])] < EPSILON))
                {
                    arc_hub1.push_back(hubs_ativos[i]);
                    arc_hub2.push_back(hubs_ativos[j]);
                }
            }
        }

        // Define o tamanho da vizinhança a ser explorada
        u = arc_hub1.size();

        // Explora a vizinhança adicionando arcos entre hubs ativos
        for (int j = 0; j < u; j++)
        {
            // Seleciona aleatoriamente um arco possível
            posicao = rand() % arc_hub1.size();
            int hub_origem = arc_hub1[posicao];
            int hub_destino = arc_hub2[posicao];

            // Adiciona o arco entre os hubs selecionados
            z[km(n, hub_origem, hub_destino)] = 1;

            // Avalia a nova solução
            f_avaliacao = Shortest_path_algorithm();

            // Verifica se a nova solução é melhor que a melhor encontrada até agora
            if (f_avaliacao > valor_solucao_linha)
            {
                z_linha = z;
                valor_solucao_linha = f_avaliacao;
            }

            // Remove o arco selecionado dos vetores de arcos possíveis
            arc_hub1.erase(arc_hub1.begin() + posicao);
            arc_hub2.erase(arc_hub2.begin() + posicao);

            // Restaura o estado original dos arcos para a próxima iteração
            z = z_copia;
        }
    }

    // Limpa os vetores auxiliares
    hubs_ativos.clear();
    arc_hub1.clear();
    arc_hub2.clear();
}


/**
 * @brief Vizinhança 4: Remoção de um arco entre hubs ativos
 *
 * Esta função implementa a quarta estrutura de vizinhança do algoritmo heurístico,
 * onde se remove um arco existente entre dois hubs ativos.
 * O objetivo é explorar soluções vizinhas à atual, removendo arcos entre hubs ativos
 * e avaliando se a nova configuração melhora o valor da função objetivo.
 */
void Heuristicas::vizinhanca4()
{
    vector<int> arc_hub1;     // Vetor para armazenar o primeiro nó dos arcos ativos entre hubs
    vector<int> arc_hub2;     // Vetor para armazenar o segundo nó dos arcos ativos entre hubs
    double f_avaliacao = 0.0; // Valor da função objetivo na iteração atual
    int u;                    // Tamanho da vizinhança que será explorada
    int posicao;              // Posição aleatória selecionada no vetor de arcos ativos

    // Salva o estado atual dos hubs e arcos
    h_linha = h;
    z_copia = z;
    z_linha = z;
    valor_solucao_linha = valor_solucao;

    // Coleta os arcos ativos entre hubs
    for (int k = 0; k < n; k++)
    {
        for (int m = 0; m < n; m++)
        {
            if (m != k)
            {
                if (z[km(n, k, m)] > (1 - EPSILON))
                {
                    arc_hub1.push_back(k);
                    arc_hub2.push_back(m);
                }
            }
        }
    }

    // Verifica se há arcos ativos para explorar
    if (!arc_hub1.empty())
    {
        u = arc_hub1.size(); // Define o tamanho da vizinhança como o número de arcos ativos

        // Explora a vizinhança removendo arcos entre hubs ativos
        for (int j = 0; j < u; j++)
        {
            // Seleciona aleatoriamente um arco ativo para remover
            posicao = rand() % arc_hub1.size();
            int hub_origem = arc_hub1[posicao];
            int hub_destino = arc_hub2[posicao];

            // Remove o arco entre os hubs selecionados
            z[km(n, hub_origem, hub_destino)] = 0;

            // Avalia a nova solução
            f_avaliacao = Shortest_path_algorithm();

            // Verifica se a nova solução é melhor que a melhor encontrada até agora
            if (f_avaliacao > valor_solucao_linha)
            {
                z_linha = z;
                valor_solucao_linha = f_avaliacao;
            }

            // Remove o arco selecionado dos vetores de arcos ativos
            arc_hub1.erase(arc_hub1.begin() + posicao);
            arc_hub2.erase(arc_hub2.begin() + posicao);

            // Restaura o estado original dos arcos para a próxima iteração
            z = z_copia;
        }
    }

    // Limpa os vetores auxiliares
    arc_hub1.clear();
    arc_hub2.clear();
}


/**
 * @brief Vizinhança 5: Troca de um hub ativo por um hub inativo
 *
 * Esta função implementa a quinta estrutura de vizinhança do algoritmo heurístico,
 * onde um hub ativo é desativado e um hub inativo é ativado em seu lugar.
 * O objetivo é explorar soluções vizinhas à atual, realizando trocas de hubs
 * e avaliando se a nova configuração melhora o valor da função objetivo.
 */
void Heuristicas::vizinhanca5()
{
    vector<int> hubs_ativos;          // Vetor para armazenar os índices dos hubs ativos
    vector<int> hubs_inativos;        // Vetor para armazenar os índices dos hubs inativos
    vector<int> hubs_inativos_copia;  // Cópia dos hubs inativos para restauração
    vector<int> arc_hub1;             // Vetor para armazenar o primeiro nó dos arcos ativos entre hubs
    vector<int> arc_hub2;             // Vetor para armazenar o segundo nó dos arcos ativos entre hubs
    double f_avaliacao = 0.0;         // Valor da função objetivo na iteração atual
    int posicao;                      // Posição aleatória selecionada nos vetores de hubs

    // Salva o estado atual dos hubs e arcos
    h_copia = h;
    h_linha = h;
    z_copia = z;
    z_linha = z;
    valor_solucao_linha = valor_solucao;

    // Coleta os hubs ativos e inativos
    for (int k = 0; k < n; k++)
    {
        if (h[k] > (1 - EPSILON))
        {
            hubs_ativos.push_back(k);
        }
        else
        {
            hubs_inativos.push_back(k);
        }
    }

    // Cópia dos hubs inativos para restauração posterior
    hubs_inativos_copia = hubs_inativos;

    // Coleta os arcos ativos entre hubs
    for (int k = 0; k < n; k++)
    {
        for (int m = 0; m < n; m++)
        {
            if (m != k && z[km(n, k, m)] > (1 - EPSILON))
            {
                arc_hub1.push_back(k);
                arc_hub2.push_back(m);
            }
        }
    }

    // Se não houver hubs ativos, ativa um hub aleatoriamente e reinicia a vizinhança
    if (hubs_ativos.empty())
    {
        posicao = rand() % hubs_inativos.size();
        h[hubs_inativos[posicao]] = 1;
        valor_solucao = Shortest_path_algorithm();
        vizinhanca5(); // Chama novamente a função após ativar um hub
    }
    else
    {
        // Explora a vizinhança trocando um hub ativo por um inativo
        for (size_t i = 0; i < hubs_ativos.size(); i++)
        {
            // Define o número de trocas a serem exploradas
            int u = hubs_inativos.size();

            for (int j = 0; j < u; j++)
            {
                // Desativa o hub ativo atual
                h[hubs_ativos[i]] = 0;

                // Seleciona aleatoriamente um hub inativo para ativar
                posicao = rand() % hubs_inativos.size();
                int hub_inativo_selecionado = hubs_inativos[posicao];
                h[hub_inativo_selecionado] = 1;

                // Remove os arcos incidentes ao hub desativado
                for (size_t k = 0; k < arc_hub1.size(); k++)
                {
                    if (arc_hub1[k] == hubs_ativos[i] || arc_hub2[k] == hubs_ativos[i])
                    {
                        z[km(n, arc_hub1[k], arc_hub2[k])] = 0;
                    }
                }

                // Avalia a nova solução
                f_avaliacao = Shortest_path_algorithm();

                // Verifica se a nova solução é melhor que a melhor encontrada até agora
                if (f_avaliacao > valor_solucao_linha)
                {
                    h_linha = h;
                    z_linha = z;
                    valor_solucao_linha = f_avaliacao;
                }

                // Remove o hub inativo selecionado da lista para não selecioná-lo novamente
                hubs_inativos.erase(hubs_inativos.begin() + posicao);

                // Restaura o estado original dos hubs e arcos para a próxima iteração
                h = h_copia;
                z = z_copia;
            }

            // Restaura a lista de hubs inativos para a próxima iteração do hub ativo
            hubs_inativos = hubs_inativos_copia;
        }
    }

    // Limpa os vetores auxiliares
    hubs_ativos.clear();
    hubs_inativos.clear();
    arc_hub1.clear();
    arc_hub2.clear();
}


/**
 * @brief Perturbação V5: Troca aleatória de um hub ativo por um hub inativo
 *
 * Esta função implementa uma perturbação baseada na vizinhança 5, onde
 * um hub ativo é desativado e um hub inativo é ativado em seu lugar.
 * É usada para perturbar a solução atual, explorando novas regiões do espaço de soluções.
 */
void Heuristicas::perturbacao_v5()
{
    vector<int> hubs_ativos;          // Vetor para armazenar os índices dos hubs ativos
    vector<int> hubs_inativos;        // Vetor para armazenar os índices dos hubs inativos
    vector<int> arc_hub1;             // Vetor para armazenar o primeiro nó dos arcos ativos entre hubs
    vector<int> arc_hub2;             // Vetor para armazenar o segundo nó dos arcos ativos entre hubs
    int posicao;                      // Posição aleatória selecionada nos vetores de hubs
    int posicao2;                     // Posição aleatória para o segundo hub

    // Coleta dos hubs ativos e inativos
    for (int k = 0; k < n; k++)
    {
        if (h[k] > (1 - EPSILON))
        {
            hubs_ativos.push_back(k);
        }
        else
        {
            hubs_inativos.push_back(k);
        }
    }

    // Coleta dos arcos ativos entre hubs
    for (int k = 0; k < n; k++)
    {
        for (int m = 0; m < n; m++)
        {
            if (m != k && z[km(n, k, m)] > (1 - EPSILON))
            {
                arc_hub1.push_back(k);
                arc_hub2.push_back(m);
            }
        }
    }

    // Se não houver hubs ativos, ativa um hub aleatoriamente e reinicia a perturbação
    if (hubs_ativos.empty())
    {
        posicao = rand() % hubs_inativos.size();
        h[hubs_inativos[posicao]] = 1;
        valor_solucao = Shortest_path_algorithm();
        perturbacao_v5(); // Chama novamente a função após ativar um hub
    }
    else
    {
        // Seleciona aleatoriamente um hub ativo para desativar
        posicao = rand() % hubs_ativos.size();
        int hub_ativo_selecionado = hubs_ativos[posicao];
        h[hub_ativo_selecionado] = 0;

        // Seleciona aleatoriamente um hub inativo para ativar
        posicao2 = rand() % hubs_inativos.size();
        int hub_inativo_selecionado = hubs_inativos[posicao2];
        h[hub_inativo_selecionado] = 1;

        // Remove os arcos incidentes ao hub desativado
        for (size_t k = 0; k < arc_hub1.size(); k++)
        {
            if (arc_hub1[k] == hub_ativo_selecionado || arc_hub2[k] == hub_ativo_selecionado)
            {
                z[km(n, arc_hub1[k], arc_hub2[k])] = 0;
            }
        }

        // Recalcula o valor da solução
        valor_solucao = Shortest_path_algorithm();
    }

    // Limpa os vetores auxiliares
    hubs_ativos.clear();
    hubs_inativos.clear();
    arc_hub1.clear();
    arc_hub2.clear();
}


/**
 * @brief Implementação do VND (Variable Neighborhood Descent)
 *
 * Esta função implementa o algoritmo VND, explorando diferentes estruturas de vizinhança
 * para melhorar a solução atual. As vizinhanças são exploradas em ordem aleatória,
 * porém o embaralhamento ocorre apenas no início.
 */
void Heuristicas::VND()
{
    // Vetor com os índices das vizinhanças
    vector<int> vetor = {0, 1, 2, 3, 4, 5};

    // Embaralha aleatoriamente o vetor de vizinhanças
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(vetor.begin(), vetor.end(), g);

    int i = 0;

    // Explora as vizinhanças até que todas tenham sido exploradas sem melhoria
    while (i < 6)
    {
        // Seleciona a vizinhança atual baseada no vetor embaralhado
        switch (vetor[i])
        {
        case 0:
            vizinhanca1();
            break;
        case 1:
            vizinhanca2();
            break;
        case 2:
            vizinhanca3();
            break;
        case 3:
            vizinhanca4();
            break;
        case 4:
            vizinhanca5();
            break;
        case 5:
            vizinhanca6();
            break;
        }

        // Verifica se a solução encontrada é melhor que a atual
        if (valor_solucao_linha > valor_solucao)
        {
            // Atualiza a solução atual com a nova solução encontrada
            h = h_linha;
            z = z_linha;
            e = e_linha;
            valor_solucao = valor_solucao_linha;

            // Reinicia a exploração das vizinhanças
            i = 0;
        }
        else
        {
            // Passa para a próxima vizinhança
            i++;
        }
    }
}


/**
 * @brief Função de Perturbação baseada na Vizinhança 5
 *
 * Aplica uma perturbação à solução atual, chamando a função `perturbacao_v5()`
 * um número de vezes igual ao parâmetro `nivel`. A perturbação troca aleatoriamente
 * hubs ativos por inativos, modificando a solução atual para explorar novas regiões
 * do espaço de soluções.
 *
 * @param nivel Número de vezes que a perturbação será aplicada
 */
void Heuristicas::Perturbacao(int nivel)
{
    // Verifica se o nível de perturbação é válido
    if (nivel <= 0)
        return;

    // Aplica a perturbação 'nivel' vezes
    for (int i = 0; i < nivel; i++)
    {
        perturbacao_v5(); // Realiza a perturbação baseada na vizinhança 5
    }
}


/**
 * @brief Implementação do Smart Iterated Local Search (SmartILS)
 *
 * Esta função implementa o algoritmo SmartILS, que é uma variação do
 * Iterated Local Search (ILS). O algoritmo realiza uma busca local usando
 * o VND (Variable Neighborhood Descent) e, em seguida, aplica perturbações
 * controladas na solução para escapar de ótimos locais. A cada iteração,
 * o algoritmo tenta melhorar a melhor solução encontrada.
 *
 * @param iter_max Número máximo de iterações sem melhoria antes de encerrar o algoritmo
 */
void Heuristicas::SmartILS(int iter_max)
{
    // ========= Inicialização ========= //
    double intermediario_CPU;   // Variável para medir o tempo de CPU intermediário

    // Realiza a busca local inicial usando o VND
    VND();

    // Salva a melhor solução encontrada até o momento
    h_star = h;
    z_star = z;
    e_star = e;
    valor_solucao_star = valor_solucao;

    // Cria cópias da solução atual para restauração futura, se necessário
    vector<double> h_copia2 = h;
    vector<double> z_copia2 = z;
    vector<double> e_copia2 = e;
    double valor_solucao_copia2 = valor_solucao;

    // Inicializa os parâmetros de controle
    int iter = 0;          // Contador de iterações sem melhoria
    int nivel = 1;         // Nível de perturbação
    int nvezes = 1;        // Número de vezes que um nível de perturbação foi aplicado
    int vezesmax = 3;      // Número máximo de vezes que um nível de perturbação pode ser aplicado antes de aumentar o nível

    // Abre o arquivo para registrar os resultados intermediários
    ofstream arq_saida3("Resultados-limitantes-tempo-novo.txt", std::ofstream::app);
    if (!arq_saida3)
    {
        cerr << "Erro ao abrir o arquivo de resultados\n";
        exit(EXIT_FAILURE);
    }

    // Registra a solução inicial no arquivo
    arq_saida3 << fixed << setprecision(2)
               << valor_solucao_star << "\t"
               << fixed << setprecision(2)
               << (clock() - inicio_CPU) / CLOCKS_PER_SEC << endl;

    // ========== Loop Principal do SmartILS ========== //
    while (iter < iter_max)
    {
        // Aplica uma perturbação na solução atual com base no nível atual
        Perturbacao(nivel);

        // Realiza a busca local a partir da solução perturbada
        VND();

        // Captura o tempo de CPU após a busca local
        intermediario_CPU = clock();

        // Registra o valor da melhor solução e o tempo atual no arquivo
        arq_saida3 << fixed << setprecision(2)
                   << valor_solucao_star << "\t"
                   << fixed << setprecision(2)
                   << (intermediario_CPU - inicio_CPU) / CLOCKS_PER_SEC << endl;

        // Verifica se a nova solução encontrada é melhor que a melhor solução global
        if (valor_solucao > valor_solucao_star)
        {
            // Atualiza a melhor solução global
            h_star = h;
            z_star = z;
            e_star = e;
            valor_solucao_star = valor_solucao;

            // Reinicia os parâmetros de controle
            iter = 0;
            nivel = 1;
            nvezes = 1;

            // Atualiza as cópias da solução atual
            h_copia2 = h;
            z_copia2 = z;
            e_copia2 = e;
            valor_solucao_copia2 = valor_solucao;
        }
        else
        {
            // Incrementa o contador de iterações sem melhoria
            iter++;

            // Verifica se o nível de perturbação deve ser aumentado
            if (nvezes >= vezesmax)
            {
                nivel++;      // Aumenta o nível de perturbação
                nvezes = 1;   // Reseta o contador de vezes no nível atual
            }
            else
            {
                nvezes++;     // Incrementa o contador de vezes no nível atual
            }

            // Restaura a solução para a melhor solução conhecida antes da perturbação
            h = h_copia2;
            z = z_copia2;
            e = e_copia2;
            valor_solucao = valor_solucao_copia2;
        }
    }

    // Fecha o arquivo de resultados
    arq_saida3.close();
}

