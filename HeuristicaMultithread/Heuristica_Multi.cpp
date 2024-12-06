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
#include <omp.h>
#include <cfloat>

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
    void ILS(int iter_max);
};

int main(int argc, char *argv[])
{
    // Cria uma instância da classe Heuristicas
    Heuristicas *ht = new Heuristicas();
    omp_set_num_threads(8);

    //========= Entrada dos dados =========//

    // Lê argumentos da linha de comando
    ht->Read_main_arg(argc, argv);
    // Lê os dados do arquivo especificado
    ht->Read_data(argv[1]);

    //========= Solução Inicial =========//

    // Inicia a contagem do tempo
    ht->inicio_CPU = clock();
    double start_time = omp_get_wtime();

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
    // Aplica o método ILS com número máximo de iterações igual a 2
    ht->ILS(2);

    // Finaliza a contagem do tempo
    double end_time = omp_get_wtime();
    ht->fim_CPU = clock();
    // Descomente a linha abaixo para imprimir o tempo de execução
    // printf("Tempo de execução = %10.2f segundos\n", (double)(ht->fim_CPU - ht->inicio_CPU) / CLOCKS_PER_SEC);
    // Descomente a linha abaixo para imprimir o valor da solução após o ILS
    // cout << "solução ILS: " << fixed << setprecision(2) << ht->valor_solucao_star << endl;

    // Abre o arquivo para escrever a configuração da rede (formato JSON)
    ofstream arq_saida2("Resultados-configuracoes.json", std::ofstream::app);
    if (!arq_saida2)
    {
        cerr << "Erro ao abrir o arquivo de configuração.\n";
        exit(0);
    }

    // Escreve os dados em formato JSON
    arq_saida2 << "{\n"
            << "  \"nome_arquivo\": \"" << argv[1] << "\",\n"
            << "  \"alpha\": " << ht->alpha << ",\n"
            << "  \"receita\": " << ht->r[0][0] << ",\n"
            << "  \"valor_inicial\": " << fixed << setprecision(2) << valor_inicial << ",\n"
            << "  \"valor_solucao\": " << fixed << setprecision(2) << ht->valor_solucao_star << ",\n"
            << "  \"tempo_processamento\": " << fixed << setprecision(2) << (double)(end_time - start_time) << ",\n";

    arq_saida2 << "  \"hubs_instalados\": [";
    bool primeiro = true;
    for (int k = 0; k < ht->n; k++)
    {
        if (ht->h_star[k] > (1 - EPSILON))
        {
            if (!primeiro) arq_saida2 << ", ";
            arq_saida2 << k;
            primeiro = false;
        }
    }
    arq_saida2 << "],\n";

    arq_saida2 << "  \"arcos_instalados\": [";
    primeiro = true;
    for (int k = 0; k < ht->n; k++)
    {
        for (int m = 0; m < ht->n; m++)
        {
            if (m != k && ht->z_star[km(ht->n, k, m)] > (1 - EPSILON))
            {
                if (!primeiro) arq_saida2 << ", ";
                arq_saida2 << "[" << k << ", " << m << "]";
                primeiro = false;
            }
        }
    }
    arq_saida2 << "]\n";

    arq_saida2 << "},\n";


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
    #pragma omp parallel for collapse(2)
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
    #pragma omp parallel for collapse(2)
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
    // Matrizes para armazenar os caminhos mais curtos
    vector<vector<double>> hub_net_sp_matrix(n, vector<double>(n, M)); // Caminhos mais curtos na rede de hubs
    vector<vector<double>> sp_matrix(n, vector<double>(n, M));         // Caminhos mais curtos para todos os nós

    double sol_value = 0.0;    // Valor da solução (lucro total)
    vector<int> hub_list;      // Lista dos hubs ativos
    e_linha = e_zero;          // Inicializa o vetor de conexões diretas

    // 1. Monta a lista de hubs ativos e subtrai os custos de instalação
    for (int k = 0; k < n; k++)
    {
        if (h[k] >= (1 - EPSILON))
        {
            hub_list.push_back(k);    // Adiciona o hub ativo à lista
            sol_value -= s[k];        // Subtrai o custo fixo de instalação do hub k
        }
    }

    // 2. Calcula os caminhos mais curtos entre os hubs na rede de hubs
    double sol_value_reduction = 0.0;

    #pragma omp parallel for collapse(2) reduction(+:sol_value_reduction)
    for (size_t m_idx = 0; m_idx < hub_list.size(); ++m_idx)
    {
        for (size_t l_idx = 0; l_idx < hub_list.size(); ++l_idx)
        {
            int m = hub_list[m_idx];
            int l = hub_list[l_idx];

            if (m != l)
            {
                if (z[km(n, m, l)] >= (1 - EPSILON))
                {
                    // Custo direto entre hubs conectados
                    hub_net_sp_matrix[m][l] = alpha * c[m][l];
                    sol_value_reduction += g[m][l];    // Acumula o custo de operação
                }
                else
                {
                    // Não há conexão direta; mantém o custo alto (M)
                    hub_net_sp_matrix[m][l] = M;
                }
            }
        }
    }

    sol_value -= sol_value_reduction; // Atualiza sol_value fora da região paralela

    // 3. Aplica o algoritmo de Floyd-Warshall
    for (size_t k_idx = 0; k_idx < hub_list.size(); ++k_idx)
    {
        int k = hub_list[k_idx];

        #pragma omp parallel for
        for (size_t m_idx = 0; m_idx < hub_list.size(); ++m_idx)
        {
            int m = hub_list[m_idx];
            for (size_t l_idx = 0; l_idx < hub_list.size(); ++l_idx)
            {
                int l = hub_list[l_idx];

                if (hub_net_sp_matrix[m][k] + hub_net_sp_matrix[k][l] < hub_net_sp_matrix[m][l])
                {
                    hub_net_sp_matrix[m][l] = hub_net_sp_matrix[m][k] + hub_net_sp_matrix[k][l];
                }
            }
        }
    }

    // 4. Calcula os caminhos mais curtos para todos os pares de nós (i, j) e o lucro
    vector<double> e_linha_local(n * n, 0); // Declarada fora da região paralela

    #pragma omp parallel reduction(+:sol_value)
    {
        // Cada thread terá sua própria cópia de e_linha_local
        vector<double> e_linha_local_private(n * n, 0);

        #pragma omp for collapse(2)
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                double sp_ij = M;  // Inicializa com um valor grande

                // Considera todos os hubs ativos como possíveis intermediários
                for (auto k : hub_list)
                {
                    // Custo passando por um único hub (k)
                    double cost_single_hub = c[i][k] + c[k][j];
                    if (cost_single_hub < sp_ij)
                    {
                        sp_ij = cost_single_hub;
                    }

                    for (auto m : hub_list)
                    {
                        if (k != m)
                        {
                            // Custo passando por dois hubs (k e m)
                            double cost_aux = c[i][k] + hub_net_sp_matrix[k][m] + c[m][j];
                            if (sp_ij > cost_aux)
                            {
                                sp_ij = cost_aux;  // Atualiza o caminho mais curto
                            }
                        }
                    }
                }

                // Calcula o lucro para cada par de nós (i, j)
                double l_ij = (r[i][j] - sp_ij) * w[i][j];
                double le_ij = (r[i][j] - c[i][j]) * w[i][j] - q[i][j];

                if (l_ij > EPSILON && l_ij > le_ij)
                {
                    sol_value += l_ij;  // Adiciona o lucro via hubs
                }
                else if (le_ij > EPSILON && le_ij > l_ij)
                {
                    sol_value += le_ij;                        // Adiciona o lucro via conexão direta
                    e_linha_local_private[km(n, i, j)] = 1;    // Marca a conexão direta no vetor local
                }
            }
        }

        // Combina e_linha_local_private em e_linha_local de forma segura
        #pragma omp critical
        {
            for (size_t idx = 0; idx < e_linha.size(); ++idx)
            {
                if (e_linha_local_private[idx] > 0)
                {
                    e_linha_local[idx] = 1;
                }
            }
        }
    }

    // Atualiza e_linha com os valores combinados
    e_linha = e_linha_local;

    return sol_value;  // Retorna o valor total da solução
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

    // Se não houver hubs ativos, ativa um hub aleatoriamente
    if (hubs_ativos.empty())
    {
        posicao = rand() % hubs_inativos.size();
        h[hubs_inativos[posicao]] = 1;
        valor_solucao = Shortest_path_algorithm();
    }
    if (!hubs_ativos.empty())
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
 * @brief Implementação do Smart Iterated Local Search (ILS)
 *
 * Esta função implementa o algoritmo ILS, que é uma variação do
 * Iterated Local Search (ILS). O algoritmo realiza uma busca local usando
 * o VND (Variable Neighborhood Descent) e, em seguida, aplica perturbações
 * controladas na solução para escapar de ótimos locais. A cada iteração,
 * o algoritmo tenta melhorar a melhor solução encontrada.
 *
 * @param iter_max Número máximo de iterações sem melhoria antes de encerrar o algoritmo
 */
void Heuristicas::ILS(int iter_max)
{
    int num_threads = omp_get_max_threads();
    h_star = h;
    z_star = z;
    e_star = e;
    valor_solucao_star = valor_solucao;

    // Cópias da melhor solução atual
    vector<double> h_copia = h_star;
    vector<double> z_copia = z_star;
    vector<double> e_copia = e_star;
    double valor_solucao_copia = valor_solucao_star;

    int iter = 0; // Contador de iterações sem melhoria

    while (iter < iter_max)
    {
        //melhores soluções encontradas por cada thread
        vector<vector<double>> h_star_vet(num_threads);
        vector<vector<double>> z_star_vet(num_threads);
        vector<vector<double>> e_star_vet(num_threads);
        vector<double> valor_solucao_star_vet(num_threads, -DBL_MAX);

        // Cada thread aplica uma perturbação diferente thread_id % 4 antes do VND
        #pragma omp parallel
        {
            int thread_id = omp_get_thread_num();

            // Cria uma cópia da melhor solução global atual para esta thread
            Heuristicas ht_thread = *this;
            ht_thread.h = h_star;
            ht_thread.z = z_star;
            ht_thread.e = e_star;
            ht_thread.valor_solucao = valor_solucao_star;

            // Aplica uma perturbação específica com base em (thread_id % 7)
            int perturbacao = (thread_id % 7) + 1;
            cout << perturbacao << endl;
            
            //1 a 7.
            ht_thread.Perturbacao(perturbacao);

            // Executa o VND após a perturbação
            ht_thread.VND();
            cout << ht_thread.valor_solucao << endl;

            // Salva a melhor solução encontrada pela thread
            h_star_vet[thread_id] = ht_thread.h;
            z_star_vet[thread_id] = ht_thread.z;
            e_star_vet[thread_id] = ht_thread.e;
            valor_solucao_star_vet[thread_id] = ht_thread.valor_solucao;
        }

        // Seleciona a melhor solução entre as threads
        double best_val_sol = valor_solucao_star;
        vector<double> best_h = h_star;
        vector<double> best_z = z_star;
        vector<double> best_e = e_star;

        for (int i = 0; i < num_threads; i++)
        {
            if (valor_solucao_star_vet[i] > best_val_sol)
            {
                best_val_sol = valor_solucao_star_vet[i];
                best_h = h_star_vet[i];
                best_z = z_star_vet[i];
                best_e = e_star_vet[i];
            }
        }

        // Verifica se houve melhoria
        if (best_val_sol > valor_solucao_star)
        {
            // Atualiza a melhor solução global
            h_star = best_h;
            z_star = best_z;
            e_star = best_e;
            valor_solucao_star = best_val_sol;

            // Reinicia o contador de iterações sem melhoria
            iter = 0;

            // Atualiza as cópias da melhor solução
            h_copia = h_star;
            z_copia = z_star;
            e_copia = e_star;
            valor_solucao_copia = valor_solucao_star;
        }
        else
        {
            cout << iter << endl;
            // Não houve melhoria
            iter++;

            // Restaura a melhor solução conhecida antes da tentativa atual
            h = h_copia;
            z = z_copia;
            e = e_copia;
            valor_solucao = valor_solucao_copia;
        }
    }
}




