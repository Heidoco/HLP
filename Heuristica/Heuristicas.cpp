#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <list>
#include <math.h>
#include <time.h>

using namespace std;

#define EPSILON 1e-07 /// Parâmetro usado para representar um número pequeno
#define M 10e10       /// Parâmetro usado para representar um número suficientemente grande

/// ==============================================
/// Recuperar índices das variáveis
/// ==============================================

inline int km(int n, int k, int m)
{ /// Função usada para recuperar os índices das variáveis z_km, m!=k
    if (k < m)
    {
        return (n - 1) * k + m - 1;
    }
    else if (k > m)
    {
        return (n - 1) * k + m;
    }
}

/// ==============================================
/// Classe
/// ==============================================

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
    vector<vector<double>> q;         // Custo operação conexao direta
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
    void heuristica_construtiva(double hub_ativo, double phi);
    void heuristica_parcialmente_gulosa(double hub_ativo, double phi);
    double Shortest_path_algorithm();

    void vizinhanca1();           // insere um hub
    void vizinhanca1_com_arcos(); // insere um hub e o conecta aos demais hubs instalados
    void vizinhanca2();           // remove um hub e os arcos incidentes sobre ele
    void vizinhanca3();           // adiciona um arco hub entre hubs ativos
    void vizinhanca3_geral();     // adiciona um arcos hubs
    void vizinhanca4();           // remove um arco hub
    void vizinhanca5();           // troca um hub ativo por um não ativo e remove os arcos hubs incidentes sobre o hub desativado

    void VND();
    void Perturbacao(int nivel);
    void ILS(int iter_max);
    void vizinhanca5_perturbacao();
    void SmartILS(int iter_max);
};

/// ==============================================
/// Main
/// ==============================================

int main(int argc, char *argv[])
{

    Heuristicas *ht = new Heuristicas();

    //=========Entrada dos dados========//

    ht->Read_main_arg(argc, argv);
    ht->Read_data(argv[1]);

    //=========Solução Inicial========//

    ht->inicio_CPU = clock();

    vector<double> vetor_um_hub;
    vetor_um_hub = ht->sol_um_hub();
    cout << "hub ativo: " << vetor_um_hub[0] << endl;
    cout << "valor da solução com um hub ativo: " << vetor_um_hub[1] << endl;

    ht->h[vetor_um_hub[0]] = 1;
    ht->valor_solucao = vetor_um_hub[1];

    double valor_inicial = ht->valor_solucao;

    // Para usar uma das heurísticas, basta comentar e descomentar as seguintes partes do código:

    //=========Heurística ILS-RVND========//

    // srand(time(NULL));
    // ht->ILS(1); //testado com 4 para igualar ao smartils
    //
    // ht->fim_CPU = clock();
    // printf("Tempo execucao = %10.2f segundos\n",(double)(ht->fim_CPU - ht->inicio_CPU)/CLOCKS_PER_SEC);
    // cout << "sol ILS: " << fixed << setprecision(2) << ht->valor_solucao_star << endl;

    //=========Heurística E-ILS-RVND========//

    srand(time(NULL));
    ht->SmartILS(4);

    ht->fim_CPU = clock();
    printf("Tempo execucao = %10.2f segundos\n", (double)(ht->fim_CPU - ht->inicio_CPU) / CLOCKS_PER_SEC);
    cout << "sol SmartILS: " << fixed << setprecision(2) << ht->valor_solucao_star << endl;

    //===========Arquivo Saída===========//

    ofstream arq_saida("Resultados-heuristica-vizinhancasall.txt", std::ofstream::app);
    if (!arq_saida)
    {
        cerr << "Erro arquivo \n";
        exit(0);
    }

    arq_saida << argv[1] << "\t" << ht->alpha << "\t" << ht->r[0][0] << "\t" << fixed << setprecision(2) << valor_inicial << "\t" << fixed << setprecision(2) << ht->valor_solucao_star << "\t" << fixed << setprecision(2) << (double)(ht->fim_CPU - ht->inicio_CPU) / CLOCKS_PER_SEC << endl;

    /*
ofstream arq_saida2("Resultados-configuracoes.txt",std::ofstream::app); //Ativar se quiser imprimir a configuração da rede
    if (!arq_saida2) {cerr << "Erro arquivo \n"; exit(0);}

    arq_saida2<< "\n \n" <<argv[1]<<"\t"<<ht->alpha<<"\t"<<ht->r[0][0]<<"\t"<<fixed << setprecision(2)<<valor_inicial<<"\t"<<fixed << setprecision(2)<<ht->valor_solucao_star<<"\t"<<fixed << setprecision(2)<<(double)(ht->fim_CPU - ht->inicio_CPU)/CLOCKS_PER_SEC<<endl;

    arq_saida2<<"hubs instalados: ";
    for(int k = 0; k < ht->n; k++){
        if(ht->h_star[k] > (1-EPSILON)){
            arq_saida2 << k << "\t";
        }
    }

    arq_saida2<< "\n" << "arcos instalados:";
    for(int k = 0; k < ht->n; k++){
        for(int m = 0; m < ht->n; m++){
            if(m!=k){
                if(ht->z_star[km(ht->n, k, m)] > (1-EPSILON)){
                    arq_saida2<< "(" << k << "," << m << ")" << "\t";
                }
            }
        }
    }
    */

    return 0;
}

void Heuristicas::Read_main_arg(int argc, char *argv[])
{
    alpha1 = (argc > 2) ? atof(argv[2]) : 0.75; // valor de alpha via linha de comando
    receita = (argc > 3) ? atof(argv[3]) : 20; // valor da receita via linha de comando
}

void Heuristicas::Read_data(char name[])
{ // Leitura do conjunto de dados AP
    ifstream arq(name);
    if (!arq)
    {
        cerr << "Erro arquivo \n";
        exit(0);
    }

    arq >> n;

    codx = vector<double>(n);                         // Coordenada x dos nós (AP)
    cody = vector<double>(n);                         // Coordenada y dos nós (AP)
    w = vector<vector<double>>(n, vector<double>(n)); // declara um vetor de vetores para representar uma matriz com n linhas e n colunas - Quantidade de demanda a ser enviada entre os nós i e j
    r = vector<vector<double>>(n, vector<double>(n)); // Receita obtida por enviar uma unidade de demanda entre os nós i e j
    c = vector<vector<double>>(n, vector<double>(n)); // Custos por enviar uma unidade de demanda entre os nós i e j
    s = vector<double>(n);                            // Custos fixos de instalação de um hub
    g = vector<vector<double>>(n, vector<double>(n)); // Custos de operação nos links inter-hubs
    h = vector<double>(n, 0);                         // Inicializando vetor hubs ativos com zero
    z = vector<double>(n * (n - 1), 0);               // Inicializando vetor arcos entre hubs ativos com zero
    e = vector<double>(n * (n - 1), 0);               // Inicializando vetor de arcos entre conexoes diretas com zero
    e_registro = vector<double>(n * (n - 1), 0.0);
    e_linha = vector<double>(n * (n - 1), 0);
    e_zero = vector<double>(n * (n - 1), 0);          // Vetor e zerado
    soma = 0;                                         // Soma para obter a média geral dos custos fixos de instalação hub - Dados AP

    alpha = alpha1; // Fator de desconto entre nós inter-hubs AP
    q = vector<vector<double>>(n, vector<double>(n, 0.0)); // Initialize q


    for (int i = 0; i < n; i++)
    { // Receita AP 20,30,50
        for (int j = 0; j < n; j++)
        {
            r[i][j] = receita;
        }
    }

    for (int i = 0; i < n; i++)
    { // Coletar as coordenadas dos nós
        arq >> codx[i];
        arq >> cody[i];
    }

    for (int i = 0; i < n; i++)
    { // Coletar custos AP
        for (int j = 0; j < n; j++)
        {
            c[i][j] = 0.001 * sqrt((codx[i] - codx[j]) * (codx[i] - codx[j]) + (cody[i] - cody[j]) * (cody[i] - cody[j]));
        }
    }

    for (int i = 0; i < n; i++)
    { // Coletar demanda AP
        for (int j = 0; j < n; j++)
        {
            arq >> w[i][j];
        }
    }

    for (int i = 0; i < n; i++)
    { // Coletar custos fixos instalação hub AP
        arq >> s[i];
        s[i] = s[i] / 10;
        soma = soma + s[i];
    }

    for (int i = 0; i < n; i++)
    { // Custos de operação link entre hubs AP
        for (int j = 0; j < n; j++)
        {
            g[i][j] = 0.1 * soma / n;
            q[i][j] = g[i][j]*0.20;
        }
    }

    arq.close();
}

vector<double> Heuristicas::sol_um_hub()
{ // Solução inicial

    double phi = 0.0;
    double phi_aux_1 = 0.0;
    double phi_aux_2 = 0.0;
    double hub_ativo = -1;
    vector<vector<double>> l(n, vector<double>(n));
    vector<vector<double>> le(n, vector<double>(n));
    vector<double> vet_sol_um_hub(2);

    for (int k = 0; k < n; k++)
    {
        phi_aux_1 = -s[k];

        e_linha = e_zero; // Reiniciar o vetor e_linha

            
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i != j){
                    l[i][j] = r[i][j] * w[i][j] - ((c[i][k] + c[k][j]) * w[i][j]);
                    // Lucro por conexao direta
                    le[i][j] = r[i][j] * w[i][j] - c[i][j] * w[i][j] - q[i][j];
                    // cout << l[i][j] << endl;
                    // cout << le[i][j] << endl;
                    if (l[i][j] > EPSILON || le[i][j] > EPSILON)
                    {
                        if (l[i][j] > le[i][j])
                        {
                            phi_aux_1 += l[i][j];
                        }
                        else
                        {
                            phi_aux_2 += le[i][j];
                            e_linha[km(n, i, j)] = 1;
                        }
                    }
                }
            }
        }
        if ((phi_aux_1 + phi_aux_2) >= phi)
        {
            e = e_linha;
            phi = std::max(0.0, phi_aux_1 + phi_aux_2); //MUDEI PRA VER, NAO FOI
            if (phi_aux_1 > 0)
            {
                hub_ativo = k;
            }
        }
    }
    vet_sol_um_hub[0] = hub_ativo; // ALTERAR DEPOIS NAO RETORNAR, MANTER O VETOR
    vet_sol_um_hub[1] = phi;

    if (hub_ativo >= 0)
        h[hub_ativo] = 1;

    return vet_sol_um_hub;
}

double Heuristicas::Shortest_path_algorithm()
{ // Calculo do custo mínimo/Avaliação da função objetivo

    vector<vector<double>> hub_net_sp_matrix = vector<vector<double>>(n, vector<double>(n)); /// Shostest path considering only the hub network
    vector<vector<double>> sp_matrix = vector<vector<double>>(n, vector<double>(n));         /// Shostest path for all nodes

    ///===================================
    /// Calculate the shortest path from each pair of hubs considering only the hub network
    ///===================================
    double sol_value = 0.0; /// Variável que computa o custo de instalação de hubs e de arcos entre hubs
    list<int> hub_list;     /// Lista de hubs

    for (int k = 0; k < n; k++)
    { /// Monta uma lista de hubs
        if (h[k] >= (1 - EPSILON))
        {
            hub_list.insert(hub_list.end(), k);
            sol_value -= s[k];
        }
    }

    for (list<int>::iterator k = hub_list.begin(); k != hub_list.end(); k++)
    { /// Calculo do caminho mais curto na rede de hubs com apenas dois hubs
        for (list<int>::iterator m = hub_list.begin(); m != hub_list.end(); m++)
        {
            if (k != m)
            {
                if (z[km(n, *k, *m)] >= (1 - EPSILON))
                {
                    hub_net_sp_matrix[*k][*m] = alpha * c[*k][*m];
                    sol_value -= g[*k][*m];
                }
                else
                {
                    hub_net_sp_matrix[*k][*m] = M;
                }
            }
        }
    }
    for (list<int>::iterator k = hub_list.begin(); k != hub_list.end(); k++)
    {
        for (list<int>::iterator m = hub_list.begin(); m != hub_list.end(); m++)
        {
            for (list<int>::iterator l = hub_list.begin(); l != hub_list.end(); l++)
            {
                if (m != l and m != k and l != k)
                {
                    double cost_aux = hub_net_sp_matrix[*m][*k] + hub_net_sp_matrix[*k][*l];
                    if ((cost_aux < hub_net_sp_matrix[*m][*l]))
                    {
                        hub_net_sp_matrix[*m][*l] = cost_aux;
                    }
                }
            }
        }
    }

    /// Compute shortest path for all pair i,j, also allowing the direct connections

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                sp_matrix[i][j] = c[i][j]; /// (Esse M não permite conexões diretas nesse caso. ???)
                for (list<int>::iterator k = hub_list.begin(); k != hub_list.end(); k++)
                {
                    if (c[i][*k] + c[*k][j] > sp_matrix[i][j])
                    {
                        e_registro[km(n, i, j)] = 1;
                    }

                    sp_matrix[i][j] = min(c[i][*k] + c[*k][j], sp_matrix[i][j]);
                    for (list<int>::iterator m = hub_list.begin(); m != hub_list.end(); m++)
                    {
                        if (k != m)
                        {
                            double cost_aux = c[i][*k] + hub_net_sp_matrix[*k][*m] + c[*m][j];
                            if (sp_matrix[i][j] > cost_aux)
                            {
                                sp_matrix[i][j] = cost_aux;
                            }
                        }
                    }
                }
                if (r[i][j] - sp_matrix[i][j] > EPSILON)
                {
                    if (e_registro[km(n, i, j)] == 1)
                    {
                        sol_value -= q[i][j];
                    }
                    sol_value += (r[i][j] - sp_matrix[i][j]) * w[i][j];
                }
            }
        }
    }

    hub_net_sp_matrix.clear();
    hub_list.clear();
    return sol_value;
}

// Estruturas de vizinhança

void Heuristicas::vizinhanca1()
{

    vector<double> hubs_inativos;
    double f_avaliacao = 0.0;
    int u; // tamanho da vizinhanca que será explorada
    int posicao;

    h_copia = h;
    h_linha = h;
    valor_solucao_linha = valor_solucao;
    z_linha = z;

    cout << "sol_linha inicial: " << valor_solucao_linha << endl;

    for (int k = 0; k < n; k++)
    { // coleta os hubs inativos
        if (h[k] < EPSILON)
        {
            hubs_inativos.push_back(k);
            // cout << "posicao inativos: " << k << endl;
        }
    }

    if (!hubs_inativos.empty())
    {
        u = 1 * (hubs_inativos.size());
        //cout << "valor de m: " << m << endl;
        for (int j = 0; j < u; j++)
        {
            if (j == 0)
            {
                posicao = rand() % hubs_inativos.size();
                cout << "posicao 0: " << hubs_inativos[posicao] << endl;
                h[hubs_inativos[posicao]] = 1;
                valor_solucao_linha = Shortest_path_algorithm();
                h_linha = h;
                hubs_inativos.erase(hubs_inativos.begin() + posicao);
                h = h_copia;
            }
            else
            {
                posicao = rand() % hubs_inativos.size();
                //  cout << "posicao depois: " << hubs_inativos[posicao] << endl;
                h[hubs_inativos[posicao]] = 1;
                f_avaliacao = Shortest_path_algorithm();
                // cout << "f_avaliacao: " << f_avaliacao << endl;
                // cout << "sol_linha: " << valor_solucao_linha << endl;
                if (f_avaliacao > valor_solucao_linha)
                {
                    h_linha = h;
                    valor_solucao_linha = f_avaliacao;
                }
                hubs_inativos.erase(hubs_inativos.begin() + posicao);
                h = h_copia;
            }
        }
    }
    hubs_inativos.clear();
}

void Heuristicas::vizinhanca1_com_arcos()
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

    vector<double> hubs_ativos;
    vector<int> arc_hub1;
    vector<int> arc_hub2;
    double f_avaliacao = 0.0;
    int u; // tamanho da vizinhanca que será explorada
    int posicao;

    h_copia = h;
    h_linha = h;
    z_copia = z;
    z_linha = z;
    valor_solucao_linha = valor_solucao;

    // cout << "sol_linha inicial vizinhanca 2: " << valor_solucao_linha << endl;

    for (int k = 0; k < n; k++)
    { // coletando os hubs ativos
        if (h[k] > (1 - EPSILON))
        {
            hubs_ativos.push_back(k);
            // cout << "hubs ativos vizinhanca 2: " << k << endl;
        }
    }
    // cout << "fim" << endl;

    for (int k = 0; k < n; k++)
    { // coletando arcos entre hubs ativos
        for (int m = 0; m < n; m++)
        {
            if (m != k)
            {
                if (z[km(n, k, m)] > (1 - EPSILON))
                {
                    arc_hub1.push_back(k);
                    arc_hub2.push_back(m);
                    // cout << "arcos ativos viz2: " << k << " , " << m << endl;
                }
            }
        }
    }

    if (hubs_ativos.size() > 1)
    {

        u = 1 * hubs_ativos.size();
        for (int j = 0; j < u; j++)
        {
            if (j == 0)
            {
                posicao = rand() % hubs_ativos.size();
                h[hubs_ativos[posicao]] = 0;
                // cout << "hub retirado: " << hubs_ativos[posicao] << endl;
                h_linha = h;
                if (!arc_hub1.empty())
                {
                    for (int i = 0; i < arc_hub1.size(); i++)
                    {
                        if ((arc_hub1[i] == hubs_ativos[posicao]) || (arc_hub2[i] == hubs_ativos[posicao]))
                        {
                            z[km(n, arc_hub1[i], arc_hub2[i])] = 0;
                            // cout << "arcos retirados: " << arc_hub1[i] << " , " << arc_hub2[i] << endl;
                        }
                    }
                }
                z_linha = z;
                valor_solucao_linha = Shortest_path_algorithm();
                hubs_ativos.erase(hubs_ativos.begin() + posicao);
                h = h_copia;
                z = z_copia;
            }
            else
            {
                posicao = rand() % hubs_ativos.size();
                h[hubs_ativos[posicao]] = 0;
                // cout << "hub retirado: " << hubs_ativos[posicao] << endl;
                if (!arc_hub1.empty())
                {
                    for (int i = 0; i < arc_hub1.size(); i++)
                    {
                        if ((arc_hub1[i] == hubs_ativos[posicao]) || (arc_hub2[i] == hubs_ativos[posicao]))
                        {
                            z[km(n, arc_hub1[i], arc_hub2[i])] = 0;
                            // cout << "arcos retirados: " << arc_hub1[i] << " , " << arc_hub2[i] << endl;
                        }
                    }
                }
                f_avaliacao = Shortest_path_algorithm();
                // cout << "f_avaliacao viz2: " << f_avaliacao << endl;
                // cout << "sol_linha viz2: " << valor_solucao_linha << endl;
                if (f_avaliacao > valor_solucao_linha)
                {
                    h_linha = h;
                    z_linha = z;
                    valor_solucao_linha = f_avaliacao;
                }
                hubs_ativos.erase(hubs_ativos.begin() + posicao);
                h = h_copia;
                z = z_copia;
            }
        }
    }
    hubs_ativos.clear();
    arc_hub1.clear();
    arc_hub2.clear();
}

void Heuristicas::vizinhanca3()
{

    vector<double> hubs_ativos;
    vector<int> arc_hub1;
    vector<int> arc_hub2;
    double f_avaliacao = 0.0;
    int u; // tamanho da vizinhanca que será explorada
    int posicao;

    h_linha = h;
    z_copia = z;
    z_linha = z;
    valor_solucao_linha = valor_solucao;

    // cout << "sol_linha inicial vizinhanca 3: " << valor_solucao_linha << endl;

    for (int k = 0; k < n; k++)
    { // coletando os hubs ativos
        if (h[k] > (1 - EPSILON))
        {
            hubs_ativos.push_back(k);
            // cout << "hubs ativos vizinhanca 3: " << k << endl;
        }
    }

    if (hubs_ativos.size() > 1)
    { // conjunto de arcos entre hubs possíveis
        for (int i = 0; i < hubs_ativos.size(); i++)
        {
            for (int j = 0; j < hubs_ativos.size(); j++)
            {
                if ((i != j) && (z[km(n, hubs_ativos[i], hubs_ativos[j])] < EPSILON))
                {
                    arc_hub1.push_back(hubs_ativos[i]);
                    arc_hub2.push_back(hubs_ativos[j]);
                    // cout << "arcos disponiveis vizinhanca 3: " << hubs_ativos[i] << "," << hubs_ativos[j] << endl;
                }
            }
        }
        u = 1 * arc_hub1.size();
        for (int j = 0; j < u; j++)
        {
            if (j == 0)
            {
                posicao = rand() % arc_hub1.size();
                z[km(n, arc_hub1[posicao], arc_hub2[posicao])] = 1;
                z_linha = z;
                valor_solucao_linha = Shortest_path_algorithm();
                // cout << "arco adicionado vizinhanca 3: " << arc_hub1[posicao] << "," << arc_hub2[posicao] << endl;
                arc_hub1.erase(arc_hub1.begin() + posicao);
                arc_hub2.erase(arc_hub2.begin() + posicao);
                z = z_copia;
            }
            else
            {
                posicao = rand() % arc_hub1.size();
                z[km(n, arc_hub1[posicao], arc_hub2[posicao])] = 1;
                f_avaliacao = Shortest_path_algorithm();
                // cout << "arco adicionado vizinhanca 3: " << arc_hub1[posicao] << "," << arc_hub2[posicao] << endl;
                // cout << "f_avaliacao vizinhanca 3: " << f_avaliacao << endl;
                // cout << "valor solucao_linha vizinhanca 3: " << valor_solucao_linha << endl;
                if (f_avaliacao > valor_solucao_linha)
                {
                    z_linha = z;
                    valor_solucao_linha = f_avaliacao;
                }
                arc_hub1.erase(arc_hub1.begin() + posicao);
                arc_hub2.erase(arc_hub2.begin() + posicao);
                z = z_copia;
            }
        }
    }
    hubs_ativos.clear();
    arc_hub1.clear();
    arc_hub2.clear();
}

void Heuristicas::vizinhanca3_geral()
{

    vector<int> arc_hub1;
    vector<int> arc_hub2;
    double f_avaliacao = 0.0;
    int u; // tamanho da vizinhanca que será explorada
    int posicao;

    h_copia = h;
    h_linha = h;
    z_copia = z;
    z_linha = z;
    valor_solucao_linha = valor_solucao;

    // cout << "sol_linha inicial vizinhanca 3 geral: " << valor_solucao_linha << endl;

    for (int k = 0; k < n; k++)
    {
        for (int m = 0; m < n; m++)
        {
            if ((m != k) && (z[km(n, k, m)] < EPSILON))
            {
                arc_hub1.push_back(k);
                arc_hub2.push_back(m);
                // cout << "arcos disponiveis vizinhanca 3 geral: " << k << "," << m << endl;
            }
        }
    }

    u = 1 * arc_hub1.size();

    for (int j = 0; j < u; j++)
    {
        if (j == 0)
        {
            posicao = rand() % arc_hub1.size();
            z[km(n, arc_hub1[posicao], arc_hub2[posicao])] = 1;
            // cout << "arco ativado: " << arc_hub1[posicao] << "," << arc_hub2[posicao] << endl;
            z_linha = z;
            h[arc_hub1[posicao]] = 1;
            h[arc_hub2[posicao]] = 1;
            h_linha = h;
            valor_solucao_linha = Shortest_path_algorithm();
            arc_hub1.erase(arc_hub1.begin() + posicao);
            arc_hub2.erase(arc_hub2.begin() + posicao);
            h = h_copia;
            z = z_copia;
        }
        else
        {
            posicao = rand() % arc_hub1.size();
            z[km(n, arc_hub1[posicao], arc_hub2[posicao])] = 1;
            // cout << "arco ativado: " << arc_hub1[posicao] << "," << arc_hub2[posicao] << endl;
            h[arc_hub1[posicao]] = 1;
            h[arc_hub2[posicao]] = 1;
            f_avaliacao = Shortest_path_algorithm();
            // cout << "valor solucao_linha: " << valor_solucao_linha << endl;
            // cout << "f_avaliacao: " << f_avaliacao << endl;
            if (f_avaliacao > valor_solucao_linha)
            {
                h_linha = h;
                z_linha = z;
                valor_solucao_linha = f_avaliacao;
            }
            arc_hub1.erase(arc_hub1.begin() + posicao);
            arc_hub2.erase(arc_hub2.begin() + posicao);
            h = h_copia;
            z = z_copia;
        }
    }
    arc_hub1.clear();
    arc_hub2.clear();
}

void Heuristicas::vizinhanca4()
{

    vector<int> arc_hub1;
    vector<int> arc_hub2;
    double f_avaliacao = 0.0;
    int u; // tamanho da vizinhanca que será explorada
    int posicao;

    h_linha = h;
    z_copia = z;
    z_linha = z;
    valor_solucao_linha = valor_solucao;

    // cout << "sol_linha inicial vizinhanca 4: " << valor_solucao_linha << endl;

    for (int k = 0; k < n; k++)
    { // coletando arcos entre hubs ativos
        for (int m = 0; m < n; m++)
        {
            if (m != k)
            {
                if (z[km(n, k, m)] > (1 - EPSILON))
                {
                    arc_hub1.push_back(k);
                    arc_hub2.push_back(m);
                    // cout << "arcos ativos viz4: " << k << " , " << m << endl;
                }
            }
        }
    }

    if (!arc_hub1.empty())
    {
        u = 1 * arc_hub1.size();
        for (int j = 0; j < u; j++)
        {
            if (j == 0)
            {
                posicao = rand() % arc_hub1.size();
                z[km(n, arc_hub1[posicao], arc_hub2[posicao])] = 0;
                // cout << "arco desativado viz4: " << arc_hub1[posicao] << " , " << arc_hub2[posicao] << endl;
                z_linha = z;
                valor_solucao_linha = Shortest_path_algorithm();
                arc_hub1.erase(arc_hub1.begin() + posicao);
                arc_hub2.erase(arc_hub2.begin() + posicao);
                z = z_copia;
            }
            else
            {
                posicao = rand() % arc_hub1.size();
                z[km(n, arc_hub1[posicao], arc_hub2[posicao])] = 0;
                // cout << "arco desativado viz4: " << arc_hub1[posicao] << " , " << arc_hub2[posicao] << endl;
                f_avaliacao = Shortest_path_algorithm();
                // cout << "f_avaliacao viz4: " << f_avaliacao << endl;
                // cout << "valor_solucao_linha viz4: " << valor_solucao_linha << endl;
                if (f_avaliacao > valor_solucao_linha)
                {
                    z_linha = z;
                    valor_solucao_linha = f_avaliacao;
                }
                arc_hub1.erase(arc_hub1.begin() + posicao);
                arc_hub2.erase(arc_hub2.begin() + posicao);
                z = z_copia;
            }
        }
    }
    arc_hub1.clear();
    arc_hub2.clear();
}

void Heuristicas::vizinhanca5()
{

    vector<double> hubs_ativos;
    vector<double> hubs_inativos;
    vector<double> hubs_inativos_copia;
    vector<int> arc_hub1;
    vector<int> arc_hub2;
    double f_avaliacao = 0.0;
    int u; // tamanho da vizinhanca que será explorada
    int posicao;

    h_copia = h;
    h_linha = h;
    z_copia = z;
    z_linha = z;
    valor_solucao_linha = valor_solucao;

    // cout << "sol_linha inicial vizinhanca 5: " << valor_solucao_linha << endl;

    for (int k = 0; k < n; k++)
    { // coletando os hubs ativos
        if (h[k] > (1 - EPSILON))
        {
            hubs_ativos.push_back(k);
            // cout << "hubs ativos vizinhanca 5: " << k << endl;
        }
        else
        {
            hubs_inativos.push_back(k);
            // cout << "hubs inativos vizinhanca 5: " << k << endl;
        }
    }

    hubs_inativos_copia = hubs_inativos;

    for (int k = 0; k < n; k++)
    { // coletando arcos entre hubs ativos
        for (int m = 0; m < n; m++)
        {
            if (m != k)
            {
                if (z[km(n, k, m)] > (1 - EPSILON))
                {
                    arc_hub1.push_back(k);
                    arc_hub2.push_back(m);
                    // cout << "arcos ativos viz5: " << k << " , " << m << endl;
                }
            }
        }
    }

    if (hubs_ativos.empty())
    { // Caso não tenha nenhum hub instalado, será ativado um hub aleatoriamente e a vizinhança5 é reiniciada
        posicao = rand() % hubs_inativos.size();
        h[hubs_inativos[posicao]] = 1;
        valor_solucao = Shortest_path_algorithm();
        vizinhanca5();
        cout << "entrou aqui" << endl;
    }
    else
    {
        for (int i = 0; i < hubs_ativos.size(); i++)
        {
            u = 1 * hubs_inativos.size();
            for (int j = 0; j < u; j++)
            {
                if (i == 0 && j == 0)
                {
                    h[hubs_ativos[i]] = 0;
                    // cout << "hub desativado viz5: " << hubs_ativos[i] << endl;
                    posicao = rand() % hubs_inativos.size();
                    h[hubs_inativos[posicao]] = 1;
                    // cout << "hub ativado viz5: " << hubs_inativos[posicao] << endl;
                    h_linha = h;
                    if (!arc_hub1.empty())
                    {
                        for (int k = 0; k < arc_hub1.size(); k++)
                        {
                            if ((arc_hub1[k] == hubs_ativos[i]) || (arc_hub2[k] == hubs_ativos[i]))
                            {
                                z[km(n, arc_hub1[k], arc_hub2[k])] = 0;
                                // cout << "arcos retirados viz5: " << arc_hub1[k] << " , " << arc_hub2[k] << endl;
                            }
                        }
                    }
                    z_linha = z;
                    valor_solucao_linha = Shortest_path_algorithm();
                    hubs_inativos.erase(hubs_inativos.begin() + posicao);
                    h = h_copia;
                    z = z_copia;
                }
                else
                {
                    h[hubs_ativos[i]] = 0;
                    // cout << "hub desativado viz5: " << hubs_ativos[i] << endl;
                    posicao = rand() % hubs_inativos.size();
                    h[hubs_inativos[posicao]] = 1;
                    // cout << "hub ativado viz5: " << hubs_inativos[posicao] << endl;
                    if (!arc_hub1.empty())
                    {
                        for (int k = 0; k < arc_hub1.size(); k++)
                        {
                            if ((arc_hub1[k] == hubs_ativos[i]) || (arc_hub2[k] == hubs_ativos[i]))
                            {
                                z[km(n, arc_hub1[k], arc_hub2[k])] = 0;
                                // cout << "arcos retirados viz5: " << arc_hub1[k] << " , " << arc_hub2[k] << endl;
                            }
                        }
                    }
                    f_avaliacao = Shortest_path_algorithm();
                    // cout << "f_avaliacao viz5: " << f_avaliacao << endl;
                    // cout << "valor_solucao_linha viz5: " << valor_solucao_linha << endl;
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
            hubs_inativos = hubs_inativos_copia;
        }
    }
    hubs_ativos.clear();
    hubs_inativos.clear();
    arc_hub1.clear();
    arc_hub2.clear();
}

void Heuristicas::vizinhanca5_perturbacao()
{

    vector<double> hubs_ativos;
    vector<double> hubs_inativos;
    vector<double> hubs_inativos_copia;
    vector<int> arc_hub1;
    vector<int> arc_hub2;
    int posicao;
    int posicao2;

    for (int k = 0; k < n; k++)
    { // coletando os hubs ativos
        if (h[k] > (1 - EPSILON))
        {
            hubs_ativos.push_back(k);
        }
        else
        {
            hubs_inativos.push_back(k);
        }
    }

    for (int k = 0; k < n; k++)
    { // coletando arcos entre hubs ativos
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

    if (hubs_ativos.empty())
    { // Caso não tenha nenhum hub instalado, será ativado um hub aleatoriamente e a vizinhança5 é reiniciada
        posicao = rand() % hubs_inativos.size();
        h[hubs_inativos[posicao]] = 1;
        valor_solucao = Shortest_path_algorithm();
        vizinhanca5_perturbacao();
    }
    else
    {
        posicao = rand() % hubs_ativos.size();
        h[hubs_ativos[posicao]] = 0;
        posicao2 = rand() % hubs_inativos.size();
        h[hubs_inativos[posicao2]] = 1;
        if (!arc_hub1.empty())
        {
            for (int k = 0; k < arc_hub1.size(); k++)
            {
                if ((arc_hub1[k] == hubs_ativos[posicao]) || (arc_hub2[k] == hubs_ativos[posicao]))
                {
                    z[km(n, arc_hub1[k], arc_hub2[k])] = 0;
                }
            }
        }
        valor_solucao = Shortest_path_algorithm();
    }
    hubs_ativos.clear();
    hubs_inativos.clear();
    arc_hub1.clear();
    arc_hub2.clear();
}

void Heuristicas::VND()
{

    //================VND=========//
    /*
    int i = 1;

    while(i <= 6){
        if(i == 1){
            vizinhanca1();
        }else if(i == 2){
            vizinhanca3();
        }else if(i == 3){
            vizinhanca4();
        }else if(i == 4){
            vizinhanca1_com_arcos();
        }else if(i == 5){
            vizinhanca2();
        }else {
            vizinhanca5();
        }
        //cout << "sol vizinhanca " << i << ":" << valor_solucao_linha << endl;
        //cout << "valor solucao " << valor_solucao << endl;
        if(valor_solucao_linha > valor_solucao){
            h = h_linha;
            z = z_linha;
            valor_solucao = valor_solucao_linha;
            i = 1;
        }else{
            i = i + 1;
        }
    }
    */

    //===============RVND=========//
    // double intermediario_CPU;
    vector<int> vetor;

    for (int i = 0; i < 6; i++)
    {
        vetor.push_back(i);
        // cout << "vetor " << vetor[i] << endl;
    }

    int aux, j1, j2;

    for (int i = 0; i < 6; i++)
    {
        j1 = rand() % (6);
        j2 = rand() % (6);
        while (j1 == j2)
            j2 = rand() % (6);
        aux = vetor[j1];
        vetor[j1] = vetor[j2];
        vetor[j2] = aux;
    }
    // cout << "sol antes: " << valor_solucao << endl;
    // cout << "ordem nova:" << endl;
    // for(int i = 0; i < 6; i++){
    //     cout << vetor[i] << endl;
    // }
    // intermediario_CPU = clock();
    // printf("Tempo antes = %10.2f segundos\n",(double)(intermediario_CPU - inicio_CPU)/CLOCKS_PER_SEC);
    // cout << "-----------" << endl;

    int i = 0;

    while (i < 6)
    {
        if (vetor[i] == 0)
        {
            vizinhanca1();
        }
        else if (vetor[i] == 1)
        {
            vizinhanca3();
        }
        else if (vetor[i] == 2)
        {
            vizinhanca4();
        }
        else if (vetor[i] == 3)
        {
            vizinhanca1_com_arcos();
        }
        else if (vetor[i] == 4)
        {
            vizinhanca2();
        }
        else
        {
            vizinhanca5();
        }
        // cout << "sol vizinhanca " << i << ":" << valor_solucao_linha << endl;
        // cout << "valor solucao " << valor_solucao << endl;
        if (valor_solucao_linha > valor_solucao)
        {
            h = h_linha;
            z = z_linha;
            valor_solucao = valor_solucao_linha;
            i = 0;
        }
        else
        {
            i = i + 1;
        }
    }
    // cout << "sol depois: " << valor_solucao << endl;
    // intermediario_CPU = clock();
    // printf("Tempo depois = %10.2f segundos\n",(double)(intermediario_CPU - inicio_CPU)/CLOCKS_PER_SEC);

} // VND e RVND: retorna h, z e valor_solucao

void Heuristicas::Perturbacao(int nivel)
{

    int i = 1;

    while (i <= nivel)
    {
        vizinhanca5_perturbacao(); // retorna novo h, z e valor_solucao
        i = i + 1;
    }

} // Perturbacao: retorna h, z e valor_solucao

void Heuristicas::ILS(int iter_max)
{

    //=========VND inicial========//
    double intermediario_CPU;
    VND();
    // cout << "sol VND inicial: " << valor_solucao << endl;
    h_star = h;
    z_star = z;
    valor_solucao_star = valor_solucao;
    vector<double> h_copia2 = h;
    vector<double> z_copia2 = z;
    double valor_solucao_copia2 = valor_solucao;
    int iter = 0;
    int nivel = 1;

    ofstream arq_saida3("Resultados-limitantes-tempo-novo.txt", std::ofstream::app);
    if (!arq_saida3)
    {
        cerr << "Erro arquivo \n";
        exit(0);
    }

    arq_saida3 << fixed << setprecision(2) << valor_solucao_star << "\t" << fixed << setprecision(2) << (clock() - inicio_CPU) / CLOCKS_PER_SEC << endl;

    //==========Loop Principal ILS==========//

    while (iter < iter_max)
    {
        Perturbacao(nivel);
        // cout << "valor solucao perturbada: " << valor_solucao << endl;
        VND();

        // Coletar limitantes e tempo
        intermediario_CPU = clock();

        arq_saida3 << fixed << setprecision(2) << valor_solucao_star << "\t" << fixed << setprecision(2) << (intermediario_CPU - inicio_CPU) / CLOCKS_PER_SEC << endl;

        if (valor_solucao > valor_solucao_star)
        {
            h_star = h;
            z_star = z;
            valor_solucao_star = valor_solucao;
            iter = 0;
            nivel = 1;
            h_copia2 = h;
            z_copia2 = z;
            valor_solucao_copia2 = valor_solucao;
            // cout << "reseta o nivel para: " << nivel << endl;
        }
        else
        {
            iter = iter + 1;
            nivel = nivel + 1;
            h = h_copia2;
            z = z_copia2;
            valor_solucao = valor_solucao_copia2;
            // cout << "aumenta o nivel para: " << nivel << endl;
        }
    }
}

void Heuristicas::SmartILS(int iter_max)
{

    //=========VND inicial========//
    double intermediario_CPU;
    VND();
    cout << "sol VND inicial: " << valor_solucao << endl;
    h_star = h;
    z_star = z;
    valor_solucao_star = valor_solucao;
    vector<double> h_copia2 = h;
    vector<double> z_copia2 = z;
    double valor_solucao_copia2 = valor_solucao;
    int iter = 0;
    int nivel = 1;
    int nvezes = 1;
    int vezesmax = 3;

    ofstream arq_saida3("Resultados-limitantes-tempo-novo.txt", std::ofstream::app);
    if (!arq_saida3)
    {
        cerr << "Erro arquivo \n";
        exit(0);
    }

    arq_saida3 << fixed << setprecision(2) << valor_solucao_star << "\t" << fixed << setprecision(2) << (clock() - inicio_CPU) / CLOCKS_PER_SEC << endl;

    //==========Loop Principal SmartILS==========//

    while (iter < iter_max)
    {
        Perturbacao(nivel);
        cout << "valor solucao perturbada: " << valor_solucao << endl;

        cout << "valor solucao star: " << valor_solucao_star << endl;
        VND();

        // Coletar limitantes e tempo
        intermediario_CPU = clock();

        arq_saida3 << fixed << setprecision(2) << valor_solucao_star << "\t" << fixed << setprecision(2) << (intermediario_CPU - inicio_CPU) / CLOCKS_PER_SEC << endl;

        if (valor_solucao > valor_solucao_star)
        {
            h_star = h;
            z_star = z;
            valor_solucao_star = valor_solucao;
            iter = 0;
            nivel = 1;
            h_copia2 = h;
            z_copia2 = z;
            valor_solucao_copia2 = valor_solucao;
            // cout << "reseta o nivel para: " << nivel << endl;
        }
        else
        {
            iter = iter + 1;
            if (nvezes >= vezesmax)
            {
                nivel = nivel + 1;
                nvezes = 1;
            }
            else
            {
                nvezes = nvezes + 1;
            }
            h = h_copia2;
            z = z_copia2;
            valor_solucao = valor_solucao_copia2;
            // cout << "aumenta o nivel para: " << nivel << endl;
        }
    }
}
