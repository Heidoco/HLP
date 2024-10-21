/*********************************************
 * Concert Model
 * Autor: Henrique, Fabricio
 * Problema - Localização de hub com max do lucro com conexao direta - FORMULAÇÃO 3
 *********************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <ilcplex/ilocplex.h>

using namespace std;
ILOSTLBEGIN

ILOMIPINFOCALLBACK3(Callback_mipinfo, double &, lb, double &, ub, double &, gap)
{							/// Callback nomeada Callback_mipinfo com 3  argumentos do tipo double
	ub = getBestObjValue(); // Melhor lower bound da árvore de B&C
	if (lb < getIncumbentObjValue())
	{
		lb = getIncumbentObjValue(); // Valor da melhor solução incumbent (solução incumbent corrente)  encontrada durante o B&C
	}
	gap = getMIPRelativeGap();
}

int main(int argc, char *argv[])
{
	try
	{

		/**============================
		 *  Leitura dos dados
		 *=============================== */
		ifstream arq(argv[1]);
		if (!arq.is_open())
		{
			cout << "Error openning file: " << argv[1] << endl;
			arq.close();
			exit(EXIT_FAILURE);
		}

		int n; // Quantidade de nós
		arq >> n;
		float alpha;									// Fator de desconto no custo de transporte em um link entre hubs
		vector<double> codx(n);							// Coordenada x dos nós (AP)
		vector<double> cody(n);							// Coordenada y dos nós (AP)
		vector<vector<double>> w(n, vector<double>(n)); // declara um vetor de vetores para representar uma matriz com n linhas e n colunas - Quantidade de demanda a ser enviada entre os nós i e j
		vector<vector<double>> r(n, vector<double>(n)); // Receita obtida por enviar uma unidade de demanda entre os nós i e j
		vector<vector<double>> c(n, vector<double>(n)); // Custos por enviar uma unidade de demanda entre os nós i e j
		vector<vector<double>> q(n, vector<double>(n)); // Custos operação de uma conexão direta entre os nós i e j
		vector<double> s(n);							// Custos fixos de instalação de um hub
		vector<vector<double>> g(n, vector<double>(n)); // Custos de operação nos links inter-hubs
		float soma = 0;									// Soma para obter a média geral dos custos fixos de instalação hub - Dados AP

		// ========================================================================================
		// COLETAR DADOS CAB (Exemplo de sintaxe para rodar: ./form3 cabn.txt 0.8 1000 150 15)
		// ========================================================================================

		for (int i = 0; i < n; i++)
		{ // Coletar demanda e custos CAB
			for (int j = 0; j < n; j++)
			{
				arq >> w[i][j];
				arq >> c[i][j];
				soma = soma + w[i][j];
			}
		}

		for (int i = 0; i < n; i++)
		{ // Escalonando demanda CAB
			for (int j = 0; j < n; j++)
			{
				w[i][j] = w[i][j] / soma;
			}
		}

		if (argc >= 3) // Coletar fator de desconto nos links inter-hubs CAB
			alpha = atof(argv[2]);
		else
			alpha = 0.2;

		for (int i = 0; i < n; i++)
		{ // Coletar receita CAB
			for (int j = 0; j < n; j++)
			{
				r[i][j] = atoi(argv[3]);
			}
		}

		for (int i = 0; i < n; i++)
		{ // Coletar custo fixo de instalação CAB
			s[i] = atoi(argv[4]);
		}

		for (int i = 0; i < n; i++)
		{ // Coletar custo de operar links inter-hubs CAB
			for (int j = 0; j < n; j++)
			{
				g[i][j] = atoi(argv[5]);
			}
		}

		for (int i = 0; i < n; i++)
		{ // Coletar custo de operar links diretos CAB
			for (int j = 0; j < n; j++)
			{
				q[i][j] = atoi(argv[6]);
			}
		}

		/**============================
		 *  Declaração do modelo
		 *=============================== */

		IloEnv env;
		IloModel mod(env);
		IloCplex cplex(mod);

		IloNumVarArray h(env, n, 0, 1, ILOBOOL); // ILOBOOL h[k] indica se um hub é localizado no nó k

		IloArray<IloNumVarArray> z(env, n); // ILOBOOL z[k][l] indica se um link inter-hub é operado entre os hubs l e k
		for (int k = 0; k < n; k++)
			z[k] = IloNumVarArray(env, n, 0, 1, ILOBOOL);

		IloArray<IloArray<IloNumVarArray>> f(env, n); // f[i][k][l] representa a quantidade de demanda originada no nó i e roteada no link inter hub k-l
		for (int i = 0; i < n; i++)
		{
			f[i] = IloArray<IloNumVarArray>(env, n);
			for (int k = 0; k < n; k++)
			{
				f[i][k] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);
			}
		}

		IloArray<IloNumVarArray> a(env, n); // a[i][k] é a fração da demanda que sai do nó i e acessa a rede pelo hub k
		for (int i = 0; i < n; i++)
			a[i] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);

		IloArray<IloArray<IloNumVarArray>> b(env, n); // b[i][j][l] é a fração da demanda entre os nós i e j que é roteada através de um caminho no qual o último hub é l
		for (int i = 0; i < n; i++)
		{
			b[i] = IloArray<IloNumVarArray>(env, n);
			for (int j = 0; j < n; j++)
			{
				b[i][j] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);
			}
		}

		IloArray<IloNumVarArray> e(env, n); // e[i][j] indica a fração da demanda que é roteada entre os nós não concentradores - conexão direta
		for (int i = 0; i < n; i++)
			e[i] = IloNumVarArray(env, n, 0, 1, ILOBOOL);

		// ====================================Formulação 3===============================================================
		// maximize sum(i in 1..n, j in 1..n, l in 1..n) r[i][j] * w[i][j] * b[i][j][l]
		// - [ sum(i in 1..n, j in 1..n, k in 1..n) c[i][k] * w[i][j] * a[i][k]
		// + sum(i in 1..n, j in 1..n, l in 1..n) c[l][j] * w[i][j] * b[i][j][l]
		// + sum(i in 1..n, k in 1..n, l in 1..n) alpha * c[k][l] * f[i][k][l] + sum(k in 1..n) s[k] * h[k]
		// + sum(k in 1..n, l in 1..n) g[k][l] * z[k][l] ];

		// ===============================================================================================================
		IloExpr expfo(env);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				expfo += r[i][j] * w[i][j] * e[i][j] - c[i][j] * w[i][j] * e[i][j];
				for (int l = 0; l < n; l++)
				{
					expfo += (r[i][j] - c[l][j]) * w[i][j] * b[i][j][l];
				}
				for (int k = 0; k < n; k++)
				{
					expfo -= c[i][k] * w[i][j] * a[i][k];
				}
			}
		}

		for (int k = 0; k < n; k++)
		{
			expfo -= s[k] * h[k];
			for (int l = 0; l < n; l++)
			{
				expfo -= g[k][l] * z[k][l];
				for (int i = 0; i < n; i++)
				{
					expfo -= alpha * c[k][l] * f[i][k][l];
				}
			}
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				expfo -= q[i][j] * e[i][j];
			}
		}

		IloAdd(mod, IloMaximize(env, expfo));
		expfo.end();

		//===========================================================================
		// forall(i in 1..n, j in 1..n)  sum(k in 1..n) a[i][k] + sum(j in 1..n) e[i][j] <= 1
		//===========================================================================
		for (int i = 0; i < n; i++)
		{
			IloExpr r1(env);
			for (int k = 0; k < n; k++)
			{
				r1 += a[i][k];
			}
			mod.add(r1 <= 1);
			r1.end();
		}

		//===========================================================================
		// forall(i in 1..n, j in 1..n)  sum(l in 1..n) b[i][j][l] + e[i][j] <= 1
		//===========================================================================
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				IloExpr r2(env);
				r2 += e[i][j];
				for (int l = 0; l < n; l++)
				{
					r2 += b[i][j][l];
				}
				mod.add(r2 <= 1);
				r2.end();
			}
		}

		//======================================================================================================================================================================================
		// forall(i in 1..n, k in 1..n)  sum(j in 1..n) w[i][j] * a[i][k] + sum(l in 1..n, l!k) f[i][l][k] == sum(j in 1..n) w[i][j] * b[i][j][k] + sum(l in 1..n, l!k) f[i][k][l]
		//======================================================================================================================================================================================
		for (int i = 0; i < n; i++)
		{
			for (int k = 0; k < n; k++)
			{
				IloExpr r3(env);
				for (int j = 0; j < n; j++)
				{
					r3 += w[i][j] * (a[i][k] - b[i][j][k]);
				}
				for (int l = 0; l < n; l++)
				{
					if (l != k)
					{
						r3 += f[i][l][k] - f[i][k][l];
					}
				}
				mod.add(r3 == 0);
				r3.end();
			}
		}

		//=============================================================================================
		// forall(i in 1..n, k in 1..n, l in 1..n, l!k)  f[i][k][l] <= sum(j in 1..n) w[i][j] * z[k][l]
		//=============================================================================================
		for (int i = 0; i < n; i++)
		{
			for (int k = 0; k < n; k++)
			{
				for (int l = 0; l < n; l++)
				{
					if (l != k)
					{
						IloExpr r4(env);
						r4 += f[i][k][l];
						for (int j = 0; j < n; j++)
						{
							r4 -= w[i][j] * z[k][l];
						}
						mod.add(r4 <= 0);
						r4.end();
					}
				}
			}
		}

		//=============================================================
		// forall(i in 1..n, j in 1..n, k in 1..n) a[i][k] <= h[k]
		//=============================================================
		for (int i = 0; i < n; i++)
		{
			for (int k = 0; k < n; k++)
			{
				mod.add(a[i][k] <= h[k]);
			}
		}

		//=============================================================
		// forall(i in 1..n, j in 1..n, l in 1..n) b[i][j][l] <= h[l]
		//=============================================================
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				for (int l = 0; l < n; l++)
				{
					mod.add(b[i][j][l] <= h[l]);
				}
			}
		}

		//===================================================
		// forall(k in 1..n, l in 1..n, l!k)  z[k][l] <= h[k]
		//===================================================
		for (int k = 0; k < n; k++)
		{
			for (int l = 0; l < n; l++)
			{
				if (l != k)
					mod.add(z[k][l] <= h[k]);
			}
		}

		//===================================================
		// forall(k in 1..n, l in 1..n, l!k)  z[k][l] <= h[l]
		//===================================================
		for (int k = 0; k < n; k++)
		{
			for (int l = 0; l < n; l++)
			{
				if (l != k)
					mod.add(z[k][l] <= h[l]);
			}
		}

		//===================================================
		// forall(i in 1..n, j in 1..n, l!k)  e[i][j] + h[i] <= 1
		//===================================================
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				mod.add(e[i][j] + h[i] <= 1);
			}
		}

		//===================================================
		// forall(i in 1..n, j in 1..n, l!k)  e[i][j] + h[j] <= 1
		//===================================================
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				mod.add(e[i][j] + h[j] <= 1);
			}
		}

		/// ==========================
		/// configurações do cplex
		/// ==========================

		cplex.setParam(IloCplex::EpGap, 0.0000001); // Definindo uma tolerancia
		cplex.setParam(IloCplex::TiLim, 18000);		// Tempo limite de resolução
		cplex.setWarning(env.getNullStream());		// Eliminar warnings
		cplex.setOut(env.getNullStream());    // Eliminar os logs do solver
													//  cplex.setParam(IloCplex::Threads, 4); // Definir a quantidade de threads
													// cplex.setParam(IloCplex::Param::Benders::Strategy, 3); // Ativar Benders do solver

		///==============================
		/// Resolvendo o problema
		///==============================

		IloTimer crono(env);						   // Variável para coletar o tempo
		double lb = 0;								   /// Callback_mipinfo
		double ub = 10e10;							   /// Callback_mipinfo
		double gap = 1;								   /// Callback_mipinfo
		cplex.use(Callback_mipinfo(env, lb, ub, gap)); /// Callback_mipinfo

		for (int i = 0; i < n; i++)
		{ // Pré-fixando variáveis
			for (int j = 0; j < n; j++)
			{
				//e[i][j].setBounds(0, 0);
				for (int l = 0; l < n; l++)
				{
					if (r[i][j] < c[l][j])
					{
						b[i][j][l].setBounds(0, 0);
					}
				}
			}
		}

		crono.start();
		cplex.solve();
		crono.stop();

		if (cplex.getStatus() == IloAlgorithm::Optimal)
		{
			lb = cplex.getObjValue();
			ub = cplex.getObjValue();
			gap = 0.0;
		}

		///=====================================
		/// Salvando os resultados - CAB
		///=====================================

		printf("%s\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.6f\n", argv[1], n, alpha, 100 * gap, lb, ub, (double)crono.getTime());

		cout << " Valor função objetivo: " << cplex.getObjValue() << endl;

		printf("Hubs Instalados:");
		for (int j = 0; j < n; j++)
		{
			if (cplex.getValue(h[j]) >= 0.1)
			{
				printf(" %d\t ", j + 1);
			}
		}
		cout << endl;

		/**=====================================
		 *  Apresenta a configuração final
		 * ====================================*/

		FILE *re;
		re = fopen("ResultadosF3NCond.txt", "aw+");
		//  fprintf(re, "\n Informações Gerais: " "%s\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f \n",  argv[1], n, alpha,  r[0][0],  s[0], g[0][0] );
		fprintf(re, "\n Valor função objetivo: " "%f\t \n", (double) cplex.getObjValue ());
		fprintf(re, "\n Tempo de CPU: " "%f\t \n", (double) crono.getTime());
		//	fprintf(re, "\n lb: " "%1.2f\t \t" "ub: " "%1.2f\t \t" "gap: " "%1.2f \n", lb, ub, 100 *  gap);

		float cont1 = 0;
		float cont2 = 0;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (cplex.getValue(e[i][j]) > 0.001)
				{
					cont2 = cont2 + 1;
				}
				for (int k = 0; k < n; k++)
				{
					if (cplex.getValue(b[i][j][k]) > 0.001)
					{
						cont1 = cont1 + 1;
					}
				}
			}
		}

		fprintf(re, "%s \t %d \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t \n", argv[1], n, alpha, r[0][0], s[0], g[0][0], q[0][0], (double)cplex.getObjValue(), (double)crono.getTime());

		fprintf(re, "\n Quantidade de pares atendidos: "
					"%f \t"
					"Demanda atendida total: "
					"%f \t \n"
					"Demanda por link direto: "
					"%f \t \n",
				cont1 + cont2, (cont1 + cont2) / 600, cont2 / 600);

	// Abrindo (ou criando) um arquivo txt para salvar os resultados no modo append
ofstream result_file("resultadosNcond.json", ios::app);

// Montando o JSON diretamente no final do código
result_file << "{\n";
result_file << "  \"instance\": \"" << argv[1] << "\",\n";
result_file << "  \"n\": " << n << ",\n";
result_file << "  \"alpha\": " << alpha << ",\n";
result_file << "  \"receita\": " << argv[3] << ",\n";
result_file << "  \"custos\": [" << argv[4] <<", " << argv[5] <<", " << argv[6] << "],\n";
result_file << "  \"gap\": " << 100 * gap << ",\n";
result_file << "  \"lower_bound\": " << lb << ",\n";
result_file << "  \"upper_bound\": " << ub << ",\n";
result_file << "  \"time\": " << (double)crono.getTime() << ",\n";
result_file << "  \"objective_value\": " << cplex.getObjValue() << ",\n";

// Hubs Instalados
result_file << "  \"installed_hubs\": [";
bool first = true;
for (int j = 0; j < n; j++)
{
    if (cplex.getValue(h[j]) >= 0.1)
    {
        if (!first)
            result_file << ", ";
        result_file << (j + 1);
        first = false;
    }
}
result_file << "],\n";


result_file << "  \"pairs_served\": " << (cont1 + cont2) << ",\n";
result_file << "  \"total_demand_served\": " << (cont1 + cont2) / 600 << ",\n";
result_file << "  \"direct_link_demand\": " << cont2 / 600 << ",\n";
// Conexões feitas (indiretas e diretas)
result_file << "  \"connections\": [\n";

// Conexões diretas
result_file << "    {\n";
result_file << "      \"direct_connections\": [";
first = true;
for (int i = 0; i < n; i++)
{
    for (int j = 0; j < n; j++)
    {
        if (cplex.getValue(e[i][j]) > 0.001) // Se houver uma conexão direta
        {
            if (!first)
                result_file << ", ";
            result_file << "[" << (i + 1) << "," << (j + 1) << "]";
            first = false;
        }
    }
}
result_file << "]\n";
result_file << "    },\n";

// Conexões indiretas (via hubs)
result_file << "    {\n";
result_file << "      \"hub_connections\": [";
first = true;
for (int i = 0; i < n; i++)
{
    for (int j = 0; j < n; j++)
    {
        for (int k = 0; k < n; k++)
        {
            if (cplex.getValue(b[i][j][k]) > 0.001) // Se houver uma conexão indireta via hub
            {
                if (!first)
                    result_file << ", ";
                result_file << "{ \"from\": " << (i + 1) << ", \"to\": " << (j + 1) << ", \"via_hub\": " << (k + 1) << " }";
                first = false;
            }
        }
    }
}
result_file << "]\n";
result_file << "    }\n";

result_file << "  ]\n";

result_file << "}," << endl;

// Fechando o arquivo
result_file.close();


	}
	catch (IloException &ex)
	{
		cerr << "Error: " << ex << endl;
	}
	return 0;
}
