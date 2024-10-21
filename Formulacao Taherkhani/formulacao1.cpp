/*********************************************
 * Concert Model
 * Autores: Fabricio Alves Oliveira e Henrique Heiderscheidt

 * Problem - Localização de hub com max do lucro - FORMULAÇÃO 1
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

		int n; 											// Quantidade de nós
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
		double soma = 0.0;								// Soma para obter a média geral dos custos fixos de instalação hub - Dados AP

		// ========================================================================================
		// COLETAR DADOS CAB (Exemplo de sintaxe para rodar: ./form1 cab.txt 0.8 1000 150 15)
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
		{
			z[k] = IloNumVarArray(env, n, 0, 1, ILOBOOL);
		}

		IloArray<IloArray<IloArray<IloNumVarArray>>> y(env, n); // ILOBOOL y[i][j][k][l] indica se a demanda entre os nós i e j é atendida através do caminho com primerio hub k e último hub l
		for (int i = 0; i < n; i++)
		{
			y[i] = IloArray<IloArray<IloNumVarArray>>(env, n);
			for (int j = 0; j < n; j++)
			{
				y[i][j] = IloArray<IloNumVarArray>(env, n);
				for (int k = 0; k < n; k++)
				{
					y[i][j][k] = IloNumVarArray(env, n, 0, 1, ILOBOOL);
				}
			}
		}

		IloArray<IloArray<IloNumVarArray>> f(env, n); // f[i][k][l] representa a quantidade de demanda originada no nó i e roteada no link inter hub k-l
		for (int i = 0; i < n; i++)
		{
			f[i] = IloArray<IloNumVarArray>(env, n);
			for (int k = 0; k < n; k++)
			{
				f[i][k] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);
			}
		}

		IloArray<IloNumVarArray> e(env, n); // e[i][j] indica a fração da demanda que é roteada entre os nós não concentradores - conexão direta
		for (int i = 0; i < n; i++)
			e[i] = IloNumVarArray(env, n, 0, 1, ILOBOOL);

		// ====================================Formulação 1=================================================
		// maximize sum(i in 1..n, j in 1..n, k in 1..n, l in 1..n) r[i][j] * w[i][j] * y[i][j][k][l] + r[i][j]*w[i][j]*e[i][j]
		// - [ sum(i in 1..n, j in 1..n, k in 1..n, l in 1..n) (c[i][k] + c[l][j])*  w[i][j] * y[i][j][k][l]
		// + sum(i in 1..n, k in 1..n, l in 1..n) alpha * c[i][j] * f[i][k][l] + sum(i in 1..n, j in 1..n) c[i][j] * w[i][j] * e[i][j] +  sum(k in 1..n) s[k] * h[k]
		// + sum(k in 1..n, l in 1..n) g[k][l] * z[k][l] + sum(i in 1..n, j in 1..n) q[i][j] * e[i][j]  ];
		// =================================================================================================

		IloExpr expfo(env);
		for (int k = 0; k < n; k++)
		{
			expfo -= s[k] * h[k];
			for (int l = 0; l < n; l++)
			{
				expfo -= g[k][l] * z[k][l];
				for (int i = 0; i < n; i++)
				{
					expfo -= alpha * c[k][l] * f[i][k][l];
					for (int j = 0; j < n; j++)
					{
						expfo += (r[i][j] - c[i][k] - c[l][j]) * w[i][j] * y[i][j][k][l];
					}
				}
			}
		}

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				expfo += (r[i][j] - c[i][j]) * w[i][j] * e[i][j] - q[i][j] * e[i][j];
			}
		}
		IloAdd(mod, IloMaximize(env, expfo));
		expfo.end();

		//===========================================================================
		// forall(i in 1..n, j in 1..n)  sum(k in 1..n, l in 1..n) y[i][j][k][l] + e[i][j] <= 1
		//===========================================================================
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				IloExpr r1(env);
				r1 += e[i][j];
				for (int k = 0; k < n; k++)
				{
					for (int l = 0; l < n; l++)
					{
						r1 += y[i][j][k][l];
					}
				}
				mod.add(r1 <= 1);
				r1.end();
			}
		}

		//===================================================================================================================
		// forall(i in 1..n, j in 1..n, k in 1..n)  sum(l in 1..n) y[i][j][k][l] + sum(l in 1..n, l!k) y[i][j][l][k]  <= h[k]
		//===================================================================================================================
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < n; k++)
				{
					IloExpr r2(env);
					for (int l = 0; l < n; l++)
					{
						r2 += y[i][j][k][l];
						if (l != k)
						{
							r2 += y[i][j][l][k];
						}
					}
					r2 -= h[k];
					mod.add(r2 <= 0);
					r2.end();
				}
			}
		}

		//========================================================================================================================================================================================================
		// forall(i in 1..n, k in 1..n)  sum(l in 1..n, l!k) f[i][l][k] + sum(j in 1..n, l in 1..n) w[i][j] * y[i][j][k][l]  == sum(l in 1..n, l!k) f[i][k][l] + sum(j in 1..n, l in 1..n) w[i][j] * y[i][j][l][k]
		//========================================================================================================================================================================================================
		for (int i = 0; i < n; i++)
		{
			for (int k = 0; k < n; k++)
			{
				IloExpr r3(env);
				for (int l = 0; l < n; l++)
				{
					if (l != k)
					{
						r3 += f[i][l][k] - f[i][k][l];
					}
					for (int j = 0; j < n; j++)
					{
						r3 += w[i][j] * y[i][j][k][l] - w[i][j] * y[i][j][l][k];
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
		cplex.setParam(IloCplex::Param::ClockType, 1);
		cplex.setParam(IloCplex::TiLim, 18000); // Tempo limite de resolução
		cplex.setWarning(env.getNullStream());	// Eliminar warnings
												//			cplex.setOut(env.getNullStream()); // Eliminar os logs do solver
												//  cplex.setParam(IloCplex::Threads, 4); // Definir a quantidade de threads
												// cplex.setParam(IloCplex::Param::Benders::Strategy, 3); // Ativar Benders do solver

		///==============================
		/// Resolvendo o problema
		///==============================

		IloTimer crono(env); // Variável para coletar o tempo
		double lb = 0;		 /// Callback_mipinfo
		double ub = 10e-10;	 /// Callback_mipinfo
		double gap;			 /// Callback_mipinfo
							 // double gap = 1; /// Callback_mipinfo
		// cplex.use(Callback_mipinfo(env, lb, ub, gap));/// Callback_mipinfo

		for (int i = 0; i < n; i++)
		{ // Pré-fixando variáveis
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < n; k++)
				{
					for (int l = 0; l < n; l++)
					{
						if (r[i][j] < (c[i][k] + c[l][j]))
							y[i][j][k][l].setBounds(0, 0);
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

		cout << argv[0] << " " << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << " " << argv[6] << "\n";
		FILE *re;
		re = fopen("ResultadosTaherkhani.txt", "aw+");
		fprintf(re, "Receita: %s\n", argv[3]);
		fprintf(re, "Custos: %s %s %s\n", argv[4], argv[5], argv[6]);
		fprintf(re, "Alpha: %s\n", argv[2]);
		fprintf(re, "Valor função objetivo: %f\n", (double)cplex.getObjValue());

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
					for (int l = 0; l < n; l++)
					{
						if (cplex.getValue(y[i][j][k][l]) > 0.001)
						{
							cont1 = cont1 + 1;
						}
					}
				}
			}
		};
		fprintf(re, "Demanda total atendida: %f\n", (cont1 + cont2) / 600);
		fprintf(re, "Demanda por link direto: %f\n", cont2 / 600);
		fprintf(re, "Quantidade de pares atendidos: %f\n", cont1 + cont2);
		fprintf(re, "Hubs: ");
		for (int j = 0; j < n; j++)
		{
			if (cplex.getValue(h[j]) >= 0.1)
			{
				fprintf(re, "%d ", j + 1);
			}
		};
		fprintf(re, "\n");
		fprintf(re, "Tempo de CPU: %f\n", (double)crono.getTime());
		fprintf(re, "======================================================================\n");
	}
	catch (IloException &ex)
	{
		cerr << "Error: " << ex << endl;
	}
	return 0;
}
