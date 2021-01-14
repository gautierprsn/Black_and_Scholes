#include "methode_abstrait.h"
#include <iostream>
#include <chrono>
#include <thread>
#include <cmath>
#include <fstream>

namespace ensiie {
	//constructeur de la classe abstraite permettant d'initialiser dT et dL
	Resolv::Resolv(float T, float L,  int N, int M) {
		dT = T / M;
		dL = L / N;
	}

	//Modele Crank Nicholson 
	void Resolv::computeCoef(float* a, float* b, float* c, float* d, int M, float sigma, float r) {
		//les coefficients ont été déterminés par calcul manuel (voir LaTex)
		for (int i = 0; i < M; i++) {
			a[i] = static_cast<float> (pow(sigma, 2)) *static_cast<float>( pow(static_cast<float>(i + 1), 2) )* static_cast<float>(get_dT()) /static_cast<float>( 4.0) - static_cast<float> (r) * static_cast<float>((static_cast<float> (i) + static_cast<float>(1)) )* static_cast<float>( get_dT()) / static_cast<float>(4.0);
			b[i] = 1 - pow(sigma * (i + 1), 2) * get_dT() / 2;
			c[i] = static_cast<float> (pow(sigma, 2)) *static_cast<float>( pow(static_cast<float>(i + 1), 2)) *static_cast<float>( get_dT()) /static_cast<float>( 4.0 )+static_cast<float> (r)  *static_cast<float>( (i+1)) *static_cast<float>( get_dT() )/ static_cast<float>(4.0);
			d[i] = 1 + r * get_dT() + pow(sigma, 2) * pow(static_cast<float>(i + 1), 2) * get_dT() / 2;
		}
	}

	void Resolv::computeMatrix(float** M1, float** M2, float* c1, float* c2, float* c3, float* c4, int M) {
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				if (i == j) {
					M1[i][j] = c2[j];
					M2[i][j] = c4[j];
				}
				else if (j == i + 1) {
					M1[i][j] = c3[i];
					M2[i][j] = -c3[i];
				}
				else {
					M1[i][j] = 0;
					M2[i][j] = 0;
				}

			}
		}
		for (int j = 0; j < M - 1; j++) {
			M1[j + 1][j] = c1[j + 1];
			M2[j + 1][j] = -c1[j + 1];
		}

	}

	//Modele Differences finies


	void Resolv::computeMatrixDiffFinies(float** M1, float* c1, float* c2, int M) {
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				if (i == j) {
					M1[i][j] = c2[j];
				}
				else if (j == i + 1) {
					M1[i][j] = -c1[i];
				}
				else {
					M1[i][j] = 0;
				}

			}
		}
		for (int j = 0; j < M - 1; j++) {
			M1[j + 1][j] = -c1[j + 1];
		}
	}

	void Resolv::computeCoefDiffFinies(float* a, float* b, int M, float mu, float sigma) {
		for (int i = 0; i < M; i++) {
			a[i] = -static_cast<float>(mu) *static_cast<float>(get_dT_tilde(sigma)) / (static_cast<float>(static_cast<float>(get_dL()) * static_cast<float>(get_dL())));
			b[i] = static_cast<float>(1.0 )+ static_cast<float>(2.0) *static_cast<float>( mu) * static_cast<float>(get_dT_tilde(sigma) )/ static_cast<float>((get_dL()) *get_dL());
		}
	}



	//Fonctions générales 

	//fonction permetant de résoudre Ax=d avec A tridiagionale.
	float* Resolv::decompLU(float** &  A , float* d, int n) {
		//AX=B revient a resoudre Ly=f puis Ux=y L=lower U=upper
		float** lower = new float* [n];
		float** upper = new float* [n];
		float* y = new float[n];
		float* x = new float[n];
		float* z = new float[n];
		for (int j = 0; j < n; j++) {
			lower[j] = new float[n];
			upper[j] = new float[n];
		}
		/*
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				lower[i][j] = 0;
				upper[i][j] =0 ;
			}
		}
		lower[0][0] = A[0][0];
		for (int i = 0; i < n; i++) upper[i][i] = 1;
		double a = 0;
		a = A[0][1] / lower[0][0];
		for (int i = 1; i < n-1 ; i++) {
			lower[i][i - 1] = A[i][i - 1];
			lower[i][i] = A[i][i] - lower[i][i - 1] * upper[i - 1][i];
			upper[i][i + 1] = A[i][i + 1] / lower[i][i];
		}
		lower[n-1][n - 2] = A[n-1][n - 2];
		lower[n-1][n-1] = A[n-1][n-1] - lower[n-1][n - 2] * upper[n - 2][n-1];
		z[0] = d[0] / lower[0][0];
		for (int i = 1; i < n; i++) {
			z[i] = (d[i] - lower[i][i - 1] * z[i - 1]) / lower[i][i];
		}
		x[n-1] = z[n-1];
		for (int i = n - 2; i >= 0; i--) {
			x[i] = z[i] - upper[i][i + 1] * x[i + 1];
		}
		*/
		// Decomposing matrix into Upper and Lower
		// triangular matrix
		
		for (int i = 0; i < n; i++)
		{
			// Upper Triangular
			for (int k = i; k < n; k++)
			{
				// Summation of L(i, j) * U(j, k)
				float sum = 0;
				for (int j = 0; j < i; j++)
					sum += (lower[i][j] * upper[j][k]);

				// Evaluating U(i, k)
				upper[i][k] = A[i][k] - sum;
			}

			// Lower Triangular
			for (int k = i; k < n; k++)
			{
				if (i == k)
					lower[i][i] = 1; // Diagonal as 1
				else
				{
					// Summation of L(k, j) * U(j, i)
					float sum = 0;
					for (int j = 0; j < i; j++)
						sum += (lower[k][j] * upper[j][i]);

					// Evaluating L(k, i)
					lower[k][i]
						= (A[k][i] - sum) / upper[i][i];
				}
			}
		}

		float* beta = new float[n];
		float* alpha = new float[n];
		alpha[0] = A[0][0];
		for (int i = 1; i < n; i++) {
			beta[i] = A[i][i - 1] / alpha[i - 1];
			alpha[i] = A[i][i] - beta[i] * A[i - 1][i];
		}
		y[0] = d[0];
		for (int i = 1; i <= n - 1; i++) {
			y[i] = d[i] - beta[i] * y[i - 1];
		}
		x[n - 1] = y[n - 1] / alpha[n - 1];
		for (int i = n - 2; i >= 0; i--) {
			x[i] = (y[i] - A[i][i + 1] * x[i + 1]) / alpha[i];
		}
		delete[]alpha;
		delete[]beta;
		delete[] z;
		for (int i = 0; i < n; i++) {
			delete upper[i];
			delete lower[i];
		}
		delete[]lower;
		delete[]upper;
		return x;
	}

	float* Resolv::product(float** &  A, float* x, int n, int m) {
		float* tmp = new float[n];
		float* res = new float[n];
		for (int i = 0; i < n; i++) {
			tmp[i] = 0;
			for (int j = 0; j < m; j++) {
				tmp[i] += A[i][j] * x[j];
			}
			res[i] = tmp[i];
		}
		delete[] tmp;
		return res;
	}


	//************************************************************************************************************
	//PUT
	//************************************************************************************************************

	//constructeur du Put
	Put::Put(float T, float L, int N, int M, int id) : Resolv(T, L, N, M), id(id)
	{}


	//Crank Nicolson 
	//Pour un put ou un call, les coefficients des matrices sont les mêmes donc on reprend la fonction issue de la classe mère!
	void Put::computeCoef(float* a, float* b, float* c, float* d, int M, float sigma, float r) {
		Resolv::computeCoef(a, b, c, d, M, sigma, r);
	}

	//de meme, les matrices pour un Put ou un Call sont identiques donc on reprend la fonction de la classe mère
	void Put::computeMatrix(float** M1, float** M2, float* c1, float* c2, float* c3, float* c4, int M) {
		Resolv::computeMatrix(M1, M2, c1, c2, c3, c4, M);
	}

	//Permet de calculer à chaque étape le nouveau vecteur V_N en connaissant V_N+1
	void Put::calculate(float** &  M1, float**&  M2, float* V, int M, int N, float* c1, float* k, float* vector1, float K, float r, float T) {
		std::ofstream fichier;
		fichier.open("result_put.txt", std::ios::out);
		for (int i = N; i > 0; i--) {
			vector1 = product(M1, vector1, M, M);
			//utilisation des conditions initiales pour un put
			k[0] = c1[0] * (K * exp(-r * (T - (i - 1) * get_dT())) + K * exp(-r * (T - i * get_dT())));
			for (int s = 1; s < M; s++) {
				k[s] = 0.0;
			}
			for (int j = 0; j < M; j++) vector1[j] = vector1[j] + k[j];
			vector1 = Put::decompLU(M2, vector1, M);
			std::cout << (i - 1) * get_dT() << "|";
			for (int j = 0; j < M; j++) {
				std::cout << vector1[j] << ";";
				fichier << vector1[j] << ";";
			}
			std::cout << "0.00 " << std::endl;
			fichier << "\n" << std::endl;
		}
		fichier.close();
	}


	//permet d'afficher le stock et le time pour un put
	void Put::print(int M, int N, float* V) const {
		std::cout << "stock for put: ";
		for (int i = 0; i < M + 1; i++) {
			std::cout << i * get_dL() << "; ";
		}
		std::cout << std::endl;

		std::cout << "Time for put: ";
		for (int i = 0; i < M + 1; i++)
			std::cout << "-";
		std::cout << std::endl;


		std::cout << std::endl;
	}

	//Permet d'initialiser le vecteur V_M, avec les conditions initiales
	float* Put::initialize(int id, float* V, int M, float K) const {
		if (id == 0) {
			for (int i = 0; i < M; i++) {
				V[i] = max(0, K - i * get_dL());
			}
		}
		else {
			for (int i = 0; i < M; i++) {
				V[i] = max(0, i * get_dL() - K);
			}
		}
		return V;
	}

	//la decomposition LU est la même pour les 2 méthodes donc on reprend celle de la fonction mère
	float* Put::decompLU(float** &  A, float* d, int n) {
		return Resolv::decompLU(A, d, n);

	}

	//Le produit entre un vecteur et une matrice est le meme 
	float* Put::product(float** &  A, float* x, int n, int m) {
		return Resolv::product(A, x, n, m);
	}




	//differences finies
	void Put::computeCoefDiffFinies(float* a, float* b, int M, float mu, float sigma) {
		Resolv::computeCoefDiffFinies(a, b, M, mu, sigma);
	}

	void Put::computeMatrixDiffFinies(float** &  M1, float* c1, float* c2, int M) {
		Resolv::computeMatrixDiffFinies(M1, c1, c2, M);
	}

	//Permet d'initialiser le vecteur V_M, avec les conditions initiales
	float* Put::initializeDiffFinies(int id, float* V, int M, float K) const {
		if (id == 0) {
			for (int i = 0; i < M; i++) {
				V[i] = max(0, K - i * get_dL());
			}
		}
		else {
			for (int i = 0; i < M; i++) {
				V[i] = max(0, i * get_dL() - K);
			}
		}
		return V;
	}

	void Put::calculateDiffFinies(float**&   M1, float* V, int M, int N, float* c1, float* k, float* vector1, float K, float r, float T, float sigma) {
		std::ofstream fichier;
		fichier.open("result_put_DF.txt", std::ios::out);
		for (int i = N; i > 0; i--) {
			vector1 = product(M1, vector1, M, M);
			//utilisation des conditions initiales pour un put
			k[0] = c1[0] * (K * exp(-r * (static_cast<float>(T*sigma*sigma*0.5) - (static_cast<float>(i) - static_cast<float>(1)) *
				static_cast<float>(get_dT_tilde(sigma)))) + K * exp(-r * (static_cast<float>(T * sigma * sigma * 0.5) - static_cast<float>(i) * static_cast<float>(get_dT_tilde(sigma)))));
			for (int s = 1; s < M; s++) {
				k[s] = 0.0;
			}
			for (int j = 0; j < M; j++) vector1[j] = vector1[j] + k[j];
			std::cout << (i - 1) * get_dT() << "|";
			for (int j = 0; j < M; j++) {
				std::cout << vector1[j] << ";";
			}
			std::cout << std::endl;
			std::cout << std::endl;
			fichier << "\n" << std::endl;
		}
		for (int j = 0; j < M; j++) {
			fichier << vector1[j] << ";";
		}
		fichier << "\n" << std::endl;
		fichier.close();
	}

	//************************************************************************************************************
	//CALL
	//************************************************************************************************************

	//Contructeur pour un call
	Call::Call(float T, float L,int N, int M, int id) : Resolv(T, L,N, M), id(id)
	{}

	//Pour un put ou un call, les coefficients des matrices sont les mêmes donc on reprend la fonction issue de la classe mère!
	void Call::computeCoef(float* a, float* b, float* c, float* d, int M, float sigma, float r) {
		Resolv::computeCoef(a, b, c, d, M, sigma, r);
	}

	//de meme, les matrices pour un Put ou un Call sont identiques donc on reprend la fonction de la classe mère
	void Call::computeMatrix(float** M1, float** M2, float* c1, float* c2, float* c3, float* c4, int M) {
		Resolv::computeMatrix(M1, M2, c1, c2, c3, c4, M);
	}

	//permet d'afficher le stock et le time pour un call
	void Call::print(int M, int N, float* V) const {
		std::cout << "stock for call: ";
		for (int i = 0; i < M + 1; i++) {
			std::cout << i * (*this).get_dL() << "; ";
		}
		std::cout << std::endl;

		for (int i = 0; i < M + 1; i++)
			std::cout << "-";
		std::cout << std::endl;


		std::cout << std::endl;
	}

	//Permet d'initialiser le vecteur V_M, avec les conditions initiales 
	float* Call::initialize(int id, float* V, int M, float K) const {
		if (id == 0) {
			for (int i = 0; i < M+1; i++) {
				V[i] = max(0, K - i * get_dL());
			}
		}
		else if (id == 1) {
			for (int i = 0; i < M+1; i++) {
				V[i] = max(0, i * get_dL() - K);
			}
		}
		return V;
	}

	//la decomposition LU est la même pour les 2 méthodes donc on reprend celle de la fonction mère
	float* Call::decompLU(float** &  A, float* d, int n) {
		return Resolv::decompLU(A, d, n);

	}

	//Le produit entre un vecteur et une matrice est le meme 
	float* Call::product(float**&   A, float* x, int n, int m) {
		return Resolv::product(A, x, n, m);
	}

	//Permet de calculer à chaque étape le nouveau vecteur V_N en connaissant V_N+1 pour un call
	void Call::calculate(float** &  M1, float** &  M2, float* V, int M, int N, float* c1, float* k, float* vector1, float K, float r, float T) {
		std::ofstream fichier;
		fichier.open("result_call.txt", std::ios::out);
		for (int i = N; i > 0; i--) {
			vector1 = product(M1, vector1, M, M);
			//utilisation des conditions initiales pour un call
			k[M -1] = c1[M-1] * (K * exp(-r * (-T + i * get_dT())) + K * exp(-r * (-T + (i -1)* get_dT())));
			for (int s = 0; s < M-1 ; s++) {
				k[s] = 0.0;
			}
			for (int j = 0; j < M; j++) vector1[j] = vector1[j] + k[j];
			vector1 = decompLU(M2, vector1, M);
			std::cout << (i -1) * get_dT() << "|";
			for (int j = 0; j < M; j++) {
				std::cout << vector1[j] << ";";
				fichier << vector1[j] << ";";
			}
			std::cout << std::endl;
			std::cout << std::endl;
			fichier << "\n" << std::endl;
		}
		fichier.close();
	}

	//differences finies
	void Call::computeCoefDiffFinies(float* a, float* b, int M, float mu, float sigma) {
		Resolv::computeCoefDiffFinies(a, b, M, mu, sigma);
	}

	void Call::computeMatrixDiffFinies(float** M1, float* c1, float* c2, int M) {
		Resolv::computeMatrixDiffFinies(M1, c1, c2, M);
	}

	//Permet d'initialiser le vecteur V_M, avec les conditions initiales
	float* Call::initializeDiffFinies(int id, float* V, int M, float K) const {
		if (id == 0) {
			for (int i = 0; i < M+1; i++) {
				V[i] = max(0, K - i * get_dL());
			}
		}
		else {
			for (int i = 0; i < M+1; i++) {
				V[i] = max(0, i * get_dL() - K);
			}
		}
		return V;
	}

	void Call::calculateDiffFinies(float**  &  M1, float* V, int M, int N, float* c1, float* k, float* vector1, float K, float r, float T, float sigma) {
		std::ofstream fichier;
		fichier.open("result_call_DF.txt", std::ios::out);
		for (int i = N; i > 0; i--) {
			vector1 = product(M1, vector1, M, M);
			//utilisation des conditions initiales pour un call
			k[M - 1] = c1[M - 1] * (K * exp(-r * (-static_cast<float>(T * sigma * sigma * 0.5) + (static_cast<float>(i) - static_cast<float>(1)) *static_cast<float>( get_dT_tilde(sigma)))) + K * exp(-r * (-static_cast<float>(T * sigma * sigma * 0.5) + static_cast<float>(i) * static_cast<float>(get_dT_tilde(sigma)))));
			for (int s = 0; s < M - 1; s++) {
				k[s] = 0.0;
			}
			for (int j = 0; j < M; j++) vector1[j] = vector1[j] + k[j];
			std::cout << (i - 1) * get_dT() << "|";
			for (int j = 0; j < M; j++) {
				std::cout << vector1[j] << ";";
				fichier << vector1[j] << ";";
			}
			std::cout << std::endl;
			std::cout << std::endl;
			fichier << "\n" << std::endl;
		}
		fichier.close();
	}

}
