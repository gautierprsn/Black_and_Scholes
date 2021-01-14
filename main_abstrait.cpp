#include "methode_abstrait.h"
#include <thread>
#include <iostream>

int main() {
	float T = 1.0;
	float r =static_cast<float>( 0.1);
	float L = 300.0;
	float sd = static_cast<float>(0.1);
	float K = 100.0;
	int N = 500;
	int M = 500;
	int demand;
	float mu = sd * sd / 2;

	
	std::cout << "Voulez vous calculer un call (rentrez 1) ou un put (rentrez 0)?" << std::endl;
	std::cin >> demand;
	
	

	if (demand == 0) {
		int demand2;
		std::cout << "Voulez vous calculer via Crank-Nicholson (rentrez 1) ou par la méthode des différences finies (rentrez 0)?" << std::endl;
		std::cin >> demand2;
		int id = 0;
		if (demand2 == 1) {


			ensiie::Put put = ensiie::Put::Put(T, L, N, M, id);

			float* c1 = new float[M];
			float* c2 = new float[M];
			float* c3 = new float[M];
			float* c4 = new float[M];


			put.ensiie::Put::computeCoef(c1, c2, c3, c4, M, sd, r);


			float** M1 = new float* [M];
			float** M2 = new float* [M];
			for (int j = 0; j < M; j++) {
				M1[j] = new float[M];
				M2[j] = new float[M];
			}

			put.computeMatrix(M1, M2, c1, c2, c3, c4, M);

			float* V = new float[M];
			float* k = new float[M];

			put.print(M, N, V);

			float* vector1 = new float[M];

			vector1 = put.initialize(put.get_id(), vector1, M, K);

			put.calculate(M1, M2, V, M, N, c1, k, vector1, K, r, T);

			put.deleted(c1, c2, c3, c4, M1, M2, M, V, k);
		}
		else {
			ensiie::Put put = ensiie::Put::Put(T, L, N, M, id);

			float* c1 = new float[M];
			float* c2 = new float[M];


			put.ensiie::Put::computeCoefDiffFinies(c1, c2, M, mu, sd);


			float** M1 = new float* [M];
			float** M2 = new float* [M];
			for (int j = 0; j < M; j++) {
				M1[j] = new float[M];
				M2[j] = new float[M];
			}

			put.computeMatrixDiffFinies(M1, c1, c2, M);

			float* V = new float[M];
			float* k = new float[M];

			put.print(M, N, V);

			float* vector1 = new float[M];

			vector1 = put.initializeDiffFinies(put.get_id(), vector1, M, K);

			put.calculateDiffFinies(M1, V, M, N, c1, k, vector1, K, r, T,sd);

			put.deletedDiffFinies(c1, c2, M1, M2, M, V, k);

		}

	}
	if (demand == 1) {
		int demand2;
		std::cout << "Voulez vous calculer via Crank-Nicholson (rentrez 1) ou par la méthode des différences finies (rentrez 0)?" << std::endl;
		std::cin >> demand2;
		int id = 1;
		if (demand2 == 1) {
			ensiie::Call call = ensiie::Call::Call(T, L, N, M, id);

			float* c1 = new float[M];
			float* c2 = new float[M];
			float* c3 = new float[M];
			float* c4 = new float[M];


			call.ensiie::Call::computeCoef(c1, c2, c3, c4, M, sd, r);


			float** M1 = new float* [M];
			float** M2 = new float* [M];
			for (int j = 0; j < M; j++) {
				M1[j] = new float[M];
				M2[j] = new float[M];
			}

			call.computeMatrix(M1, M2, c1, c2, c3, c4, M);
			
			float* V = new float[M];
			float* k = new float[M];

			call.print(M, N, V);

			float* vector1 = new float[M];

			vector1 = call.initialize(call.get_id(), vector1, M, K);

			call.calculate(M1, M2, V, M, N, c1, k, vector1, K, r, T);


			call.deleted(c1, c2, c3, c4, M1, M2, M, V, k);
		}
		if (demand2 == 0) {
			ensiie::Call call= ensiie::Call::Call(T, L, N, M, id);

			float* c1 = new float[M];
			float* c2 = new float[M];


			call.ensiie::Call::computeCoefDiffFinies(c1, c2, M, mu, sd);


			float** M1 = new float* [M];
			float** M2 = new float* [M];
			for (int j = 0; j < M; j++) {
				M1[j] = new float[M];
				M2[j] = new float[M];
			}

			call.computeMatrixDiffFinies(M1, c1, c2, M);

			float* V = new float[M];
			float* k = new float[M];

			call.print(M, N, V);

			float* vector1 = new float[M];

			vector1 = call.initializeDiffFinies(call.get_id(), vector1, M, K);

			call.calculateDiffFinies(M1, V, M, N, c1, k, vector1, K, r, T, sd);

			call.deletedDiffFinies(c1, c2, M1, M2, M, V, k);

		}
	}
	else { std::cout << "Rentrez 0 ou 1" << std::endl; }

	return 0; 
}