#include "tables.h"
#include <stdlib.h>

unsigned long long **Binomial;
unsigned long long *Pow2;
unsigned long long *Pow3;

// Fill the global tables for a chain of N sites.
// Binomial[k][n] is the binomial coefficient n choose k.
// Pow2 and Pow3 hold the powers of 2 and 3 up to N.
void CreateTables(int N){
	Binomial = malloc((2 * N + 1) * sizeof(unsigned long long *));
	for (int k = 0; k < 2 * N + 1; k++){
		Binomial[k] = malloc((N + 1) * sizeof(unsigned long long));
		for (int n = 0; n < N + 1; n++){
			if (k == 0)
				Binomial[k][n] = 1;
			else if (n < k)
				Binomial[k][n] = 0;
			else if (n == k)
				Binomial[k][n] = 1;
			else
				Binomial[k][n] = Binomial[k][n - 1] + Binomial[k - 1][n - 1];
		}
	}
	Pow2 = malloc((N + 1) * sizeof(unsigned long long));
	Pow3 = malloc((N + 1) * sizeof(unsigned long long));
	Pow2[0] = 1ULL;
	Pow3[0] = 1ULL;
	for (int k = 1; k <= N; k++){
		Pow2[k] = Pow2[k - 1] * 2ULL;
		Pow3[k] = Pow3[k - 1] * 3ULL;
	}
}

// Free the global tables.
void FreeTables(int N){
	for (int k = 0; k < 2 * N + 1; k++)
		free(Binomial[k]);
	free(Binomial);
	free(Pow2);
	free(Pow3);
}

// Return the dimension of a magnetization sector.
// Sz is the sum of the local states. Each local state is 0 or 1 for spin 1/2 and 0, 1 or 2 for spin 1.
// For spin 1 we count the sites in state 2 with k and place the remaining charge on the other sites.
// Without the Sz symmetry the sector is the full Hilbert space.
unsigned long long SectorSize(int Sz, int N, int SpinOne, int SzSym){
	unsigned long long Size = 0;
	if (SpinOne == 1){
		if (SzSym == 0)
			return Pow3[N];
		for (int k = 0; k < Sz / 2 + 1; k++){
			if (Sz - 2 * k > N - k)
				continue;
			Size += Binomial[k][N] * Binomial[Sz - 2 * k][N - k];
		}
	}
	else{
		if (SzSym == 0)
			return Pow2[N];
		Size = Binomial[Sz][N];
	}
	return Size;
}
