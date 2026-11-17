#include "entropy.h"
#include "io.h"
#include "indexing.h"
#include "tables.h"
#include "parallel.h"
#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define READ_CHUNK (1ULL << 20)

// Return the von Neumann entropy of one block from its singular values.
static double BlockEntropy(double *Matrix, unsigned long long rows, unsigned long long cols){
	unsigned long long nsv = (rows < cols) ? rows : cols;
	double *S = malloc(nsv * sizeof(double));
	double *superb = malloc((nsv > 1 ? nsv : 2) * sizeof(double));
	LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'N', 'N', (lapack_int)rows, (lapack_int)cols, Matrix, (lapack_int)rows, S, NULL, 1, NULL, 1, superb);
	double entropy = 0.0;
	for (unsigned long long i = 0; i < nsv; i++){
		double p = S[i] * S[i];
		if (p > 1e-15)
			entropy -= p * log(p);
	}
	free(S);
	free(superb);
	return entropy;
}

// Return the entanglement entropy of a ground state file for the subsystem A.
// SitesA lists the sites of A. Sites start at 0 and do not need to be contiguous.
// The ground state lives in one Sz sector. The reduced density matrix is block diagonal.
// Each block has a fixed charge qA in A and Sz - qA in the rest B.
// We build the matrix psi(a,b) of each block and add the entropies of their singular values.
// Without the Sz symmetry there is only one block of size d^NA x d^NB.
// The blocks are shared between the processes. Every process calls this function.
double EntanglementEntropy(const char *path, const int *SitesA, int SizeA){
	GroundStateHeader h;
	ReadGroundStateHeader(path, &h);
	int N = h.N;
	int SizeB = N - SizeA;
	int maxlocal = h.SpinOne + 1;
	unsigned long long *Pow = (h.SpinOne == 1) ? Pow3 : Pow2;

	int *inA = calloc(N, sizeof(int));
	for (int k = 0; k < SizeA; k++){
		if (SitesA[k] < 0 || SitesA[k] >= N || inA[SitesA[k]])
			Abort("The subsystem A has a wrong or repeated site.");
		inA[SitesA[k]] = 1;
	}
	if (SizeA == 0 || SizeB == 0){
		free(inA);
		return 0.0;
	}

	int qmin = 0;
	int nblocks = 1;
	if (h.SzSym == 1){
		qmin = h.Sz - SizeB * maxlocal;
		if (qmin < 0) qmin = 0;
		int qmax = (h.Sz < SizeA * maxlocal) ? h.Sz : SizeA * maxlocal;
		nblocks = qmax - qmin + 1;
	}

	int rank = Rank();
	int nranks = NumberOfRanks();
	unsigned long long *rows = malloc(nblocks * sizeof(unsigned long long));
	unsigned long long *cols = malloc(nblocks * sizeof(unsigned long long));
	double **blocks = calloc(nblocks, sizeof(double *));
	for (int b = 0; b < nblocks; b++){
		if (h.SzSym == 1){
			rows[b] = SectorSize(qmin + b, SizeA, h.SpinOne, 1);
			cols[b] = SectorSize(h.Sz - qmin - b, SizeB, h.SpinOne, 1);
		}
		else{
			rows[b] = Pow[SizeA];
			cols[b] = Pow[SizeB];
		}
		if (b % nranks == rank){
			blocks[b] = calloc(rows[b] * cols[b], sizeof(double));
			if (!blocks[b])
				Abort("Not enough memory for a Schmidt block.");
		}
	}

	FILE *f = fopen(path, "rb");
	fseek(f, HEADER_BYTES, SEEK_SET);
	double *buffer = malloc(READ_CHUNK * sizeof(double));
	unsigned long long done = 0;
	while (done < h.Size){
		unsigned long long count = h.Size - done;
		if (count > READ_CHUNK) count = READ_CHUNK;
		if (fread(buffer, sizeof(double), count, f) != count)
			Abort("The ground state file is too short.");
		for (unsigned long long k = 0; k < count; k++){
			unsigned long long state = IndexToState(done + k, N, h.Sz, h.SzSym, h.SpinOne);
			unsigned long long stateA = 0;
			unsigned long long stateB = 0;
			int posA = 0;
			int posB = 0;
			int qA = 0;
			for (int site = 0; site < N; site++){
				int digit = GetSpin(state, site, h.SpinOne);
				if (inA[site]){
					stateA += digit * Pow[posA++];
					qA += digit;
				}
				else{
					stateB += digit * Pow[posB++];
				}
			}
			int b = (h.SzSym == 1) ? qA - qmin : 0;
			if (!blocks[b])
				continue;
			unsigned long long row = StateToIndex(stateA, SizeA, qA, h.SzSym, h.SpinOne);
			unsigned long long col = StateToIndex(stateB, SizeB, h.Sz - qA, h.SzSym, h.SpinOne);
			blocks[b][row + col * rows[b]] = buffer[k];
		}
		done += count;
	}
	fclose(f);
	free(buffer);

	double entropy = 0.0;
	for (int b = 0; b < nblocks; b++){
		if (!blocks[b])
			continue;
		entropy += BlockEntropy(blocks[b], rows[b], cols[b]);
		free(blocks[b]);
	}
	entropy = SumDouble(entropy);

	free(blocks);
	free(rows);
	free(cols);
	free(inA);
	return entropy;
}
