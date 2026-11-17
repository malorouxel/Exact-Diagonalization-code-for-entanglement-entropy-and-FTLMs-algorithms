#include "observables.h"
#include "hamiltonian.h"
#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Write the FTLM data of one random vector.
// Diagonalize the tridiagonal Lanczos matrix.
// Each line holds one Ritz value and the squared first component of its eigenvector.
// The file is folder/<prefix><Sz>_<r>.txt.
void WriteSpectralWeights(const LanczosRun *run, int Sz, int r, const char *folder, const char *prefix){
	int K = run->k;
	double *d = malloc(K * sizeof(double));
	double *e = malloc(K * sizeof(double));
	double *z = malloc((size_t)K * K * sizeof(double));
	memcpy(d, run->alphas, K * sizeof(double));
	memcpy(e, run->betas, K * sizeof(double));
	LAPACKE_dstevd(LAPACK_COL_MAJOR, 'V', K, d, e, z, K);

	char path[512];
	snprintf(path, sizeof(path), "%s/%s%d_%d.txt", folder, prefix, Sz, r);
	FILE *f = fopen(path, "w");
	for (int j = 0; j < K; j++){
		double w = z[(size_t)j * K];
		fprintf(f, "%.15f  %.15f\n", d[j], w * w);
	}
	fclose(f);
	free(d);
	free(e);
	free(z);
}

// Write the full spectrum of a sector solved by exact diagonalization.
// Each eigenvalue gets the weight 1/Size. The file has the same format as the FTLM files.
// The file is folder/<prefix><Sz>_exact.txt.
void WriteExactSpectrum(const double *EigVals, unsigned long long Size, int Sz, const char *folder, const char *prefix){
	char path[512];
	snprintf(path, sizeof(path), "%s/%s%d_exact.txt", folder, prefix, Sz);
	FILE *f = fopen(path, "w");
	for (unsigned long long i = 0; i < Size; i++)
		fprintf(f, "%.15f  %.15f\n", EigVals[i], 1.0 / Size);
	fclose(f);
}

// Diagonalize the full matrix of a sector held by one process.
// Use jobz equal to 'N' for the eigenvalues only and 'V' to also get the eigenvectors.
// The eigenvectors replace the matrix. They are stored by columns.
void ExactDiagonalization(Sector *s, double *EigVals, double *Matrix, char jobz){
	BuildDenseMatrix(s, Matrix);
	LAPACKE_dsyev(LAPACK_COL_MAJOR, jobz, 'U', (lapack_int)s->size, Matrix, (lapack_int)s->size, EigVals);
}
