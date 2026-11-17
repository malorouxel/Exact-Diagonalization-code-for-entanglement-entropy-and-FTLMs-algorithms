#include "lanczos.h"
#include "hamiltonian.h"
#include "parallel.h"
#include "common.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// Allocate the Lanczos coefficients.
void LanczosRunInit(LanczosRun *run, int capacity){
	run->alphas   = malloc(capacity * sizeof(double));
	run->betas    = malloc(capacity * sizeof(double));
	run->k        = 0;
	run->capacity = capacity;
	run->seed     = 0;
}

// Copy the coefficients and the seed of a Lanczos run.
void LanczosRunCopy(LanczosRun *dst, const LanczosRun *src){
	memcpy(dst->alphas, src->alphas, src->k * sizeof(double));
	memcpy(dst->betas, src->betas, src->k * sizeof(double));
	dst->k    = src->k;
	dst->seed = src->seed;
}

// Free the Lanczos coefficients.
void LanczosRunFree(LanczosRun *run){
	free(run->alphas);
	free(run->betas);
	run->alphas = NULL;
	run->betas  = NULL;
}

// Return a pseudo random 64 bit number. This is the SplitMix64 generator.
static unsigned long long SplitMix(unsigned long long x){
	x += 0x9e3779b97f4a7c15ULL;
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
	x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
	return x ^ (x >> 31);
}

// Return a uniform number in ]0,1[ for a seed and a counter.
static double Uniform(unsigned long long seed, unsigned long long counter){
	unsigned long long x = SplitMix(SplitMix(seed) + counter);
	return ((x >> 11) + 0.5) / 9007199254740992.0;
}

// Return the dot product of two distributed vectors.
double Dot(Sector *s, const double *x, const double *y){
	return SumDouble(cblas_ddot((int)s->nlocal, x, 1, y, 1));
}

// Fill the local slice of a normalized Gaussian random vector.
// Each component only depends on the seed and its global index.
// The vector does not depend on the number of processes.
void RandomVector(Sector *s, unsigned long long seed, double *v){
	for (unsigned long long local = 0; local < s->nlocal; local++){
		unsigned long long g = s->start + local;
		double u1 = Uniform(seed, 2 * g);
		double u2 = Uniform(seed, 2 * g + 1);
		v[local] = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
	}
	double norm = sqrt(Dot(s, v, v));
	cblas_dscal((int)s->nlocal, 1.0 / norm, v, 1);
}

// Remove from v its components along the orthogonal vectors.
static void Project(Sector *s, double **OrthoVecs, int northo, double *v){
	for (int d = 0; d < northo; d++){
		double projection = Dot(s, v, OrthoVecs[d]);
		cblas_daxpy((int)s->nlocal, -projection, OrthoVecs[d], 1, v, 1);
	}
}

// Return the lowest eigenvalue of the tridiagonal Lanczos matrix.
double LowestRitzValue(const LanczosRun *run){
	double *d = malloc(run->k * sizeof(double));
	double *e = malloc(run->k * sizeof(double));
	memcpy(d, run->alphas, run->k * sizeof(double));
	memcpy(e, run->betas, run->k * sizeof(double));
	LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', run->k, d, e, NULL, 1);
	double E = d[0];
	free(d);
	free(e);
	return E;
}

// Run the Lanczos recursion from the random vector of run->seed.
// Store the coefficients alphas and betas in run.
// Keep every Lanczos vector orthogonal to OrthoVecs when northo is not zero.
// Stop at NIterationMax steps or when the lowest Ritz value converges if LanczosConv is 1.
// Stop early when the Krylov space is invariant. Return the lowest Ritz value.
static double LanczosRecursion(Sector *s, LanczosRun *run, double **OrthoVecs, int northo, int LanczosConv, int NIterationMax){
	int n = (int)s->nlocal;
	double *work = malloc(3 * (n + 1) * sizeof(double));
	double *V1 = work;
	double *V2 = work + (n + 1);
	double *W  = work + 2 * (n + 1);
	if (NIterationMax > run->capacity)
		NIterationMax = run->capacity;

	RandomVector(s, run->seed, V1);
	if (northo > 0){
		Project(s, OrthoVecs, northo, V1);
		double norm = sqrt(Dot(s, V1, V1));
		if (norm < EPSILON)
			Abort("The random start vector lies in the orthogonal space.");
		cblas_dscal(n, 1.0 / norm, V1, 1);
	}

	HamiltonianApply(s, V1, V2);
	double alpha = Dot(s, V1, V2);
	run->alphas[0] = alpha;
	cblas_daxpy(n, -alpha, V1, 1, V2, 1);

	int k = 1;
	double E = 1e300;
	int converged = 0;
	while (!converged && k < NIterationMax){
		if (northo > 0)
			Project(s, OrthoVecs, northo, V2);
		double beta = sqrt(Dot(s, V2, V2));
		if (beta < EPSILON)
			break;
		run->betas[k - 1] = beta;
		cblas_dscal(n, 1.0 / beta, V2, 1);

		HamiltonianApply(s, V2, W);
		cblas_daxpy(n, -beta, V1, 1, W, 1);
		alpha = Dot(s, V2, W);
		run->alphas[k] = alpha;
		cblas_daxpy(n, -alpha, V2, 1, W, 1);

		double *tmp = V1;
		V1 = V2;
		V2 = W;
		W  = tmp;
		k++;

		if (LanczosConv == 1){
			run->k = k;
			double Enew = LowestRitzValue(run);
			if (fabs(E - Enew) < CONV_THRESHOLD)
				converged = 1;
			E = Enew;
		}
	}
	run->k = k;
	free(work);
	return LowestRitzValue(run);
}

// Build the lowest Ritz vector of a finished Lanczos run.
// Run the recursion again from the same start vector with the stored coefficients.
// Add the Lanczos vectors weighted by the lowest eigenvector of the tridiagonal matrix.
void OrthogonalEigenvector(Sector *s, double **OrthoVecs, int northo, const LanczosRun *run, double *out){
	int n = (int)s->nlocal;
	int K = run->k;

	double *d = malloc(K * sizeof(double));
	double *e = malloc(K * sizeof(double));
	double *z = malloc((size_t)K * K * sizeof(double));
	memcpy(d, run->alphas, K * sizeof(double));
	memcpy(e, run->betas, K * sizeof(double));
	LAPACKE_dstevd(LAPACK_COL_MAJOR, 'V', K, d, e, z, K);

	double *work = malloc(3 * (n + 1) * sizeof(double));
	double *V1 = work;
	double *V2 = work + (n + 1);
	double *W  = work + 2 * (n + 1);

	RandomVector(s, run->seed, V1);
	if (northo > 0){
		Project(s, OrthoVecs, northo, V1);
		double norm = sqrt(Dot(s, V1, V1));
		cblas_dscal(n, 1.0 / norm, V1, 1);
	}
	for (int i = 0; i < n; i++)
		out[i] = z[0] * V1[i];

	HamiltonianApply(s, V1, V2);
	cblas_daxpy(n, -run->alphas[0], V1, 1, V2, 1);

	for (int k = 1; k < K; k++){
		if (northo > 0)
			Project(s, OrthoVecs, northo, V2);
		cblas_dscal(n, 1.0 / run->betas[k - 1], V2, 1);
		cblas_daxpy(n, z[k], V2, 1, out, 1);
		if (k == K - 1)
			break;

		HamiltonianApply(s, V2, W);
		cblas_daxpy(n, -run->betas[k - 1], V1, 1, W, 1);
		cblas_daxpy(n, -run->alphas[k], V2, 1, W, 1);

		double *tmp = V1;
		V1 = V2;
		V2 = W;
		W  = tmp;
	}

	double norm = sqrt(Dot(s, out, out));
	cblas_dscal(n, 1.0 / norm, out, 1);

	free(d);
	free(e);
	free(z);
	free(work);
}

// Run a Lanczos recursion from a new random vector. Return the lowest Ritz value.
// The seed is increased by one at each call.
double FirstLanczos(Sector *s, unsigned long long *seed, LanczosRun *run, int LanczosConv, int NIterationMax){
	(*seed)++;
	run->seed = *seed;
	return LanczosRecursion(s, run, NULL, 0, LanczosConv, NIterationMax);
}

// Build the ground state from the coefficients of a previous FirstLanczos run.
void SecondLanczos(Sector *s, const LanczosRun *run, double *GroundState){
	OrthogonalEigenvector(s, NULL, 0, run, GroundState);
}

// Run a Lanczos recursion in the space orthogonal to OrthoVecs. Return the lowest Ritz value.
double OrthogonalLanczos(Sector *s, double **OrthoVecs, int northo, unsigned long long *seed, LanczosRun *run, int LanczosConv, int NIterationMax){
	(*seed)++;
	run->seed = *seed;
	return LanczosRecursion(s, run, OrthoVecs, northo, LanczosConv, NIterationMax);
}

// Find all the states of the sector with energy Eref.
// OrthoVecs holds one ground state on entry.
// Search the lowest state orthogonal to the known ones. Keep it if its energy is Eref.
// Stop at the first state with a higher energy.
// Return the degeneracy.
int FindDegenerateStates(Sector *s, double Eref, double ***OrthoVecs, unsigned long long *seed, int LanczosConv, int NIterationMax){
	int degen = 1;
	LanczosRun run;
	LanczosRunInit(&run, NIterationMax);
	while ((unsigned long long)degen < s->size){
		double E = OrthogonalLanczos(s, *OrthoVecs, degen, seed, &run, LanczosConv, NIterationMax);
		if (fabs(Eref - E) >= CONV_THRESHOLD)
			break;
		*OrthoVecs = realloc(*OrthoVecs, (degen + 1) * sizeof(double *));
		(*OrthoVecs)[degen] = malloc((s->nlocal + 1) * sizeof(double));
		OrthogonalEigenvector(s, *OrthoVecs, degen, &run, (*OrthoVecs)[degen]);
		degen++;
	}
	LanczosRunFree(&run);
	return degen;
}
