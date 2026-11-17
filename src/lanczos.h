#ifndef LANCZOS_H
#define LANCZOS_H

#include "sector.h"

typedef struct {
	double *alphas;
	double *betas;
	int k;
	int capacity;
	unsigned long long seed;
} LanczosRun;

void LanczosRunInit(LanczosRun *run, int capacity);
void LanczosRunCopy(LanczosRun *dst, const LanczosRun *src);
void LanczosRunFree(LanczosRun *run);

void RandomVector(Sector *s, unsigned long long seed, double *v);
double Dot(Sector *s, const double *x, const double *y);
double LowestRitzValue(const LanczosRun *run);

double FirstLanczos(Sector *s, unsigned long long *seed, LanczosRun *run, int LanczosConv, int NIterationMax);
void SecondLanczos(Sector *s, const LanczosRun *run, double *GroundState);
double OrthogonalLanczos(Sector *s, double **OrthoVecs, int northo, unsigned long long *seed, LanczosRun *run, int LanczosConv, int NIterationMax);
void OrthogonalEigenvector(Sector *s, double **OrthoVecs, int northo, const LanczosRun *run, double *out);
int FindDegenerateStates(Sector *s, double Eref, double ***OrthoVecs, unsigned long long *seed, int LanczosConv, int NIterationMax);

#endif
