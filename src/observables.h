#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include "sector.h"
#include "lanczos.h"

void WriteSpectralWeights(const LanczosRun *run, int Sz, int r, const char *folder, const char *prefix);
void WriteExactSpectrum(const double *EigVals, unsigned long long Size, int Sz, const char *folder, const char *prefix);
void ExactDiagonalization(Sector *s, double *EigVals, double *Matrix, char jobz);

#endif
