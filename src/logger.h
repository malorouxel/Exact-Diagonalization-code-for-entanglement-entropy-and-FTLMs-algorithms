#ifndef LOGGER_H
#define LOGGER_H

#include <stdio.h>

double WallTime(void);
void LogStep(FILE *flog, int Sz, const char *name, double elapsed);
void LogRAM(FILE *flog, int Sz, double bytes);
double EstimationRAM(unsigned long long nlocal, int nvec, int dense);

#endif
