#include "logger.h"
#include <time.h>

// Return the wall clock time in seconds.
double WallTime(void){
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return t.tv_sec + t.tv_nsec * 1e-9;
}

// Write the duration of one step in the log file.
void LogStep(FILE *flog, int Sz, const char *name, double elapsed){
	if (!flog)
		return;
	fprintf(flog, "Sz = %d : %s : %.6f s\n", Sz, name, elapsed);
	fflush(flog);
}

// Write the memory estimate of one sector in the log file.
void LogRAM(FILE *flog, int Sz, double bytes){
	if (!flog)
		return;
	fprintf(flog, "Sz = %d : RAM estimation per process : %.3f MB\n", Sz, bytes / 1.0e6);
	fflush(flog);
}

// Estimate the memory of the vectors of one process in bytes.
// Add the dense matrix when the sector is solved by exact diagonalization.
double EstimationRAM(unsigned long long nlocal, int nvec, int dense){
	double bytes = (double)nvec * (double)nlocal * sizeof(double);
	if (dense)
		bytes += (double)nlocal * (double)nlocal * sizeof(double);
	return bytes;
}
