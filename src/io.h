#ifndef IO_H
#define IO_H

#include "hamiltonian.h"

#define HEADER_BYTES 32

typedef struct {
	int SpinOne;
	int N;
	unsigned long long EDThreshold;
	int SzSym;
	int R;
	int LanczosConv;
	int NIterationMax;
	int KnownDege;
	int Dege;
	int BipartitionSize;
	int ResultsE0;
	int ResultsGS;
	int ResultsDege;
	int ResultsEntropy;
	int ResultsORFTLM;
	int ResultsRFTLM;
	int ResultsFTLM;
	double logTmin;
	double logTmax;
	int NT;
	int RandomSeed;
	unsigned long long Seed;
	char FolderName[256];
} Parameters;

typedef struct {
	int N;
	int SpinOne;
	int SzSym;
	int Sz;
	unsigned long long Size;
	double E0;
} GroundStateHeader;

void ReadParameters(const char *path, Parameters *P);
void ReadInteractions(const char *path, Interactions *H);
void WriteGroundState(const char *path, const GroundStateHeader *header, const double *values, unsigned long long start, unsigned long long count);
void ReadGroundStateHeader(const char *path, GroundStateHeader *header);

#endif
