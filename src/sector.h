#ifndef SECTOR_H
#define SECTOR_H

#include "hamiltonian.h"
#include <limits.h>

#define HASH_EMPTY ULLONG_MAX

typedef struct Sector {
	int N;
	int Sz;
	int SzSym;
	int SpinOne;
	Interactions *H;

	unsigned long long size;
	unsigned long long start;
	unsigned long long nlocal;
	unsigned long long *states;

	unsigned long long *hash_keys;
	unsigned long long *hash_values;
	unsigned long long hash_mask;

	int distributed;
	unsigned long long nghost;
	unsigned long long *ghost_index;
	double *ghost_values;
	unsigned long long nsend;
	unsigned long long *send_index;
	double *send_values;
	int *sendcounts;
	int *senddispls;
	int *recvcounts;
	int *recvdispls;
} Sector;

void SectorBuild(Sector *s, int N, int Sz, int SzSym, int SpinOne, Interactions *H, int distributed);
void SectorFree(Sector *s);
unsigned long long SectorState(const Sector *s, unsigned long long local);
unsigned long long SectorIndex(const Sector *s, unsigned long long state);
void SectorExchange(Sector *s, const double *in);
double SectorValue(const Sector *s, const double *in, unsigned long long index);

#endif
