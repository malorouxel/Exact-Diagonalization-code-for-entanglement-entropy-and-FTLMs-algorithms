#include "sector.h"
#include "indexing.h"
#include "tables.h"
#include "parallel.h"
#include "common.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

// Mix the bits of a basis state. Use it as the hash function of the lookup table.
static inline unsigned long long Hash(unsigned long long x){
	x ^= x >> 33;
	x *= 0xff51afd7ed558ccdULL;
	x ^= x >> 33;
	x *= 0xc4ceb9fe1a85ec53ULL;
	x ^= x >> 33;
	return x;
}

// Build a hash table from basis state to index.
// It replaces the combinatorial ranking when one process holds the whole sector.
// It is faster than the ranking. It costs about four vectors of memory.
static void BuildHashTable(Sector *s){
	unsigned long long capacity = 1;
	while (capacity < 2 * s->size)
		capacity <<= 1;
	s->hash_mask   = capacity - 1;
	s->hash_keys   = malloc(capacity * sizeof(unsigned long long));
	s->hash_values = malloc(capacity * sizeof(unsigned long long));
	if (!s->hash_keys || !s->hash_values)
		Abort("Not enough memory for the hash table.");
	for (unsigned long long h = 0; h < capacity; h++)
		s->hash_keys[h] = HASH_EMPTY;
	for (unsigned long long i = 0; i < s->size; i++){
		unsigned long long h = Hash(s->states[i]) & s->hash_mask;
		while (s->hash_keys[h] != HASH_EMPTY)
			h = (h + 1) & s->hash_mask;
		s->hash_keys[h]   = s->states[i];
		s->hash_values[h] = i;
	}
}

// Find the basis indices that this process needs and does not own.
// These are the ghost indices. Their values come from other processes before each H|v>.
// The communication pattern does not change during a Lanczos run. We build it once.
static void BuildHalo(Sector *s){
	int nranks = NumberOfRanks();
	unsigned long long nwords = (s->size + 63) / 64;
	uint64_t *needed = calloc(nwords, sizeof(uint64_t));
	if (!needed)
		Abort("Not enough memory for the ghost bitmap.");

	unsigned long long end = s->start + s->nlocal;

	#pragma omp parallel for schedule(static)
	for (unsigned long long local = 0; local < s->nlocal; local++){
		unsigned long long states[MAX_ELEMENTS];
		double coeffs[MAX_ELEMENTS];
		int digits[64];
		unsigned long long state = SectorState(s, local);
		DecodeState(state, s->N, s->SpinOne, digits);
		for (int t = 0; t < s->H->n; t++){
			int n = TermElements(&s->H->terms[t], state, digits, s->SpinOne, states, coeffs);
			for (int e = 0; e < n; e++){
				if (states[e] == state)
					continue;
				unsigned long long j = SectorIndex(s, states[e]);
				if (j < s->start || j >= end){
					#pragma omp atomic
					needed[j / 64] |= (1ULL << (j % 64));
				}
			}
		}
	}

	s->nghost = 0;
	for (unsigned long long w = 0; w < nwords; w++)
		s->nghost += __builtin_popcountll(needed[w]);

	s->ghost_index  = malloc((s->nghost + 1) * sizeof(unsigned long long));
	s->ghost_values = malloc((s->nghost + 1) * sizeof(double));
	s->recvcounts   = calloc(nranks, sizeof(int));
	s->recvdispls   = calloc(nranks, sizeof(int));
	s->sendcounts   = calloc(nranks, sizeof(int));
	s->senddispls   = calloc(nranks, sizeof(int));

	unsigned long long g = 0;
	for (unsigned long long w = 0; w < nwords; w++){
		uint64_t word = needed[w];
		while (word){
			int bit = __builtin_ctzll(word);
			unsigned long long j = w * 64 + bit;
			s->ghost_index[g++] = j;
			s->recvcounts[SliceOwner(j, s->size, nranks)]++;
			word &= word - 1;
		}
	}
	free(needed);

	ExchangeCounts(s->recvcounts, s->sendcounts);
	s->nsend = 0;
	for (int p = 0; p < nranks; p++){
		s->recvdispls[p] = (p == 0) ? 0 : s->recvdispls[p - 1] + s->recvcounts[p - 1];
		s->senddispls[p] = (p == 0) ? 0 : s->senddispls[p - 1] + s->sendcounts[p - 1];
		s->nsend += s->sendcounts[p];
	}

	s->send_index  = malloc((s->nsend + 1) * sizeof(unsigned long long));
	s->send_values = malloc((s->nsend + 1) * sizeof(double));
	ExchangeIndices(s->ghost_index, s->recvcounts, s->recvdispls, s->send_index, s->sendcounts, s->senddispls);
	for (unsigned long long k = 0; k < s->nsend; k++)
		s->send_index[k] -= s->start;
}

// Prepare one magnetization sector.
// With distributed equal to 1 the basis is split between the processes.
// With distributed equal to 0 this process holds the whole sector.
void SectorBuild(Sector *s, int N, int Sz, int SzSym, int SpinOne, Interactions *H, int distributed){
	memset(s, 0, sizeof(Sector));
	s->N       = N;
	s->Sz      = Sz;
	s->SzSym   = SzSym;
	s->SpinOne = SpinOne;
	s->H       = H;
	s->size    = SectorSize(Sz, N, SpinOne, SzSym);

	int nranks = distributed ? NumberOfRanks() : 1;
	int rank   = distributed ? Rank() : 0;
	s->distributed = (nranks > 1);
	s->start  = SliceStart(s->size, rank, nranks);
	s->nlocal = SliceStart(s->size, rank + 1, nranks) - s->start;

	if (SzSym == 1){
		s->states = malloc((s->nlocal + 1) * sizeof(unsigned long long));
		if (!s->states)
			Abort("Not enough memory for the basis states.");
		#pragma omp parallel for schedule(static)
		for (unsigned long long local = 0; local < s->nlocal; local++)
			s->states[local] = IndexToState(s->start + local, N, Sz, SzSym, SpinOne);
		if (!s->distributed)
			BuildHashTable(s);
	}

	if (s->distributed)
		BuildHalo(s);
}

// Free the memory of a sector.
void SectorFree(Sector *s){
	free(s->states);
	free(s->hash_keys);
	free(s->hash_values);
	free(s->ghost_index);
	free(s->ghost_values);
	free(s->send_index);
	free(s->send_values);
	free(s->sendcounts);
	free(s->senddispls);
	free(s->recvcounts);
	free(s->recvdispls);
	memset(s, 0, sizeof(Sector));
}

// Return the basis state of a local index.
unsigned long long SectorState(const Sector *s, unsigned long long local){
	if (s->states)
		return s->states[local];
	return s->start + local;
}

// Return the global index of a basis state.
// Use the hash table when it exists. Use the combinatorial ranking without it.
unsigned long long SectorIndex(const Sector *s, unsigned long long state){
	if (s->SzSym == 0)
		return state;
	if (s->hash_keys){
		unsigned long long h = Hash(state) & s->hash_mask;
		while (s->hash_keys[h] != state)
			h = (h + 1) & s->hash_mask;
		return s->hash_values[h];
	}
	return StateToIndex(state, s->N, s->Sz, s->SzSym, s->SpinOne);
}

// Send the local values that other processes need and receive the ghost values.
// Call it before reading ghost values of a vector.
void SectorExchange(Sector *s, const double *in){
	if (!s->distributed)
		return;
	for (unsigned long long k = 0; k < s->nsend; k++)
		s->send_values[k] = in[s->send_index[k]];
	ExchangeValues(s->send_values, s->sendcounts, s->senddispls, s->ghost_values, s->recvcounts, s->recvdispls);
}

// Return the value of a vector at a global index.
// Read the local slice for an owned index. Search the ghost values for the other indices.
double SectorValue(const Sector *s, const double *in, unsigned long long index){
	if (index >= s->start && index - s->start < s->nlocal)
		return in[index - s->start];
	unsigned long long low = 0;
	unsigned long long high = s->nghost;
	while (high - low > 1){
		unsigned long long mid = (low + high) / 2;
		if (s->ghost_index[mid] <= index)
			low = mid;
		else
			high = mid;
	}
	return s->ghost_values[low];
}
