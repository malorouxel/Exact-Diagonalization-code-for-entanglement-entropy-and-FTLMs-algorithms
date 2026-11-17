#ifndef INDEXING_H
#define INDEXING_H

#include "tables.h"

// Return the local state of one site.
// A basis state is a number written in base 2 for spin 1/2 and in base 3 for spin 1.
// Each digit is the local state of one site.
static inline int GetSpin(unsigned long long state, int site, int SpinOne){
	if (SpinOne == 1)
		return (int)((state / Pow3[site]) % 3);
	return (int)((state / Pow2[site]) % 2);
}

// Return the coefficient of S+ acting on a local state.
// Return 0 when the local state is already maximal.
static inline double SpCoeff(int value, int SpinOne){
	if (SpinOne == 1)
		return (value == 2) ? 0.0 : 1.4142135623730951;
	return (value == 1) ? 0.0 : 1.0;
}

// Return the coefficient of S- acting on a local state.
// Return 0 when the local state is already minimal.
static inline double SmCoeff(int value, int SpinOne){
	if (SpinOne == 1)
		return (value == 0) ? 0.0 : 1.4142135623730951;
	return (value == 0) ? 0.0 : 1.0;
}

// Return the eigenvalue of Sz for a local state.
static inline double SzValue(int value, int SpinOne){
	if (SpinOne == 1)
		return value - 1.0;
	return value - 0.5;
}


unsigned long long SpinHalfIndexToState(unsigned long long index, int Charge, int N);
unsigned long long SpinHalfStateToIndex(unsigned long long state, int N);
unsigned long long SpinOneIndexToState(unsigned long long index, int Charge, int N);
unsigned long long SpinOneStateToIndex(unsigned long long state, int Charge, int N);

unsigned long long StateToIndex(unsigned long long state, int N, int Sz, int SzSym, int SpinOne);
unsigned long long IndexToState(unsigned long long index, int N, int Sz, int SzSym, int SpinOne);

#endif
