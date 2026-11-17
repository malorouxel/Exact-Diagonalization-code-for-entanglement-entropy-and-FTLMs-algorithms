#include "hamiltonian.h"
#include "sector.h"
#include "indexing.h"
#include "common.h"
#include "tables.h"
#include <stdlib.h>
#include <string.h>

// Add one elementary term to the Hamiltonian.
void PushTerm(Interactions *H, double J, OpType op, int site_i, int site_j){
	H->terms = realloc(H->terms, (H->n + 1) * sizeof(Term));
	H->terms[H->n] = (Term){ .J = J, .op = op, .site_i = site_i, .site_j = site_j };
	H->n++;
}

// Free the list of terms.
void FreeInteractions(Interactions *H){
	free(H->terms);
	H->terms = NULL;
	H->n = 0;
}

// Move the spin on site i by di and the spin on site j by dj.
// The new state is found by adding powers of the base. No division is needed.
static unsigned long long Flip(unsigned long long state, int si, int sj, int di, int dj, int SpinOne){
	const unsigned long long *Pow = (SpinOne == 1) ? Pow3 : Pow2;
	return state + (unsigned long long)di * Pow[si] + (unsigned long long)dj * Pow[sj];
}

// Write the local state of every site of a basis state in digits.
void DecodeState(unsigned long long state, int N, int SpinOne, int *digits){
	int base = SpinOne + 2;
	for (int site = 0; site < N; site++){
		digits[site] = (int)(state % base);
		state /= base;
	}
}

// Apply S-_i S+_j or S+_i S-_j to a state and add the terms of (S_i.S_j)^2 that start with it.
// sign is -1 for S-_i S+_j and +1 for S+_i S-_j.
static int SS2Branch(const Term *t, unsigned long long state, int vi, int vj, int sign, int SpinOne, unsigned long long *states, double *coeffs, int n){
	int si = t->site_i;
	int sj = t->site_j;
	double m1 = SzValue(vi, SpinOne);
	double m2 = SzValue(vj, SpinOne);
	double ci = (sign < 0) ? SmCoeff(vi, SpinOne) : SpCoeff(vi, SpinOne);
	double cj = (sign < 0) ? SpCoeff(vj, SpinOne) : SmCoeff(vj, SpinOne);
	if (ci == 0.0 || cj == 0.0)
		return n;

	unsigned long long state2 = Flip(state, si, sj, sign, -sign, SpinOne);
	int vi2 = vi + sign;
	int vj2 = vj - sign;
	double nm1 = SzValue(vi2, SpinOne);
	double nm2 = SzValue(vj2, SpinOne);
	states[n] = state2;
	coeffs[n] = t->J * 0.5 * (m1 * m2 + nm1 * nm2) * ci * cj;
	n++;

	double a = SmCoeff(vi2, SpinOne);
	double b = SpCoeff(vj2, SpinOne);
	if (a != 0.0 && b != 0.0){
		states[n] = Flip(state2, si, sj, -1, 1, SpinOne);
		coeffs[n] = t->J * 0.25 * a * b * ci * cj;
		n++;
	}
	a = SpCoeff(vi2, SpinOne);
	b = SmCoeff(vj2, SpinOne);
	if (a != 0.0 && b != 0.0){
		states[n] = Flip(state2, si, sj, 1, -1, SpinOne);
		coeffs[n] = t->J * 0.25 * a * b * ci * cj;
		n++;
	}
	return n;
}

// Apply one term to a basis state.
// digits holds the local states of the basis state.
// Store the reached states and the matrix elements <state'|term|state>.
// Return the number of elements. It is never larger than MAX_ELEMENTS.
int TermElements(const Term *t, unsigned long long state, const int *digits, int SpinOne, unsigned long long *states, double *coeffs){
	int si = t->site_i;
	int sj = t->site_j;
	int n = 0;

	if (t->op == OP_Sz){
		states[0] = state;
		coeffs[0] = t->J * SzValue(digits[si], SpinOne);
		return 1;
	}

	int vi = digits[si];
	int vj = digits[sj];
	double m1 = SzValue(vi, SpinOne);
	double m2 = SzValue(vj, SpinOne);

	if (t->op == OP_SzSz){
		states[0] = state;
		coeffs[0] = t->J * m1 * m2;
		n = 1;
	}
	else if (t->op == OP_SpSm){
		double ci = SpCoeff(vi, SpinOne);
		double cj = SmCoeff(vj, SpinOne);
		if (ci != 0.0 && cj != 0.0){
			states[0] = Flip(state, si, sj, 1, -1, SpinOne);
			coeffs[0] = t->J * ci * cj;
			n = 1;
		}
	}
	else if (t->op == OP_SmSp){
		double ci = SmCoeff(vi, SpinOne);
		double cj = SpCoeff(vj, SpinOne);
		if (ci != 0.0 && cj != 0.0){
			states[0] = Flip(state, si, sj, -1, 1, SpinOne);
			coeffs[0] = t->J * ci * cj;
			n = 1;
		}
	}
	else if (t->op == OP_SS2){
		states[0] = state;
		coeffs[0] = t->J * m1 * m1 * m2 * m2;
		n = 1;
		n = SS2Branch(t, state, vi, vj, -1, SpinOne, states, coeffs, n);
		n = SS2Branch(t, state, vi, vj, 1, SpinOne, states, coeffs, n);
	}
	return n;
}

// Compute out = H in on the local slice of the vectors.
// Each process computes its own rows. Row i is the sum of <j|H|i> in_j over the states j reached from i.
// This uses the symmetry of H. The values in_j owned by other processes are the ghost values.
void HamiltonianApply(Sector *s, const double *in, double *out){
	SectorExchange(s, in);

	#pragma omp parallel for schedule(static)
	for (unsigned long long local = 0; local < s->nlocal; local++){
		unsigned long long states[MAX_ELEMENTS];
		double coeffs[MAX_ELEMENTS];
		int digits[64];
		unsigned long long state = SectorState(s, local);
		DecodeState(state, s->N, s->SpinOne, digits);
		double sum = 0.0;
		for (int t = 0; t < s->H->n; t++){
			int n = TermElements(&s->H->terms[t], state, digits, s->SpinOne, states, coeffs);
			for (int e = 0; e < n; e++){
				if (states[e] == state)
					sum += coeffs[e] * in[local];
				else
					sum += coeffs[e] * SectorValue(s, in, SectorIndex(s, states[e]));
			}
		}
		out[local] = sum;
	}
}

// Fill the full matrix of H in a sector held by one process.
// The matrix is stored by columns.
void BuildDenseMatrix(Sector *s, double *Matrix){
	unsigned long long Size = s->size;
	memset(Matrix, 0, Size * Size * sizeof(double));
	for (unsigned long long i = 0; i < Size; i++){
		unsigned long long states[MAX_ELEMENTS];
		double coeffs[MAX_ELEMENTS];
		int digits[64];
		unsigned long long state = SectorState(s, i);
		DecodeState(state, s->N, s->SpinOne, digits);
		for (int t = 0; t < s->H->n; t++){
			int n = TermElements(&s->H->terms[t], state, digits, s->SpinOne, states, coeffs);
			for (int e = 0; e < n; e++)
				Matrix[SectorIndex(s, states[e]) + i * Size] += coeffs[e];
		}
	}
}
