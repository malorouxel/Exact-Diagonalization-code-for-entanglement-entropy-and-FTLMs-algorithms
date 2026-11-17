#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

typedef enum {
	OP_SzSz,
	OP_SpSm,
	OP_SmSp,
	OP_Sz,
	OP_SS2
} OpType;

typedef struct {
	double J;
	OpType op;
	int site_i;
	int site_j;
} Term;

typedef struct {
	Term *terms;
	int n;
} Interactions;

struct Sector;

void PushTerm(Interactions *H, double J, OpType op, int site_i, int site_j);
void FreeInteractions(Interactions *H);
void DecodeState(unsigned long long state, int N, int SpinOne, int *digits);
int TermElements(const Term *t, unsigned long long state, const int *digits, int SpinOne, unsigned long long *states, double *coeffs);
void HamiltonianApply(struct Sector *s, const double *in, double *out);
void BuildDenseMatrix(struct Sector *s, double *Matrix);

#endif
