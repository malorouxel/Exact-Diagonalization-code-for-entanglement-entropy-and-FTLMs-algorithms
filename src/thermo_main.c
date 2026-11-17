#include <stdio.h>
#include <stdlib.h>

#include "tables.h"
#include "thermodynamics.h"
#include "io.h"

// Compute the thermodynamics again from the files of a finished run.
// An optional number R uses only the first R random vectors.
// Usage. thermo [parameters file] [R]
int main(int argc, char **argv){
	const char *params_path = (argc > 1) ? argv[1] : "InputsParameters.txt";

	Parameters P;
	ReadParameters(params_path, &P);
	CreateTables(P.N);
	if (argc > 2 && atoi(argv[2]) > 0)
		P.R = atoi(argv[2]);

	if (P.ResultsFTLM || P.ResultsRFTLM || P.ResultsORFTLM)
		ComputeThermodynamics(&P);
	else
		fprintf(stderr, "No FTLM method in the Results line.\n");

	FreeTables(P.N);
	return 0;
}
