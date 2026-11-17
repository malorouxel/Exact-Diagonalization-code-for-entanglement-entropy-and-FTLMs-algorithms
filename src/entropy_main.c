#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parallel.h"
#include "tables.h"
#include "entropy.h"
#include "io.h"

// Read one subsystem from a command line argument.
// The argument is a list of sites separated by commas. A range a-b adds every site from a to b.
// Example. 0,2,4-6 gives the sites 0 2 4 5 6.
// Return the number of sites.
static int ParseSites(const char *arg, int *sites, int N){
	int n = 0;
	char *copy = strdup(arg);
	for (char *token = strtok(copy, ","); token; token = strtok(NULL, ",")){
		int a, b;
		if (sscanf(token, "%d-%d", &a, &b) != 2)
			b = a = atoi(token);
		for (int k = a; k <= b; k++){
			if (n >= N)
				Abort("The subsystem has more sites than the chain.");
			sites[n++] = k;
		}
	}
	free(copy);
	return n;
}

// Compute the entanglement entropy of a saved ground state for several subsystems.
// Usage. entropy <ground state file> <sites of A> [<sites of A> ...]
// Each subsystem gives one line on the standard output.
int main(int argc, char **argv){
	ParallelInit(&argc, &argv);
	if (argc < 3){
		if (Rank() == 0)
			fprintf(stderr, "Usage. %s <ground state file> <sites of A> [<sites of A> ...]\nExample. %s GroundState.bin 0-4 0,2,4,6\n", argv[0], argv[0]);
		ParallelFinalize();
		return 1;
	}

	GroundStateHeader h;
	ReadGroundStateHeader(argv[1], &h);
	CreateTables(h.N);
	if (Rank() == 0)
		fprintf(stderr, "N = %d : Spin = %s : Sz sector = %d : Size = %llu : E0 = %.15f\n", h.N, h.SpinOne ? "One" : "Half", h.Sz, h.Size, h.E0);

	int *sites = malloc(h.N * sizeof(int));
	for (int a = 2; a < argc; a++){
		int n = ParseSites(argv[a], sites, h.N);
		double S = EntanglementEntropy(argv[1], sites, n);
		if (Rank() == 0)
			printf("A = %s : S = %.15f\n", argv[a], S);
	}

	free(sites);
	FreeTables(h.N);
	ParallelFinalize();
	return 0;
}
