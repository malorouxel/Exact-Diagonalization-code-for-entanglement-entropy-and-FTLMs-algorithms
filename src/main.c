#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#include "common.h"
#include "parallel.h"
#include "tables.h"
#include "sector.h"
#include "hamiltonian.h"
#include "lanczos.h"
#include "observables.h"
#include "entropy.h"
#include "thermodynamics.h"
#include "io.h"
#include "logger.h"

// Solve one small sector by exact diagonalization on rank 0.
// Write the spectrum files and the degeneracy when they are requested.
// Return the lowest energy of the sector on every process.
static double DenseSector(Parameters *P, Interactions *H, int Sz, FILE *fd, FILE *flog){
	double Esector = 0.0;
	if (Rank() == 0){
		Sector s;
		SectorBuild(&s, P->N, Sz, P->SzSym, P->SpinOne, H, 0);
		LogRAM(flog, Sz, EstimationRAM(s.size, 1, 1));
		double *EigVals = malloc(s.size * sizeof(double));
		double *Matrix  = malloc(s.size * s.size * sizeof(double));
		if (!EigVals || !Matrix)
			Abort("Not enough memory for exact diagonalization. Lower the ED threshold.");

		double t0 = WallTime();
		ExactDiagonalization(&s, EigVals, Matrix, 'N');
		LogStep(flog, Sz, "ExactDiagonalization", WallTime() - t0);
		Esector = EigVals[0];

		if (P->ResultsFTLM || P->ResultsRFTLM)
			WriteExactSpectrum(EigVals, s.size, Sz, P->FolderName, "");
		if (P->ResultsDege || P->ResultsORFTLM){
			int degen = 0;
			for (unsigned long long i = 0; i < s.size; i++)
				if (fabs(EigVals[i] - EigVals[0]) < EPSILON)
					degen++;
			fprintf(fd, "Sz = %d : Dege = %d\n", Sz, degen);
			fflush(fd);
		}
		if (P->ResultsORFTLM)
			WriteExactSpectrum(EigVals, s.size, Sz, P->FolderName, "O");

		free(EigVals);
		free(Matrix);
		SectorFree(&s);
	}
	BroadcastDouble(&Esector);
	return Esector;
}

// Solve one large sector with Lanczos. The basis is split between the processes.
// Compute the ground state energy and the requested FTLM data.
// run receives the coefficients of the ground state run.
static double LanczosSector(Parameters *P, Interactions *H, int Sz, unsigned long long *seed, LanczosRun *run, FILE *fd, FILE *flog){
	int root = (Rank() == 0);
	Sector s;
	double t0 = WallTime();
	SectorBuild(&s, P->N, Sz, P->SzSym, P->SpinOne, H, 1);
	if (root){
		LogStep(flog, Sz, "SectorBuild", WallTime() - t0);
		int nvec = (P->ResultsDege || P->ResultsORFTLM) ? 6 : 3;
		LogRAM(flog, Sz, EstimationRAM(s.nlocal, nvec, 0) + s.nghost * 16.0);
	}

	t0 = WallTime();
	double Esector = FirstLanczos(&s, seed, run, P->LanczosConv, P->NIterationMax);
	if (root)
		LogStep(flog, Sz, "FirstLanczos", WallTime() - t0);

	LanczosRun tmp;
	LanczosRunInit(&tmp, P->NIterationMax);

	double **OrthoVecs = NULL;
	int degen = 0;
	if (P->ResultsDege || P->ResultsORFTLM){
		OrthoVecs = malloc(sizeof(double *));
		OrthoVecs[0] = malloc((s.nlocal + 1) * sizeof(double));
		t0 = WallTime();
		SecondLanczos(&s, run, OrthoVecs[0]);
		degen = FindDegenerateStates(&s, Esector, &OrthoVecs, seed, P->LanczosConv, P->NIterationMax);
		if (root){
			LogStep(flog, Sz, "Degeneracy", WallTime() - t0);
			fprintf(fd, "Sz = %d : Dege = %d\n", Sz, degen);
			fflush(fd);
		}
	}

	if (P->ResultsFTLM || P->ResultsRFTLM){
		t0 = WallTime();
		for (int r = 0; r < P->R; r++){
			FirstLanczos(&s, seed, &tmp, P->LanczosConv, P->NIterationMax);
			if (root)
				WriteSpectralWeights(&tmp, Sz, r, P->FolderName, "");
		}
		if (root)
			LogStep(flog, Sz, "FTLM", WallTime() - t0);
	}

	if (P->ResultsORFTLM && (unsigned long long)degen < s.size){
		t0 = WallTime();
		for (int r = 0; r < P->R; r++){
			OrthogonalLanczos(&s, OrthoVecs, degen, seed, &tmp, P->LanczosConv, P->NIterationMax);
			if (root)
				WriteSpectralWeights(&tmp, Sz, r, P->FolderName, "O");
		}
		if (root)
			LogStep(flog, Sz, "OFTLM", WallTime() - t0);
	}

	for (int d = 0; d < degen; d++)
		free(OrthoVecs[d]);
	free(OrthoVecs);
	LanczosRunFree(&tmp);
	SectorFree(&s);
	return Esector;
}

// Compute the ground state of the best sector again and write it to folder/GroundState.bin.
// Use exact diagonalization on rank 0 for a small sector. Use Lanczos for a large sector.
static void WriteBestGroundState(Parameters *P, Interactions *H, int bestSz, int bestIsED, LanczosRun *best, double E0, const char *path, FILE *flog){
	GroundStateHeader header = { P->N, P->SpinOne, P->SzSym, bestSz, SectorSize(bestSz, P->N, P->SpinOne, P->SzSym), E0 };
	double t0 = WallTime();
	if (bestIsED){
		double *Matrix = NULL;
		if (Rank() == 0){
			Sector s;
			SectorBuild(&s, P->N, bestSz, P->SzSym, P->SpinOne, H, 0);
			double *EigVals = malloc(s.size * sizeof(double));
			Matrix = malloc(s.size * s.size * sizeof(double));
			ExactDiagonalization(&s, EigVals, Matrix, 'V');
			free(EigVals);
			SectorFree(&s);
		}
		WriteGroundState(path, &header, Matrix, 0, (Rank() == 0) ? header.Size : 0);
		free(Matrix);
	}
	else{
		Sector s;
		SectorBuild(&s, P->N, bestSz, P->SzSym, P->SpinOne, H, 1);
		double *GroundState = malloc((s.nlocal + 1) * sizeof(double));
		SecondLanczos(&s, best, GroundState);
		WriteGroundState(path, &header, GroundState, s.start, s.nlocal);
		free(GroundState);
		SectorFree(&s);
	}
	if (Rank() == 0)
		LogStep(flog, bestSz, "GroundState", WallTime() - t0);
}

// Run the full calculation described by a parameters file and an interactions file.
// Usage. lanczos [parameters file] [interactions file]
int main(int argc, char **argv){
	ParallelInit(&argc, &argv);
	int root = (Rank() == 0);

	const char *params_path = (argc > 1) ? argv[1] : "InputsParameters.txt";
	const char *inter_path  = (argc > 2) ? argv[2] : "Interactions.txt";

	Parameters P;
	ReadParameters(params_path, &P);
	CreateTables(P.N);
	Interactions H = { NULL, 0 };
	ReadInteractions(inter_path, &H);

	double clock_seed = (double)time(NULL);
	BroadcastDouble(&clock_seed);
	unsigned long long seed = P.RandomSeed ? (unsigned long long)clock_seed : P.Seed;

	FILE *flog = NULL;
	FILE *fr = NULL;
	FILE *fd = NULL;
	char path[512];
	if (root){
		mkdir(P.FolderName, 0755);
		snprintf(path, sizeof(path), "%s/Run.log", P.FolderName);
		flog = fopen(path, "w");
		fprintf(flog, "Processes = %d : Seed = %llu\n", NumberOfRanks(), seed);
		snprintf(path, sizeof(path), "%s/Results.txt", P.FolderName);
		fr = fopen(path, "w");
		if (P.ResultsDege || P.ResultsORFTLM){
			snprintf(path, sizeof(path), "%s/dege.txt", P.FolderName);
			fd = fopen(path, "w");
		}
	}
	Barrier();

	LanczosRun run, best;
	LanczosRunInit(&run, P.NIterationMax);
	LanczosRunInit(&best, P.NIterationMax);

	int nsectors = NumberOfSectors(&P);
	double center = P.N * (0.5 + 0.5 * P.SpinOne);
	double E0global = 0.0;
	int bestSz = 0;
	int bestIsED = 0;

	for (int Sz = 0; Sz < nsectors; Sz++){
		unsigned long long Size = SectorSize(Sz, P.N, P.SpinOne, P.SzSym);
		int dense = (Size < P.EDThreshold);
		double Esector;
		if (dense)
			Esector = DenseSector(&P, &H, Sz, fd, flog);
		else
			Esector = LanczosSector(&P, &H, Sz, &seed, &run, fd, flog);

		if (root){
			fprintf(fr, "E0 Sz = %d : %.15f\n", Sz, Esector);
			fflush(fr);
		}

		int newbest = (Sz == 0) || (Esector < E0global - EPSILON);
		if (!newbest && fabs(Esector - E0global) < EPSILON && fabs(Sz - center) < fabs(bestSz - center))
			newbest = 1;
		if (newbest){
			E0global = Esector;
			bestSz = Sz;
			bestIsED = dense;
			if (!dense)
				LanczosRunCopy(&best, &run);
		}
	}

	if (root){
		fprintf(fr, "Ground state : Sz = %d : E0 = %.15f\n", bestSz, E0global);
		fflush(fr);
	}

	if (P.ResultsGS || P.ResultsEntropy){
		snprintf(path, sizeof(path), "%s/GroundState.bin", P.FolderName);
		WriteBestGroundState(&P, &H, bestSz, bestIsED, &best, E0global, path, flog);
		if (P.ResultsEntropy){
			int *sites = malloc((P.BipartitionSize + 1) * sizeof(int));
			for (int k = 0; k < P.BipartitionSize; k++)
				sites[k] = k;
			double t0 = WallTime();
			double S = EntanglementEntropy(path, sites, P.BipartitionSize);
			if (root){
				LogStep(flog, bestSz, "Entropy", WallTime() - t0);
				fprintf(fr, "Entropy = %.15f\n", S);
			}
			free(sites);
		}
	}

	if (root){
		if (fd) fclose(fd);
		fclose(fr);
		fclose(flog);
		if (P.ResultsFTLM || P.ResultsRFTLM || P.ResultsORFTLM)
			ComputeThermodynamics(&P);
	}

	LanczosRunFree(&run);
	LanczosRunFree(&best);
	FreeInteractions(&H);
	FreeTables(P.N);
	ParallelFinalize();
	return 0;
}
