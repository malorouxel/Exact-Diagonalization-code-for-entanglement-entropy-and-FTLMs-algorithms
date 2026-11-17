#include "thermodynamics.h"
#include "tables.h"
#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Return the number of Sz sectors to compute.
// With the Sz symmetry we only keep the sectors with Sz <= 0. The others follow by spin flip.
int NumberOfSectors(const Parameters *P){
	double spin = 0.5 + 0.5 * P->SpinOne;
	return (int)(P->N * spin * P->SzSym) + 1;
}

// Read a spectral file. Each line holds an energy and a weight.
// Return the number of lines or -1 if the file does not exist.
static int ReadSpectralFile(const char *path, double **energies, double **weights){
	FILE *f = fopen(path, "r");
	if (!f)
		return -1;
	int capacity = 256;
	int n = 0;
	double *e = malloc(capacity * sizeof(double));
	double *w = malloc(capacity * sizeof(double));
	double ev, wv;
	while (fscanf(f, "%lf %lf", &ev, &wv) == 2){
		if (n == capacity){
			capacity *= 2;
			e = realloc(e, capacity * sizeof(double));
			w = realloc(w, capacity * sizeof(double));
		}
		e[n] = ev;
		w[n] = wv;
		n++;
	}
	fclose(f);
	*energies = e;
	*weights = w;
	return n;
}

// Add the Boltzmann sums of a list of levels to Z, E and E2.
// factor multiplies every weight. Energies are shifted by E0min to avoid overflow.
static void AddLevels(double *Z, double *E, double *E2, const double *beta, int NT, const double *energies, const double *weights, int start, int n, double factor, double E0min){
	for (int t = 0; t < NT; t++){
		for (int j = start; j < n; j++){
			double boltz = factor * weights[j] * exp(-beta[t] * (energies[j] - E0min));
			Z[t]  += boltz;
			E[t]  += boltz * energies[j];
			E2[t] += boltz * energies[j] * energies[j];
		}
	}
}

// Compute the specific heat from the spectral files of a finished run.
// Sectors solved by exact diagonalization enter exactly.
// Sectors solved by Lanczos enter through the R random vectors.
// SFTLM is the standard FTLM. RFTLM treats the lowest level exactly with its known degeneracy.
// OFTLM uses Lanczos runs orthogonal to the ground states and adds the ground states exactly.
// The error bars come from a jackknife over the random vectors.
// Write one file folder/<method>_R<R>_thermo.txt per method.
void ComputeThermodynamics(const Parameters *P){
	int N = P->N;
	int R = P->R;
	int NT = P->NT;
	int nsectors = NumberOfSectors(P);
	double spin = 0.5 + 0.5 * P->SpinOne;
	char path[512];

	double *T = malloc(NT * sizeof(double));
	double *beta = malloc(NT * sizeof(double));
	for (int t = 0; t < NT; t++){
		double exponent = (NT == 1) ? P->logTmin : P->logTmin + (P->logTmax - P->logTmin) * t / (NT - 1);
		T[t] = pow(10.0, exponent);
		beta[t] = 1.0 / T[t];
	}

	double *E0sector = calloc(nsectors, sizeof(double));
	snprintf(path, sizeof(path), "%s/Results.txt", P->FolderName);
	FILE *fres = fopen(path, "r");
	if (fres){
		char line[1024];
		while (fgets(line, sizeof(line), fres)){
			int sz;
			double e0;
			if (sscanf(line, "E0 Sz = %d : %lf", &sz, &e0) == 2 && sz >= 0 && sz < nsectors)
				E0sector[sz] = e0;
		}
		fclose(fres);
	}
	double E0min = E0sector[0];
	for (int s = 1; s < nsectors; s++)
		if (E0sector[s] < E0min)
			E0min = E0sector[s];

	int *DegeSector = malloc(nsectors * sizeof(int));
	for (int s = 0; s < nsectors; s++)
		DegeSector[s] = 1;
	int DegeFileAvailable = 0;
	snprintf(path, sizeof(path), "%s/dege.txt", P->FolderName);
	FILE *fdeg = fopen(path, "r");
	if (fdeg){
		char line[1024];
		while (fgets(line, sizeof(line), fdeg)){
			int sz, dg;
			if (sscanf(line, "Sz = %d : Dege = %d", &sz, &dg) == 2 && sz >= 0 && sz < nsectors){
				DegeSector[sz] = dg;
				DegeFileAvailable = 1;
			}
		}
		fclose(fdeg);
	}

	double nspin = N * spin;
	int CenterIsSector = (fabs(nspin - round(nspin)) < 1e-10);

	int methods[3] = { P->ResultsFTLM, P->ResultsRFTLM, P->ResultsORFTLM };
	const char *prefixes[3] = { "", "", "O" };
	char names[3][32];
	snprintf(names[0], 32, "SFTLM_R%d", R);
	snprintf(names[1], 32, "RFTLM_R%d", R);
	snprintf(names[2], 32, "OFTLM_R%d", R);

	for (int m = 0; m < 3; m++){
		if (!methods[m])
			continue;
		if (m == 1 && !DegeFileAvailable && !P->KnownDege){
			fprintf(stderr, "RFTLM skipped. The degeneracy is unknown. Use dege in Results or give it in the parameters.\n");
			continue;
		}
		const char *prefix = prefixes[m];

		double *Zexact  = calloc(NT, sizeof(double));
		double *Eexact  = calloc(NT, sizeof(double));
		double *E2exact = calloc(NT, sizeof(double));
		double **Zr  = malloc(R * sizeof(double *));
		double **Er  = malloc(R * sizeof(double *));
		double **E2r = malloc(R * sizeof(double *));
		for (int r = 0; r < R; r++){
			Zr[r]  = calloc(NT, sizeof(double));
			Er[r]  = calloc(NT, sizeof(double));
			E2r[r] = calloc(NT, sizeof(double));
		}

		for (int Sz = 0; Sz < nsectors; Sz++){
			unsigned long long Size = SectorSize(Sz, N, P->SpinOne, P->SzSym);
			double multiplicity = 2.0;
			if (P->SzSym == 0 || (CenterIsSector && Sz == nsectors - 1))
				multiplicity = 1.0;

			double *energies, *weights;
			snprintf(path, sizeof(path), "%s/%s%d_exact.txt", P->FolderName, prefix, Sz);
			int nlines = ReadSpectralFile(path, &energies, &weights);
			if (nlines > 0){
				AddLevels(Zexact, Eexact, E2exact, beta, NT, energies, weights, 0, nlines, multiplicity * Size, E0min);
				free(energies);
				free(weights);
				continue;
			}

			if (m == 2){
				double one = 1.0;
				AddLevels(Zexact, Eexact, E2exact, beta, NT, &E0sector[Sz], &one, 0, 1, multiplicity * DegeSector[Sz], E0min);
			}
			double factor = (m == 2) ? (double)(Size - DegeSector[Sz]) : (double)Size;

			for (int r = 0; r < R; r++){
				snprintf(path, sizeof(path), "%s/%s%d_%d.txt", P->FolderName, prefix, Sz, r);
				nlines = ReadSpectralFile(path, &energies, &weights);
				if (nlines <= 0)
					continue;
				int start = 0;
				if (m == 1){
					double E0file = energies[0];
					while (start < nlines && fabs(energies[start] - E0file) < 1e-10)
						start++;
					int dege = DegeFileAvailable ? DegeSector[Sz] : P->Dege;
					double one = 1.0;
					AddLevels(Zr[r], Er[r], E2r[r], beta, NT, &E0file, &one, 0, 1, multiplicity * dege, E0min);
				}
				AddLevels(Zr[r], Er[r], E2r[r], beta, NT, energies, weights, start, nlines, multiplicity * factor, E0min);
				free(energies);
				free(weights);
			}
		}

		double *Zavg  = calloc(NT, sizeof(double));
		double *Eavg  = calloc(NT, sizeof(double));
		double *E2avg = calloc(NT, sizeof(double));
		for (int r = 0; r < R; r++){
			for (int t = 0; t < NT; t++){
				Zavg[t]  += Zr[r][t] / R;
				Eavg[t]  += Er[r][t] / R;
				E2avg[t] += E2r[r][t] / R;
			}
		}

		double *Zsum  = malloc(NT * sizeof(double));
		double *Emean = malloc(NT * sizeof(double));
		double *Cfull = malloc(NT * sizeof(double));
		for (int t = 0; t < NT; t++){
			Zsum[t] = Zexact[t] + Zavg[t];
			Emean[t] = (Eexact[t] + Eavg[t]) / Zsum[t];
			double E2mean = (E2exact[t] + E2avg[t]) / Zsum[t];
			Cfull[t] = (E2mean - Emean[t] * Emean[t]) / (T[t] * T[t]);
		}

		double *CJK = calloc(NT, sizeof(double));
		double *sigma = calloc(NT, sizeof(double));
		if (R >= 2){
			double *Cbar = calloc(NT, sizeof(double));
			double **Creps = malloc(R * sizeof(double *));
			for (int k = 0; k < R; k++){
				Creps[k] = malloc(NT * sizeof(double));
				for (int t = 0; t < NT; t++){
					double Zk  = Zexact[t]  + (R * Zavg[t]  - Zr[k][t])  / (R - 1);
					double Ek  = Eexact[t]  + (R * Eavg[t]  - Er[k][t])  / (R - 1);
					double E2k = E2exact[t] + (R * E2avg[t] - E2r[k][t]) / (R - 1);
					double em = Ek / Zk;
					Creps[k][t] = (E2k / Zk - em * em) / (T[t] * T[t]);
					Cbar[t] += Creps[k][t] / R;
				}
			}
			for (int t = 0; t < NT; t++){
				CJK[t] = R * Cfull[t] - (R - 1) * Cbar[t];
				double var = 0.0;
				for (int k = 0; k < R; k++)
					var += (Creps[k][t] - Cbar[t]) * (Creps[k][t] - Cbar[t]);
				sigma[t] = sqrt((double)(R - 1) / R * var);
			}
			for (int k = 0; k < R; k++)
				free(Creps[k]);
			free(Creps);
			free(Cbar);
		}
		else{
			for (int t = 0; t < NT; t++)
				CJK[t] = Cfull[t];
		}

		snprintf(path, sizeof(path), "%s/%s_thermo.txt", P->FolderName, names[m]);
		FILE *fout = fopen(path, "w");
		fprintf(fout, "# %s specific heat - Jackknife error bars\n", names[m]);
		fprintf(fout, "# N = %d  R = %d  E0_min = %.15f\n", N, R, E0min);
		fprintf(fout, "# T  Z  E/N  C_full/N  C_JK/N  C_low1  C_high1  C_low2  C_high2\n");
		for (int t = 0; t < NT; t++){
			long double Zreal = (long double)Zsum[t] * expl(-(long double)beta[t] * (long double)E0min);
			double cf = Cfull[t] / N;
			double cj = CJK[t] / N;
			double sn = sigma[t] / N;
			fprintf(fout, "%.8e  %.8Le  %.8e  %.8e  %.8e  %.8e  %.8e  %.8e  %.8e\n", T[t], Zreal, Emean[t] / N, cf, cj, cj - sn, cj + sn, cj - 2.0 * sn, cj + 2.0 * sn);
		}
		fclose(fout);

		free(Zexact); free(Eexact); free(E2exact);
		free(Zavg); free(Eavg); free(E2avg);
		free(Zsum); free(Emean); free(Cfull);
		free(CJK); free(sigma);
		for (int r = 0; r < R; r++){
			free(Zr[r]); free(Er[r]); free(E2r[r]);
		}
		free(Zr); free(Er); free(E2r);
	}

	free(T);
	free(beta);
	free(E0sector);
	free(DegeSector);
}
