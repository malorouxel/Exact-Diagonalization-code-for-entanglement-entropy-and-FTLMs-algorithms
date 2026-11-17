#include "io.h"
#include "parallel.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Skip the label line and read the value line that follows it.
// Return 0 at the end of the file.
static int ReadValueLine(FILE *file, char *line, int length){
	if (!fgets(line, length, file))
		return 0;
	if (!fgets(line, length, file))
		return 0;
	return 1;
}

// Read the run parameters.
// The file alternates a label line and a value line. The order of the entries is fixed.
// The last entry Seed is optional. Without it the seed comes from the clock.
void ReadParameters(const char *path, Parameters *P){
	char line[1024];
	FILE *file = fopen(path, "r");
	if (!file)
		Abort("Cannot open the parameters file.");

	ReadValueLine(file, line, sizeof(line));
	P->SpinOne = (strncmp(line, "One", 3) == 0) ? 1 : 0;

	ReadValueLine(file, line, sizeof(line));
	sscanf(line, "%d", &P->N);

	ReadValueLine(file, line, sizeof(line));
	sscanf(line, "%llu", &P->EDThreshold);

	ReadValueLine(file, line, sizeof(line));
	P->SzSym = (strncmp(line, "Yes", 3) == 0) ? 1 : 0;

	ReadValueLine(file, line, sizeof(line));
	sscanf(line, "%d", &P->R);

	ReadValueLine(file, line, sizeof(line));
	if (strncmp(line, "Conv", 4) == 0){
		P->LanczosConv = 1;
		P->NIterationMax = 1000;
	}
	else{
		P->LanczosConv = 0;
		sscanf(line, "%d", &P->NIterationMax);
	}

	ReadValueLine(file, line, sizeof(line));
	if (strncmp(line, "Unknown", 7) == 0){
		P->KnownDege = 0;
		P->Dege = 1;
	}
	else{
		P->KnownDege = 1;
		sscanf(line, "%d", &P->Dege);
	}

	ReadValueLine(file, line, sizeof(line));
	P->ResultsE0      = (strstr(line, "E0")      != NULL);
	P->ResultsGS      = (strstr(line, "GS")      != NULL);
	P->ResultsDege    = (strstr(line, "dege")    != NULL);
	P->ResultsEntropy = (strstr(line, "Entropy") != NULL);
	P->ResultsORFTLM  = (strstr(line, "OFTLM")   != NULL);
	P->ResultsRFTLM   = (strstr(line, "RFTLM")   != NULL);
	P->ResultsFTLM    = (strstr(line, "SFTLM")   != NULL);

	ReadValueLine(file, line, sizeof(line));
	sscanf(line, "%d", &P->BipartitionSize);

	ReadValueLine(file, line, sizeof(line));
	char folder[200];
	sscanf(line, "%199s", folder);
	snprintf(P->FolderName, sizeof(P->FolderName), "%s_N%d", folder, P->N);

	ReadValueLine(file, line, sizeof(line));
	sscanf(line, "%lf %lf %d", &P->logTmin, &P->logTmax, &P->NT);

	P->RandomSeed = 1;
	P->Seed = 0;
	if (ReadValueLine(file, line, sizeof(line)) && strncmp(line, "Random", 6) != 0){
		if (sscanf(line, "%llu", &P->Seed) == 1)
			P->RandomSeed = 0;
	}

	fclose(file);
}

// Read the list of interactions.
// The first line is a header. Each other line has the form J;OP;SITES.
// OP is SS for S_i.S_j, SS2 for (S_i.S_j)^2 and Sz for a field on one site.
// SITES is i-j for SS and SS2 and i for Sz. Sites start at 0.
// SS is split into Sz Sz plus one half of S+ S- and S- S+.
void ReadInteractions(const char *path, Interactions *H){
	char line[1024];
	FILE *f = fopen(path, "r");
	if (!f)
		Abort("Cannot open the interactions file.");
	if (!fgets(line, sizeof(line), f))
		Abort("The interactions file is empty.");
	while (fgets(line, sizeof(line), f)){
		double J;
		char op[16];
		char sites[32];
		if (sscanf(line, "%lf;%15[^;];%31s", &J, op, sites) != 3)
			continue;
		if (strcmp(op, "SS") == 0){
			int si, sj;
			sscanf(sites, "%d-%d", &si, &sj);
			PushTerm(H, J, OP_SzSz, si, sj);
			PushTerm(H, J / 2.0, OP_SpSm, si, sj);
			PushTerm(H, J / 2.0, OP_SmSp, si, sj);
		}
		else if (strcmp(op, "SS2") == 0){
			int si, sj;
			sscanf(sites, "%d-%d", &si, &sj);
			PushTerm(H, J, OP_SS2, si, sj);
		}
		else if (strcmp(op, "Sz") == 0){
			int si;
			sscanf(sites, "%d", &si);
			PushTerm(H, J, OP_Sz, si, -1);
		}
	}
	fclose(f);
}

// Write a ground state to a binary file.
// The file starts with a header of 32 bytes. The Size amplitudes follow in basis order.
// Every process writes its own slice. Every process calls this function.
void WriteGroundState(const char *path, const GroundStateHeader *header, const double *values, unsigned long long start, unsigned long long count){
	if (Rank() == 0){
		FILE *f = fopen(path, "wb");
		if (!f)
			Abort("Cannot create the ground state file.");
		fwrite(&header->N, sizeof(int), 1, f);
		fwrite(&header->SpinOne, sizeof(int), 1, f);
		fwrite(&header->SzSym, sizeof(int), 1, f);
		fwrite(&header->Sz, sizeof(int), 1, f);
		fwrite(&header->Size, sizeof(unsigned long long), 1, f);
		fwrite(&header->E0, sizeof(double), 1, f);
		fclose(f);
	}
	Barrier();
	WriteVectorSlice(path, HEADER_BYTES + (long long)(start * sizeof(double)), values, count);
	Barrier();
}

// Read the header of a ground state file.
void ReadGroundStateHeader(const char *path, GroundStateHeader *header){
	FILE *f = fopen(path, "rb");
	if (!f)
		Abort("Cannot open the ground state file.");
	int ok = 1;
	ok &= fread(&header->N, sizeof(int), 1, f) == 1;
	ok &= fread(&header->SpinOne, sizeof(int), 1, f) == 1;
	ok &= fread(&header->SzSym, sizeof(int), 1, f) == 1;
	ok &= fread(&header->Sz, sizeof(int), 1, f) == 1;
	ok &= fread(&header->Size, sizeof(unsigned long long), 1, f) == 1;
	ok &= fread(&header->E0, sizeof(double), 1, f) == 1;
	fclose(f);
	if (!ok)
		Abort("The ground state file is too short.");
}
