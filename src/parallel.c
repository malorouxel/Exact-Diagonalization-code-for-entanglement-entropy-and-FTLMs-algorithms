#include "parallel.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

#define IO_CHUNK (1ULL << 26)

// Start MPI. Do nothing in the serial build.
void ParallelInit(int *argc, char ***argv){
#ifdef USE_MPI
	MPI_Init(argc, argv);
#else
	(void)argc;
	(void)argv;
#endif
}

// Stop MPI. Do nothing in the serial build.
void ParallelFinalize(void){
#ifdef USE_MPI
	MPI_Finalize();
#endif
}

// Return the rank of this process. Return 0 in the serial build.
int Rank(void){
#ifdef USE_MPI
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	return rank;
#else
	return 0;
#endif
}

// Return the number of processes. Return 1 in the serial build.
int NumberOfRanks(void){
#ifdef USE_MPI
	int nranks;
	MPI_Comm_size(MPI_COMM_WORLD, &nranks);
	return nranks;
#else
	return 1;
#endif
}

// Print an error message and stop every process.
void Abort(const char *message){
	fprintf(stderr, "[rank %d] Error. %s\n", Rank(), message);
#ifdef USE_MPI
	MPI_Abort(MPI_COMM_WORLD, 1);
#endif
	exit(1);
}

// Wait for every process.
void Barrier(void){
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
}

// Sum a number over all processes.
double SumDouble(double value){
#ifdef USE_MPI
	double total;
	MPI_Allreduce(&value, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return total;
#else
	return value;
#endif
}

// Sum an array over all processes. The result replaces the array.
void SumDoubleArray(double *values, int count){
#ifdef USE_MPI
	MPI_Allreduce(MPI_IN_PLACE, values, count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
	(void)values;
	(void)count;
#endif
}

// Copy a number from rank 0 to every process.
void BroadcastDouble(double *value){
#ifdef USE_MPI
	MPI_Bcast(value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
	(void)value;
#endif
}

// Return the first index owned by a rank.
// Each rank owns a contiguous slice of the basis.
// The first Size % nranks ranks get one more element.
unsigned long long SliceStart(unsigned long long Size, int rank, int nranks){
	unsigned long long base  = Size / nranks;
	unsigned long long extra = Size % nranks;
	unsigned long long r = rank;
	return r * base + (r < extra ? r : extra);
}

// Return the rank that owns a basis index.
int SliceOwner(unsigned long long index, unsigned long long Size, int nranks){
	unsigned long long base  = Size / nranks;
	unsigned long long extra = Size % nranks;
	unsigned long long limit = extra * (base + 1);
	if (index < limit)
		return (int)(index / (base + 1));
	return (int)(extra + (index - limit) / base);
}

// Send one count to each process and receive one count from each process.
void ExchangeCounts(const int *sendcounts, int *recvcounts){
#ifdef USE_MPI
	MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);
#else
	recvcounts[0] = sendcounts[0];
#endif
}

// Send blocks of basis indices to each process.
void ExchangeIndices(const unsigned long long *sendbuf, const int *sendcounts, const int *senddispls, unsigned long long *recvbuf, const int *recvcounts, const int *recvdispls){
#ifdef USE_MPI
	MPI_Alltoallv(sendbuf, sendcounts, senddispls, MPI_UNSIGNED_LONG_LONG, recvbuf, recvcounts, recvdispls, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
#else
	memcpy(recvbuf + recvdispls[0], sendbuf + senddispls[0], sendcounts[0] * sizeof(unsigned long long));
	(void)recvcounts;
#endif
}

// Send blocks of vector values to each process.
void ExchangeValues(const double *sendbuf, const int *sendcounts, const int *senddispls, double *recvbuf, const int *recvcounts, const int *recvdispls){
#ifdef USE_MPI
	MPI_Alltoallv(sendbuf, sendcounts, senddispls, MPI_DOUBLE, recvbuf, recvcounts, recvdispls, MPI_DOUBLE, MPI_COMM_WORLD);
#else
	memcpy(recvbuf + recvdispls[0], sendbuf + senddispls[0], sendcounts[0] * sizeof(double));
	(void)recvcounts;
#endif
}

// Write the local slice of a vector into an existing file.
// The slice starts at a byte offset. Every process calls this function.
void WriteVectorSlice(const char *path, long long offset, const double *values, unsigned long long count){
#ifdef USE_MPI
	MPI_File fh;
	if (MPI_File_open(MPI_COMM_WORLD, path, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh) != MPI_SUCCESS)
		Abort("Cannot open the output file.");
	unsigned long long done = 0;
	while (done < count){
		unsigned long long n = count - done;
		if (n > IO_CHUNK) n = IO_CHUNK;
		MPI_File_write_at(fh, (MPI_Offset)(offset + done * sizeof(double)), values + done, (int)n, MPI_DOUBLE, MPI_STATUS_IGNORE);
		done += n;
	}
	MPI_File_close(&fh);
#else
	FILE *f = fopen(path, "r+b");
	if (!f)
		Abort("Cannot open the output file.");
	fseek(f, offset, SEEK_SET);
	fwrite(values, sizeof(double), count, f);
	fclose(f);
#endif
}
