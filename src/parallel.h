#ifndef PARALLEL_H
#define PARALLEL_H

void ParallelInit(int *argc, char ***argv);
void ParallelFinalize(void);
int Rank(void);
int NumberOfRanks(void);
void Abort(const char *message);
void Barrier(void);

double SumDouble(double value);
void SumDoubleArray(double *values, int count);
void BroadcastDouble(double *value);

unsigned long long SliceStart(unsigned long long Size, int rank, int nranks);
int SliceOwner(unsigned long long index, unsigned long long Size, int nranks);

void ExchangeCounts(const int *sendcounts, int *recvcounts);
void ExchangeIndices(const unsigned long long *sendbuf, const int *sendcounts, const int *senddispls, unsigned long long *recvbuf, const int *recvcounts, const int *recvdispls);
void ExchangeValues(const double *sendbuf, const int *sendcounts, const int *senddispls, double *recvbuf, const int *recvcounts, const int *recvdispls);

void WriteVectorSlice(const char *path, long long offset, const double *values, unsigned long long count);

#endif
