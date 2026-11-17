#ifndef TABLES_H
#define TABLES_H

extern unsigned long long **Binomial;
extern unsigned long long *Pow2;
extern unsigned long long *Pow3;

void CreateTables(int N);
void FreeTables(int N);
unsigned long long SectorSize(int Sz, int N, int SpinOne, int SzSym);

#endif
