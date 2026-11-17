#ifndef COMMON_H
#define COMMON_H

#ifdef USE_MKL
#include <mkl.h>
#else
#include <cblas.h>
#include <lapacke.h>
#endif

#define PI                  3.14159265358979323846
#define EPSILON             1e-12
#define CONV_THRESHOLD      1e-6
#define MAX_ELEMENTS        8

#endif
