#ifndef QDFT_H
#define QDFT_H

#define USE_MATH_DEFINES

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#define double_complex double _Complex

#define V(i) (i / 10.)

typedef struct {
	int Ne; // # of electrons
	int N; // # of spin-orbitals (2*Ne)
	int M; // # of states (2^N)
} environment;

typedef struct {
	int nnz; // # of nonzero elements in H
	int *row; // computational basis |I> = | 1up 2up ... Neup 1dn 2dn ... Nedn >
	int *col; // computational basis |J>
	int *col_csr; // beginning of each column in H
	double_complex *val; // H_IJ = <I|H|J>
} hamiltonian;

#endif
