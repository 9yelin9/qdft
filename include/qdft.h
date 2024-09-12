#ifndef QDFT_H
#define QDFT_H

#define USE_MATH_DEFINES

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#define double_complex double _Complex

typedef struct {
	int Ni; // # of sites
	int Ne; // # of electrons
	int N; // # of spin-orbitals (2*Ni)
	int M; // # of states
	double U; // on-site Coulomb interaction
} environment;

typedef struct {
	int nnz; // # of nonzero elements in H
	int *row; // computational basis |I> = | 1up 2up ... Neup 1dn 2dn ... Nedn >
	int *col; // computational basis |J>
	int *col_csc; // beginning of each column in H
	double_complex *val; // H_IJ = <I|H|J>
} hamiltonian;

void print_H(environment *env, hamiltonian *H);
double laeig(environment *env, hamiltonian *H); // LAPACK_zheev eigensolver (full matrix)

#ifdef __cplusplus
extern "C" {
#endif
	double areig(environment *env, hamiltonian *H); // ARPACK AREig eigensolver (csc matrix)
#ifdef __cplusplus
}
#endif

#endif
