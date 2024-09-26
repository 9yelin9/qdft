#ifndef HUB1D_H
#define HUB1D_H

#define USE_MATH_DEFINES

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <sys/stat.h>

#define double_complex double _Complex

#define SYSTEM "hub1d"
#define BETA 1.9999592596

typedef struct {
	int Ni; // # of sites
	int Ne; // # of electrons
	int N; // # of spin-orbitals (2 * Ni)
	int Nb; // # of bases
	double U; // on-site Coulomb interaction
	int verbose;
} params;

typedef struct {
	int idx; // index in bases
	int val; // | val (10) > = | 1up 2up ... Niup 1dn 2dn ... Nidn (2) >
} basis;

typedef struct {
	int nnz; // # of nonzero elements in H
	int *row; // state (10)
	int *col; // state (10)
	int *col_csc; // beginning of each column in H (csc: compressed spares columns)
	double_complex *val; // < row | H | col >
} hamiltonian;

void dec2bin(int len, int dec, char *bin);
void gen_dir_output(params *pm, char *dir_output);
void print_basis(params *pm, basis *b, char *dir_output, char *method);
void print_H(params *pm, hamiltonian *H, char *dir_output, char *method);
double laeig(params *pm, hamiltonian *H); // LAPACK_zheev eigensolver (full matrix)

#ifdef __cplusplus
extern "C" {
#endif
	double areig(params *pm, hamiltonian *H); // ARPACK AREig eigensolver (csc matrix)
#ifdef __cplusplus
}
#endif

#endif
