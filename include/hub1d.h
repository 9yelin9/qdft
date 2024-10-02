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
#define square_complex(c) (creal(c)*creal(c) + cimag(c)*cimag(c))

#define SYSTEM "hub1d"

typedef struct {
	int Ni; // # of sites
	int Ne; // # of electrons
	int N; // # of spin-orbitals (2 * Ni)
	int Nb; // # of bases
	double U; // on-site Coulomb interaction
	double beta; // bethe ansatz parameter
} params;

typedef struct {
	int idx; // index in bases
	int val; // | state (10) > = | 1up 2up ... Niup 1dn 2dn ... Nidn (2) >
} basis;

typedef struct {
	int nnz; // # of nonzero elements in H
	int *row; // state (10)
	int *col; // state (10)
	int *col_csc; // beginning of each column in H (csc: compressed spares columns)
	double e_grd; // ground-state energy
	double_complex *val; // < row | H | col >
} hamiltonian;

int compare_basis(const void *key, const void *p);
int check_hamiltonian_hermitian(params *pm, hamiltonian *h);
int check_coo2csc(params *pm, hamiltonian *h);
void dec2bin(int len, int dec, char *bin);
void gen_dir_output(params *pm, char *dir_output);
void print_basis(params *pm, basis *b, char *dir_output, char *method, int verbose);
void print_hamiltonian(params *pm, hamiltonian *H, char *dir_output, char *method, int verbose);

double lapack_eigval(params *pm, hamiltonian *H); // LAPACK_zheev eigensolver (full matrix)
void lapack_eig(params *pm, hamiltonian *H, double *eigval, double_complex *eigvec); // LAPACK_zheev eigensolver (full matrix)

#ifdef __cplusplus
extern "C" {
#endif
	double arpack_eigval(params *pm, hamiltonian *H); // ARPACK AREig eigensolver (csc matrix)
#ifdef __cplusplus
}
#endif

#endif
