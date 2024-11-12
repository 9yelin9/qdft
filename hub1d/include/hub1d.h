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

typedef struct {
	int N; // # of sites
	int Nx; // # of spin-orbitals (2 * N)
	int Ne; // # of electrons
	int Nb; // # of bases (fci: Nx_C_Ne, soft: N)
	double U; // on-site Coulomb interaction
	double beta; // bethe ansatz parameter
} params;

typedef struct {
	int idx; // index in bases
	int val; // | 1up 2up ... Nup 1dn 2dn ... Ndn >
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
void save_basis(params *pm, basis *b, char *method);
void save_hamiltonian(params *pm, hamiltonian *H, char *method);
double gen_beta(double U, int verbose);

void lapack_eig(params *pm, hamiltonian *H, double *eigval, double_complex *eigvec); // LAPACK_zheev eigensolver (full matrix)

#ifdef __cplusplus
extern "C" {
#endif
	void arpack_eig(params *pm, hamiltonian *H, double *eigval, double_complex *eigvec); // ARPACK AREig eigensolver (csc matrix)
#ifdef __cplusplus
}
#endif

#endif
