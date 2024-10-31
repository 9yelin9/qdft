#include <lapack.h>
#include "hub1d.h"

void lapack_eig(params *pm, hamiltonian *h, double *eigval, double_complex *eigvec) {
	int i, ln=pm->Nb, lda=ln, lwork=2*ln-1, info;
	char jobz='V', uplo='U';
	double rwork[3*ln-2];
	double_complex work[lwork];

	memset(eigvec, 0, sizeof(double_complex) * ln*ln);
	for(i=0; i<h->nnz; i++) eigvec[pm->Nb*h->row[i] + h->col[i]] = h->val[i];

	LAPACK_zheev(&jobz, &uplo, &ln, eigvec, &lda, eigval, work, &lwork, rwork, &info);
	if(info) {
		printf("\nERROR: LAPACK_zheev fail (info = %d)\n", info);
		exit(1);
	}
}
