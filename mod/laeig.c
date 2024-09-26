#include <lapack.h>
#include "hub1d.h"

double laeig(params *pm, hamiltonian *H) {
	int i, ln=pm->Nb, lda=ln, lwork=2*ln-1, info;
	char jobz='N', uplo='U';
	double w[ln], rwork[3*ln-2];
	double_complex a[ln*ln], work[lwork];

	memset(a, 0, sizeof(double_complex) * ln*ln);
	for(i=0; i<H->nnz; i++) a[pm->Nb*H->row[i] + H->col[i]] = H->val[i];

	LAPACK_zheev(&jobz, &uplo, &ln, a, &lda, w, work, &lwork, rwork, &info);
	if(info) {
		printf("ERROR: LAPACK_zheev fail (info = %d)\n", info);
		exit(1);
	}

	return w[0];
}
