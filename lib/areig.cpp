#include <arcomp.h>
#include "areig.h"
#include "acompsol.h"
#include "qdft.h"

/*
int acompreg(environment *env, hamiltonian *H) {
	int n=env->M*env->M, nnz=env->H_len, nconv;

	nconv = AREig(EigVal, EigVec, n, nnz, A, irow, pcol, nev);
	Solution(nconv, n, nnz, A, irow, pcol, EigVal, EigVec);

	return nconv < nev ? EXIT_FAILURE : EXIT_SUCCESS;
}
*/
