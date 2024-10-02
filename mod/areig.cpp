#include <arcomp.h>
#include <arlnsmat.h>
#include <arlscomp.h>
#include "hub1d.h"

extern "C" double arpack_eigval(params *pm, hamiltonian *h) {
	int i, n=pm->Nb, nnz=h->nnz, nev=8;
	arcomplex<double> A[nnz], EigVal[nev], *EigVal_p=EigVal;

	for(i=0; i<nnz; i++) A[i] = arcomplex<double>(h->val[i]);

	ARluNonSymMatrix<arcomplex<double>, double> matrix(n, nnz, A, h->row, h->col_csc);
	ARluCompStdEig<double> prob(nev, matrix, "SR");
	prob.Eigenvalues(EigVal_p);

	return EigVal[nev-1].real();
}
