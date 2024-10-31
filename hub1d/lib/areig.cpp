#include <arcomp.h>
#include <arlnsmat.h>
#include <arlscomp.h>
#include "hub1d.h"

extern "C" void arpack_eig(params *pm, hamiltonian *h, double *eigval, double_complex *eigvec) {
	int i, n=pm->Nb, nnz=h->nnz, nev=8;
	arcomplex<double> A[nnz], EigVal[nev], EigVec[nev*n], *EigVal_p=EigVal, *EigVec_p=EigVec;

	for(i=0; i<nnz; i++) A[i] = arcomplex<double>(h->val[i]);

	ARluNonSymMatrix<arcomplex<double>, double> matrix(n, nnz, A, h->row, h->col_csc);
	ARluCompStdEig<double> prob(nev, matrix, "SR");
	prob.EigenValVectors(EigVec_p, EigVal_p);
	
	*eigval = EigVal[0].real();
	for(i=0; i<n; i++) eigvec[i] = EigVec[n*0 + i].real() + I*EigVec[n*0 + i].imag();
}
