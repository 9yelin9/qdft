#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include "hub1d.h"

#define METHOD "soft"
#define E_BA(occ, beta) (-(2 * beta / M_PI) * sin(M_PI * occ / beta))

double bethe_ansatz_r_integrand(double x, void *v) {
	double U=*(double*)v, num, den;
	num = gsl_sf_bessel_J0(x) * gsl_sf_bessel_J1(x);
	den = x * (1. + exp(U * x / 2.));
	return num / den;
}

double bethe_ansatz(double x, void *v) {
	double bethe_ansatz_l, bethe_ansatz_r;
	bethe_ansatz_l = -(2. * x / M_PI) * sin(M_PI / x);
	bethe_ansatz_r = *(double*)v;
	return bethe_ansatz_l - bethe_ansatz_r;
}

double bethe_ansatz_deriv(double x, void *v) {
	double bethe_ansatz_deriv_l, bethe_ansatz_deriv_r;
	bethe_ansatz_deriv_l = -(2. / M_PI) * sin(M_PI / x) + (2. / x) * cos(M_PI / x);
	bethe_ansatz_deriv_r = 0.;
	return bethe_ansatz_deriv_l - bethe_ansatz_deriv_r;
}

void bethe_ansatz_fdf(double x, void *v, double *y, double *dy) {
	double bethe_ansatz_l, bethe_ansatz_r;
	bethe_ansatz_l = -(2. * x / M_PI) * sin(M_PI / x);
	bethe_ansatz_r = *(double*)v;
	*y = bethe_ansatz_l - bethe_ansatz_r;

	double bethe_ansatz_deriv_l, bethe_ansatz_deriv_r;
	bethe_ansatz_deriv_l = -(2. / M_PI) * sin(M_PI / x) + (2. / x) * cos(M_PI / x);
	bethe_ansatz_deriv_r = 0.;
	*dy = bethe_ansatz_deriv_l - bethe_ansatz_deriv_r;
}

double gen_beta(params *pm) {
	double bethe_ansatz_r, abserr;
	gsl_function F = {
		.function = &bethe_ansatz_r_integrand,
		.params = &(pm->U),
	};
	gsl_integration_workspace *w=gsl_integration_workspace_alloc(1024);
	gsl_integration_qagiu(&F, 0., 0., 1e-4, 1024, w, &bethe_ansatz_r, &abserr);
	gsl_integration_workspace_free(w);
	bethe_ansatz_r *= -4;

	double beta=1, beta_old=1;
	gsl_function_fdf FDF = {
		.f = &bethe_ansatz,
		.df = &bethe_ansatz_deriv,
		.fdf = &bethe_ansatz_fdf,
		.params = &bethe_ansatz_r,
	};
	const gsl_root_fdfsolver_type *T=gsl_root_fdfsolver_secant;
	gsl_root_fdfsolver *s=gsl_root_fdfsolver_alloc(T);
	gsl_root_fdfsolver_set(s, &FDF, beta);

	if(!pm->verbose) freopen("/dev/null", "w", stdout);

	printf("---------- gen_beta (secant method) ----------\n");
	printf("%8s%16s%16s\n", "itr", "beta", "delta");

	int itr, max_itr=100, status;
	for(itr=0; itr<max_itr; itr++) {
		status = gsl_root_fdfsolver_iterate(s);
		beta = gsl_root_fdfsolver_root(s);
		status = gsl_root_test_delta(beta, beta_old, 0, 1e-6);

		printf("%8d%16f%16f\n", itr, beta, beta - beta_old);

		if(status == GSL_SUCCESS) break;
		else beta_old = beta;
	}
	gsl_root_fdfsolver_free(s);
	printf("----------------------------------------------\n\n");

	if(!pm->verbose) freopen("/dev/tty", "w", stdout);

	return beta;
}

void gen_basis_soft(params *pm, basis *b) {
	int n; 
	for(n=0; n<pm->Nb; n++) {
		b[n].idx = n;
		b[n].val = 1 << n;
	}
}

void gen_H_soft(params *pm, hamiltonian *H, double *occ, double beta) {
	/*
	int i, j, nnz_tmp=0;
	double e_xc[pm->Ni];

	for(i=0; i<pm->Ni; i++) e_xc[i] = E_BA(occ[i], beta) - E_BA(occ[i], 2.) - U * (occ[i] / 2.)*(occ[i] / 2.);

	for(i=0; i<pm->Ni; i++) {
		for(j=0; j<pm->Ni; j++) {
			H->row[nnz_tmp] = i;
			H->col[nnz_tmp] = j;
			H->val[nnz_tmp] = ;
			nnz_tmp++;
		}
	}
	*/
}

void run_soft(params *pm, hamiltonian *H) {
	int i, itr, max_itr=100;
	double beta=gen_beta(pm), occ[pm->Ni];

	for(i=0; i<pm->Ni; i++) occ[i] = pm->Ne / pm->Ni; // init occ

	for(itr=0; itr<max_itr; itr++) {
		gen_H_soft(pm, H, occ, beta);
	}
}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("Usage: %s <Ni> <Ne> <U> [verbose=0/1]\n\n", argv[0]);
		exit(1);
	}

	params pm = {
		.Ni = atoi(argv[1]),
		.Ne = atoi(argv[2]),
		.N  = pm.Ni,
		.Nb = pm.Ni,
		.U = atof(argv[3]),
		.verbose = argv[4] == NULL ? 0 : 1,
	};
	printf("Ni=%d\nNe=%d\nU=%f\n\n", pm.Ni, pm.Ne, pm.U);

	char dir_output[1024];
	gen_dir_output(&pm, dir_output);

	basis b[pm.Ni];
	gen_basis_soft(&pm, b);
	print_basis(&pm, b, dir_output, METHOD);

	hamiltonian H = {
		.nnz = pm.Ni * pm.Ni,
		.row = (int*)malloc(sizeof(int) * H.nnz),
		.col = (int*)malloc(sizeof(int) * H.nnz),
		.val = (double_complex*)malloc(sizeof(double_complex) * H.nnz),
	};

	run_soft(&pm, &H);
	print_H(&pm, &H, dir_output, METHOD);

	free(H.row);
	free(H.col);
	free(H.col_csc);
	free(H.val);

	return 0;
}
