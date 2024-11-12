#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include "hub1d.h"

#define BETHE_ANSATZ(x)       (-(2. * x / M_PI) * sin(M_PI / x))
#define BETHE_ANSATZ_DERIV(x) (-(2. / M_PI) * sin(M_PI / x) + (2. / x) * cos(M_PI / x))

#define METHOD "soft"

#define V(i) (i / 10.)

#define E_H(occ, U)       (U * (occ / 2.) * (occ / 2.))
#define E_H_DERIV(occ, U) (U * (occ / 2.))

#define E_BA(occ, beta)       (-(2 * beta / M_PI) * sin(M_PI * occ / beta))
#define E_BA_DERIV(occ, beta) (-2 * cos(M_PI * occ / beta))

#define E_XC(occ, U, beta)       (E_BA(occ, beta) - E_BA(occ, 2.) - E_H(occ, U))
#define E_XC_DERIV(occ, U, beta) (E_BA_DERIV(occ, beta) - E_BA_DERIV(occ, 2.) - E_H_DERIV(occ, U))

double bethe_ansatz_r_integrand(double x, void *v) {
	double U=*(double*)v, num, den;
	num = gsl_sf_bessel_J0(x) * gsl_sf_bessel_J1(x);
	den = x * (1. + exp(U * x / 2.));
	return num / den;
}

double bethe_ansatz(double x, void *v) {
	double bethe_ansatz_l, bethe_ansatz_r;
	bethe_ansatz_l = BETHE_ANSATZ(x);
	bethe_ansatz_r = *(double*)v;
	return bethe_ansatz_l - bethe_ansatz_r;
}

double bethe_ansatz_deriv(double x, void *v) {
	double bethe_ansatz_deriv_l, bethe_ansatz_deriv_r;
	bethe_ansatz_deriv_l = BETHE_ANSATZ_DERIV(x); 
	bethe_ansatz_deriv_r = 0.;
	return bethe_ansatz_deriv_l - bethe_ansatz_deriv_r;
}

void bethe_ansatz_fdf(double x, void *v, double *y, double *dy) {
	double bethe_ansatz_l, bethe_ansatz_r;
	bethe_ansatz_l = BETHE_ANSATZ(x);
	bethe_ansatz_r = *(double*)v;
	*y = bethe_ansatz_l - bethe_ansatz_r;

	double bethe_ansatz_deriv_l, bethe_ansatz_deriv_r;
	bethe_ansatz_deriv_l = BETHE_ANSATZ_DERIV(x);
	bethe_ansatz_deriv_r = 0.;
	*dy = bethe_ansatz_deriv_l - bethe_ansatz_deriv_r;
}

double gen_beta(double U, int verbose) {
	double bethe_ansatz_r, abserr;
	gsl_function F = {
		.function = &bethe_ansatz_r_integrand,
		.params = &U,
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

	if(!verbose) freopen("/dev/null", "w", stdout);

	printf("---------- gen_beta (secant method) ----------\n");
	printf("%8s%16s%16s\n", "itr", "beta", "delta");

	int itr, itr_max=100, status;
	for(itr=0; itr<itr_max; itr++) {
		status = gsl_root_fdfsolver_iterate(s);
		beta = gsl_root_fdfsolver_root(s);
		status = gsl_root_test_delta(beta, beta_old, 0, 1e-6);

		printf("%8d%16f%16f\n", itr, beta, beta - beta_old);

		if(status == GSL_SUCCESS) break;
		else beta_old = beta;
	}
	gsl_root_fdfsolver_free(s);
	printf("----------------------------------------------\n\n");

	if(!verbose) freopen("/dev/tty", "w", stdout);

	return beta;
}

void gen_basis_soft(params *pm, basis *b) {
	int n; 
	for(n=0; n<pm->Nb; n++) {
		b[n].idx = n;
		b[n].val = 1 << n;
	}
}

void gen_hamiltonian_soft(params *pm, basis *b, hamiltonian *h, double *occ) {
	int n, p, q, i, j, key, d_ij=1;
	basis *find;

	h->nnz = 0;
	for(n=0; n<pm->Nb; n++) {
		i = __builtin_ctz(b[n].val);

		h->row[h->nnz] = n;
		h->col[h->nnz] = n;
		h->val[h->nnz] = V(i) + E_XC_DERIV(occ[i], pm->U, pm->beta) + E_H_DERIV(occ[i], pm->U);
		//h->val[h->nnz] = 0;
		h->nnz++;

		for(i=0; i<pm->N; i++) {
			j = (i + d_ij) % pm->N;

			p = (b[n].val >> i) & 1;
			q = (b[n].val >> j) & 1;

			if(p ^ q) {
				key = b[n].val ^ ((1 << i) | (1 << j));
				find = bsearch(&key, b, pm->Nb, sizeof(basis), compare_basis);

				if(find != NULL) {
					h->row[h->nnz] = n;
					h->col[h->nnz] = find->idx;
					h->val[h->nnz] = ((i + d_ij) / pm->N) & 1 ? 1 : -1;
					h->nnz++;
				}
				else {
					printf("\nERROR: bsearch fail (n = %d)\n", n);
					exit(1);
				}
			}
		}
	}
}

void run_soft(params *pm, basis *b, hamiltonian *h, int verbose) {
	int i, k, itr, itr_max=100;
	double e_grd_old=100, occ_mix=0.1, occ_tot, occ[pm->N], occ_old[pm->N], eigval[pm->Nb];
	double_complex eigvec[pm->Nb*pm->Nb];

	if(!verbose) freopen("/dev/null", "w", stdout);

	printf("------------------------------------------------ kohn-sham self-consistent field ------------------------------------------------\n");
	printf("%8s%12s", "itr", "e_grd"); for(i=0; i<pm->N; i++) printf("%10s%02d", "occ", i); printf("\n");

	for(i=0; i<pm->N; i++) occ[i] = (double)pm->Ne / pm->N; // init occ
	printf("%8d%12s", 0, "-"); for(i=0; i<pm->N; i++) printf("%12f", occ[i]); printf("\n");

	for(itr=1; itr<=itr_max; itr++) {
		gen_hamiltonian_soft(pm, b, h, occ);
		lapack_eig(pm, h, eigval, eigvec);

		for(i=0; i<pm->N; i++) occ_old[i] = occ[i];
		memset(occ, 0, sizeof(occ));
		for(k=0; k<pm->Ne; k++) for(i=0; i<pm->N; i++) occ[i] += square_complex(eigvec[pm->Nb*(k/2) + i]);

		h->e_grd = 0;
		for(k=0; k<pm->Ne; k++) h->e_grd += eigval[k/2];
		for(i=0; i<pm->N;  i++) h->e_grd += E_XC(occ[i], pm->U, pm->beta) + E_H(occ[i], pm->U);
		for(i=0; i<pm->N;  i++) h->e_grd -= (E_XC_DERIV(occ_old[i], pm->U, pm->beta) + E_H_DERIV(occ_old[i], pm->U)) * occ[i];

		printf("%8d%12f", itr, h->e_grd); for(i=0; i<pm->N; i++) printf("%12f", occ[i]); printf("\n");

		occ_tot = 0;
		for(i=0; i<pm->N; i++) occ_tot += occ[i];
		if(fabs(occ_tot - pm->Ne) > 1e-6) {
			printf("\nERROR: occ_tot(%f) != Ne(%d)\n", occ_tot, pm->Ne);
			exit(1);
		}

		if(fabs(h->e_grd - e_grd_old) < 1e-6) break;

		e_grd_old = h->e_grd;
		for(i=0; i<pm->N; i++) occ[i] = (1 - occ_mix) * occ_old[i] + occ_mix * occ[i];
	}
	printf("---------------------------------------------------------------------------------------------------------------------------------\n\n");

	if(!verbose) freopen("/dev/tty", "w", stdout);
}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("Usage: %s <N> <Ne> <U> [verbose]\n\n", argv[0]);
		exit(1);
	}

	params pm = {
		.N = atoi(argv[1]),
		.Nx = pm.N,
		.Ne = atoi(argv[2]),
		.Nb = pm.N,
		.U = atof(argv[3]),
	};
	int verbose = argv[4] == NULL ? 0 : 1;
	pm.beta = gen_beta(pm.U, verbose);
	printf("N = %d\nNe = %d\nNb = %d\nU = %f\nbeta = %f\n\n", pm.N, pm.Ne, pm.Nb, pm.U, pm.beta);

	basis b[pm.N];
	gen_basis_soft(&pm, b);
	save_basis(&pm, b, METHOD);

	hamiltonian h = {
		.nnz = pm.N * pm.N,
		.row = (int*)malloc(sizeof(int) * h.nnz),
		.col = (int*)malloc(sizeof(int) * h.nnz),
		.val = (double_complex*)malloc(sizeof(double_complex) * h.nnz),
	};
	run_soft(&pm, b, &h, 1);
	save_hamiltonian(&pm, &h, METHOD);

	free(h.row);
	free(h.col);
	free(h.col_csc);
	free(h.val);

	return 0;
}
