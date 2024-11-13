#include "hub1d.h"

#define METHOD "fci"

#define V(i) (i / 10.)

typedef struct {
	int row;
	double_complex val;
} hamiltonian_n;

int combination(int n, int r) {
	int i;
	long num=1, den=1;
	for(i=1; i<=r; i++) {
		num *= n - r + i;
		den *= i;
	}
	return num / den;
}

int compare_hamiltonian_n(const void *p, const void *q) {
	int p_row = ((hamiltonian_n*)p)->row;
	int q_row = ((hamiltonian_n*)q)->row;
	return (p_row > q_row) - (p_row < q_row);
}

void gen_basis_fci(params *pm, basis *b) {
	int i, cnt=0;
	for(i=0; i<(1 << pm->Nx); i++) {
		if(__builtin_popcount(i) == pm->Ne) {
			b[cnt].idx = cnt;
			b[cnt].val = i;
			cnt++;
		}
	}
	if(cnt != pm->Nb) {
		printf("\nERROR: Too many/few bases (%d != %d)\n", cnt, pm->Nb);
		exit(1);
	}
}

void gen_hamiltonian_fci(params *pm, basis *b, hamiltonian *h) {
	int n, p, q, i, j, ii, jj, sp, key, nnz_n, d_ij=1;
	basis *find;
	hamiltonian_n h_n[h->nnz];

	h->nnz = 0;
	for(n=0; n<pm->Nb; n++) {
		h->col_csc[n] = h->nnz;
		nnz_n = 0;

		h_n[nnz_n].row = n;
		h_n[nnz_n].val = 0;
		for(i=0; i<pm->N; i++) {
			p = (b[n].val >> i) & 1;
			q = (b[n].val >> (i + pm->N)) & 1;
			h_n[nnz_n].val += V(i) * (p + q) + pm->U * (p * q);
			//h_n[nnz_n].val += 0; 
		}
		nnz_n++;

		for(i=0; i<pm->N; i++) {
			j = (i + d_ij) % pm->N;
			for(sp=0; sp<2; sp++) {
				ii = i + pm->N * sp;
				jj = j + pm->N * sp;

				p = (b[n].val >> ii) & 1;
				q = (b[n].val >> jj) & 1;

				if(p ^ q) {
					key = b[n].val ^ ((1 << ii) | (1 << jj));
					find = bsearch(&key, b, pm->Nb, sizeof(basis), compare_basis);

					if(find != NULL) {
						h_n[nnz_n].row = find->idx;
						h_n[nnz_n].val = ((i + d_ij) / pm->N) & 1 ? 1 : -1;
						nnz_n++;
					}
					else {
						printf("\nERROR: bsearch fail (n = %d)\n", n);
						exit(1);
					}
				}
			}
		}

		qsort(h_n, nnz_n, sizeof(hamiltonian_n), compare_hamiltonian_n);
		for(i=0; i<nnz_n; i++) {
			h->col[h->nnz] = n;
			h->row[h->nnz] = h_n[i].row;
			h->val[h->nnz] = h_n[i].val;
			h->nnz++;
		}
	}
	h->col_csc[n] = h->nnz;
}

void run_fci(params *pm, basis *b, hamiltonian *h) {
	int i, n;
	double occ_tot, occ[pm->N], eigval[pm->Nb];
	double_complex eigvec[pm->Nb*pm->Nb];

	gen_hamiltonian_fci(pm, b, h);
	lapack_eig(pm, h, eigval, eigvec);
	h->e_grd  = eigval[0];
	for(i=0; i<16; i++) printf("%.2f\t", eigval[i]);
	printf("\n");

	save_hamiltonian(pm, h, METHOD);

	printf("------------------------------------------------ full configuration interaction ------------------------------------------------\n");
	printf("%8s%12s", " ", "e_grd"); for(i=0; i<pm->N; i++) printf("%10s%02d", "occ", i); printf("\n");

	memset(occ, 0, sizeof(occ));
	for(n=0; n<pm->Nb; n++) for(i=0; i<pm->Nx; i++) occ[i % pm->N] += square_complex(eigvec[n]) * ((b[n].val >> i) & 1);

	occ_tot = 0;
	for(i=0; i<pm->N; i++) occ_tot += occ[i];
	if(fabs(occ_tot - pm->Ne) > 1e-6) {
		printf("\nERROR: occ_tot(%f) != Ne(%d)\n", occ_tot, pm->Ne);
		exit(1);
	}

	printf("%8s%12f", " ", h->e_grd); for(i=0; i<pm->N; i++) printf("%12f", occ[i]); printf("\n");
	printf("--------------------------------------------------------------------------------------------------------------------------------\n\n");
}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("Usage: %s <N> <Ne> <U>\n\n", argv[0]);
		exit(1);
	}

	params pm = {
		.N = atoi(argv[1]),
		.Nx = 2 * pm.N,
		.Ne = atoi(argv[2]),
		.Nb = combination(pm.Nx, pm.Ne),
		.U = atof(argv[3]),
	};
	printf("N = %d\nNe = %d\nNx = %d\nNb = %d\nU = %f\n\n", pm.N, pm.Ne, pm.Nx, pm.Nb, pm.U);

	basis b[pm.Nb];
	gen_basis_fci(&pm, b);
	save_basis(&pm, b, METHOD);

	hamiltonian h = {
		.nnz = pm.Nb * pm.Nb, // init nnz
		.row = (int*)malloc(sizeof(int) * h.nnz),
		.col = (int*)malloc(sizeof(int) * h.nnz),
		.col_csc = (int*)malloc(sizeof(int) * (pm.Nb + 1)),
		.val = (double_complex*)malloc(sizeof(double_complex) * h.nnz),
	};
	run_fci(&pm, b, &h);

	free(h.row);
	free(h.col);
	free(h.col_csc);
	free(h.val);

	return 0;
}
