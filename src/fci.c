#include "hub1d.h"

#define METHOD "fci"
#define V(i) (i / 10.)

typedef struct {
	int row;
	double_complex val;
} hamiltonian_n;

typedef struct {
	int idx;
	double val;
} probability;

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

int compare_probability(const void *p, const void *q) {
	double p_val = ((probability*)p)->val;
	double q_val = ((probability*)q)->val;
	return (q_val > p_val) - (q_val < p_val);
}

void gen_basis_fci(params *pm, basis *b) {
	int n, cnt=0;
	for(n=0; n<(1 << pm->N); n++) {
		if(__builtin_popcount(n) == pm->Ne) {
			b[cnt].idx = cnt;
			b[cnt].val = n;
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
	hamiltonian_n h_n[pm->N];

	h->nnz = 0;
	for(n=0; n<pm->Nb; n++) {
		h->col_csc[n] = h->nnz;
		nnz_n = 0;

		h_n[nnz_n].row = n;
		h_n[nnz_n].val = 0;
		for(i=0; i<pm->Ni; i++) {
			p = b[n].val >> i;
			q = b[n].val >> (i + pm->Ni);
			h_n[nnz_n].val += V(i) * ((p & 1) + (q & 1)) + pm->U * ((p & 1) * (q & 1));
		}
		nnz_n++;

		for(i=0; i<pm->Ni; i++) {
			j = (i + d_ij) % pm->Ni;
			for(sp=0; sp<2; sp++) {
				ii = i + pm->Ni * sp;
				jj = j + pm->Ni * sp;

				p = b[n].val >> ii;
				q = b[n].val >> jj;

				if((!(p & 1) && (q & 1)) || ((p & 1) && !(q & 1))) {
					key = b[n].val ^ ((1 << ii) | (1 << jj));
					find = bsearch(&key, b, pm->Nb, sizeof(basis), compare_basis);

					if(find != NULL) {
						h_n[nnz_n].row = find->idx;
						//h_n[nnz_n].val = ((i + d_ij) / pm->Ni) & 1 ? 1 : -1;
						h_n[nnz_n].val = -1;
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

void run_fci(params *pm, basis *b, hamiltonian *h, int verbose) {
	int i, idx_max=10;
	char buf[pm->N+2];
	double_complex eigvec[pm->Nb];
	probability prob[pm->Nb];

	if(!verbose) freopen("/dev/null", "w", stdout);

	printf("--------- full configuration interaction ---------\n");

	gen_hamiltonian_fci(pm, b, h);
	arpack_eig(pm, h, &(h->e_grd), eigvec);

	printf("e_grd = %f\n", h->e_grd);
	printf("%8s%12s%22s\n", "order", "prob", "basis(2)");

	for(i=0; i<pm->Nb; i++) {
		prob[i].idx = i;
		prob[i].val = square_complex(eigvec[i]);
	}
	qsort(prob, pm->Nb, sizeof(probability), compare_probability);

	for(i=0; i<idx_max; i++) {
		dec2bin(pm->N, b[prob[i].idx].val, buf);
		printf("%8d%12f%22s\n", i, prob[i].val, buf);
	}
	printf("--------------------------------------------------\n");

	if(!verbose) freopen("/dev/tty", "w", stdout);
}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("Usage: %s <Ni> <Ne> <U> [verbose=0/1]\n\n", argv[0]);
		exit(1);
	}

	params pm = {
		.Ni = atoi(argv[1]),
		.Ne = atoi(argv[2]),
		.N = 2 * pm.Ni,
		.Nb = combination(pm.N, pm.Ne),
		.U = atof(argv[3]),
	};
	int verbose = argv[4] == NULL ? 0 : 1;
	printf("Ni = %d\nNe = %d\nN = %d\nNb = %d\nU = %f\n\n", pm.Ni, pm.Ne, pm.N, pm.Nb, pm.U);

	char dir_output[1024];
	gen_dir_output(&pm, dir_output);

	basis b[pm.Nb];
	gen_basis_fci(&pm, b);
	print_basis(&pm, b, dir_output, METHOD, verbose);

	hamiltonian h = {
		.nnz = pm.Nb * pm.Nb, // init nnz
		.row = (int*)malloc(sizeof(int) * h.nnz),
		.col = (int*)malloc(sizeof(int) * h.nnz),
		.col_csc = (int*)malloc(sizeof(int) * (pm.Nb + 1)),
		.val = (double_complex*)malloc(sizeof(double_complex) * h.nnz),
	};
	run_fci(&pm, b, &h, 1);
	print_hamiltonian(&pm, &h, dir_output, METHOD, verbose);

	free(h.row);
	free(h.col);
	free(h.col_csc);
	free(h.val);

	return 0;
}
