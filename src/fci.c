#include "hub1d.h"

#define METHOD "fci"
#define V(i) (i / 10.)

typedef struct {
	int row;
	double_complex val;
} hamiltonian_tmp;

int combination(int n, int r) {
	int i;
	long num=1, den=1;
	for(i=1; i<=r; i++) {
		num *= n - r + i;
		den *= i;
	}
	return num / den;
}

int cnt_ones(int n) {
	int cnt=0;
	while(n) {
		cnt += n & 1;
		n >>= 1;
	}
	return cnt;
}

int compare_hamiltonian_tmp(const void *p, const void *q) {
	int p_row = ((hamiltonian_tmp*)p)->row;
	int q_row = ((hamiltonian_tmp*)q)->row;
	return (p_row > q_row) - (p_row < q_row);
}

int compare_basis(const void *key, const void *p) {
	int key_val = *(int*)key;
	int p_val = ((basis*)p)->val;
	return (key_val > p_val) - (key_val < p_val);
}

void gen_basis_fci(params *pm, basis *b) {
	int n, cnt=0;
	for(n=0; n<(1 << pm->N); n++) {
		if(cnt_ones(n) == pm->Ne) {
			b[cnt].idx = cnt;
			b[cnt].val = n;
			cnt++;
		}
	}
	if(cnt != pm->Nb) {
		printf("ERROR: Too many/few bases (%d != %d)\n", cnt, pm->Nb);
		exit(1);
	}
}

void gen_H_fci(params *pm, basis *b, hamiltonian *H) {
	int n, p, q, i, j, ii, jj, d_ij=1, sp, key, nnz_tmp;
	basis *find;
	hamiltonian_tmp H_tmp[pm->N];

	H->nnz = 0;
	for(n=0; n<pm->Nb; n++) {
		H->col_csc[n] = H->nnz;

		nnz_tmp = 0;
		H_tmp[nnz_tmp].row = n;
		H_tmp[nnz_tmp].val = 0;
		for(i=0; i<pm->Ni; i++) {
			p = b[n].val >> i;
			q = b[n].val >> (i + pm->Ni);

			H_tmp[nnz_tmp].val += V(i) * ((p & 1) + (q & 1)) + pm->U * ((p & 1) * (q & 1));
		}
		nnz_tmp++;

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
						H_tmp[nnz_tmp].row = find->idx;
						H_tmp[nnz_tmp].val = ((i + d_ij) / pm->Ni) & 1 ? 1 : -1;
						nnz_tmp++;
					}
					else {
						printf("ERROR: bsearch fail (n = %d)\n", n);
						exit(1);
					}
				}
			}
		}

		qsort(H_tmp, nnz_tmp, sizeof(hamiltonian_tmp), compare_hamiltonian_tmp);
		for(i=0; i<nnz_tmp; i++) {
			H->col[H->nnz] = n;
			H->row[H->nnz] = H_tmp[i].row;
			H->val[H->nnz] = H_tmp[i].val;
			H->nnz++;
		}
	}
	H->col_csc[n] = H->nnz;
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
		.verbose = argv[4] == NULL ? 0 : 1,
	};
	printf("Ni=%d\nNe=%d\nN=%d\nNb=%d\nU=%f\n\n", pm.Ni, pm.Ne, pm.N, pm.Nb, pm.U);

	char dir_output[1024];
	gen_dir_output(&pm, dir_output);

	basis b[pm.Nb];
	gen_basis_fci(&pm, b);
	print_basis(&pm, b, dir_output, METHOD);

	hamiltonian H = {
		.nnz = pm.Nb * pm.Nb, // init nnz
		.row = (int*)malloc(sizeof(int) * H.nnz),
		.col = (int*)malloc(sizeof(int) * H.nnz),
		.col_csc = (int*)malloc(sizeof(int) * (pm.Nb + 1)),
		.val = (double_complex*)malloc(sizeof(double_complex) * H.nnz),
	};
	gen_H_fci(&pm, b, &H);
	print_H(&pm, &H, dir_output, METHOD);

	printf("\n");
	printf("FCI(LAPACK): %f\n", laeig(&pm, &H));
	printf("FCI(ARPACK): %f\n", areig(&pm, &H));

	free(H.row);
	free(H.col);
	free(H.col_csc);
	free(H.val);

	return 0;
}
