#include "qdft.h"

#define V(i) (i / 10.)

typedef struct {
	int row;
	double_complex val;
} hamiltonian_tmp;

int cnt_one(int n) {
	int cnt=0;
	while(n) {
		cnt += n & 1;
		n >>= 1;
	}
	return cnt;
}

int compare(const void *a, const void *b) {
	int a_row = ((hamiltonian_tmp*)a)->row;
	int b_row = ((hamiltonian_tmp*)b)->row;
	return (a_row > b_row) - (a_row < b_row);
}

void gen_H(environment *env, hamiltonian *H) {
	int n, a, b, i, j, ii, jj, sp, d_ij=1, nnz_tmp;
	hamiltonian_tmp H_tmp[env->N];

	H->nnz = 0;
	for(n=0; n<env->M; n++) {
		if(cnt_one(n) == env->Ne) {
			H->col_csc[n] = H->nnz;

			nnz_tmp = 0;
			H_tmp[nnz_tmp].row = n;
			H_tmp[nnz_tmp].val = 0;
			for(i=0; i<env->Ne; i++) {
				a = n >> i;
				b = n >> (i + env->Ne);

				H_tmp[nnz_tmp].val += V(i) * ((a & 1) + (b & 1)) + env->U * ((a & 1) * (b & 1));
			}
			nnz_tmp++;

			for(i=0; i<env->Ne; i++) {
				j = (i + d_ij) % env->Ne;
				for(sp=0; sp<2; sp++) {
					ii = i + env->Ne * sp;
					jj = j + env->Ne * sp;

					a = n >> ii;
					b = n >> jj;

					if((!(a & 1) && (b & 1)) || ((a & 1) && !(b & 1))) {
						H_tmp[nnz_tmp].row = n ^ ((1 << ii) | (1 << jj));
						H_tmp[nnz_tmp].val = ((i + d_ij) / env->Ne) & 1 ? 1 : -1;
						nnz_tmp++;
					}
				}
			}

			qsort(H_tmp, nnz_tmp, sizeof(hamiltonian_tmp), compare);
			for(i=0; i<nnz_tmp; i++) {
				H->col[H->nnz] = n;
				H->row[H->nnz] = H_tmp[i].row;
				H->val[H->nnz] = H_tmp[i].val;
				H->nnz++;
			}
		}
	}
	H->col_csc[n] = H->nnz;
}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("Usage: %s <Ni> <Ne> <U>\n", argv[0]);
		exit(1);
	}

	environment env = {
		.Ni = atoi(argv[1]),
		.Ne = atoi(argv[2]),
		.N = 2 * env.Ni,
		.M = 1 << env.N,
		.U = atof(argv[3]),
	};
	printf("Ni=%d\tNe=%d\tN=%d\tM=%d\tU=%f\n\n", env.Ni, env.Ne, env.N, env.M, env.U);

	hamiltonian H = {
		.nnz = env.N * env.M, // init nnz
		.row = (int*)malloc(sizeof(int) * H.nnz),
		.col = (int*)malloc(sizeof(int) * H.nnz),
		.col_csc = (int*)malloc(sizeof(int) * (env.M + 1)),
		.val = (double_complex*)malloc(sizeof(double_complex) * H.nnz),
	};

	gen_H(&env, &H);
	print_H(&env, &H);
	printf("FCI(LAPACK): %f\n", laeig(&env, &H));
	printf("FCI(ARPACK): %f\n", areig(&env, &H));

	free(H.row);
	free(H.col);
	free(H.col_csc);
	free(H.val);

	return 0;
}
