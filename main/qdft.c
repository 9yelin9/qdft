#include "qdft.h"

void dec2bin(int len, int dec, char *bin) {
	int i;
	for(i=0; i<len; i++) {
		if((dec >> i) & 1) bin[i] = '1';
		else               bin[i] = '0';
	}
}

void print_H(environment *env, hamiltonian *H) {
	int i;
	char row_bin[env->N], col_bin[env->N];

	printf("--------------------------------- Hamiltonian ---------------------------------\n");
	printf("%8s%8s%8s%16s%16s%16s\n", "i", "row", "col", "row(2)", "col(2)", "val");
	for(i=0; i<H->nnz; i++) {
		dec2bin(env->N, H->row[i], row_bin);
		dec2bin(env->N, H->col[i], col_bin);
		printf("%8d%8d%8d%16s%16s%16f\n", i, H->row[i], H->col[i], row_bin, col_bin, creal(H->val[i]));
	}
	printf("-------------------------------------------------------------------------------\n\n");

	printf("------- check coo2csr -------\n");
	printf("%8s%8s%8s\n", "col_csr", "", "col_coo");
	for(i=0; i<env->M+1; i++) printf("%8d%8s%8d\n", H->col_csr[i], "", i);
	printf("-----------------------------\n\n");
}

void gen_H(environment *env, hamiltonian *H, double U) {
	int n, a, b, i, j, ii, jj, sp, d_ij=1;

	H->nnz = 0;
	for(n=0; n<env->M; n++) {
		H->col_csr[n] = H->nnz;

		H->col[H->nnz] = H->row[H->nnz] = n;
		H->val[H->nnz] = 0;
		for(i=0; i<env->Ne; i++) {
			a = n >> i;
			b = n >> (i + env->Ne);

			H->val[H->nnz] += V(i) * ((a & 1) + (b & 1)) + U * ((a & 1) * (b & 1));
		}
		H->nnz++;

		for(i=0; i<env->Ne; i++) {
			j = (i + d_ij) % env->Ne;
			for(sp=0; sp<2; sp++) {
				ii = i + env->Ne * sp;
				jj = j + env->Ne * sp;

				a = n >> ii;
				b = n >> jj;

				if((!(a & 1) && (b & 1)) || ((a & 1) && !(b & 1))) {
					H->col[H->nnz] = n;
					H->row[H->nnz] = n ^ ((1 << ii) | (1 << jj));
					H->val[H->nnz] = ((i + d_ij) / env->Ne) & 1 ? 1 : -1;
					H->nnz++;
				}
			}
		}
	}
	H->col_csr[n] = H->nnz;
}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("Usage: %s <Ne> <U>\n", argv[0]);
		exit(1);
	}

	environment env = {
		.Ne = atoi(argv[1]),
		.N = 2 * env.Ne,
		.M = 1 << env.N,
	};
	double U = atof(argv[2]);
	printf("Ne=%d\tN=%d\tM=%d\tU=%.1f\n\n", env.Ne, env.N, env.M, U);

	hamiltonian H = {
		.nnz = env.N * env.M, // init nnz
		.row = (int*)malloc(sizeof(int) * H.nnz),
		.col = (int*)malloc(sizeof(int) * H.nnz),
		.col_csr = (int*)malloc(sizeof(int) * (env.M + 1)),
		.val = (double_complex*)malloc(sizeof(double_complex) * H.nnz),
	};
	gen_H(&env, &H, U); print_H(&env, &H);

	free(H.row);
	free(H.col);
	free(H.col_csr);
	free(H.val);

	return 0;
}
