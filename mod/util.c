#include "qdft.h"

void dec2bin(int len, int dec, char *bin) {
	int i;
	for(i=0; i<len; i++) {
		if((dec >> i) & 1) bin[i] = '1';
		else               bin[i] = '0';
	}
	bin[i] = '\0';
}

void print_H(environment *env, hamiltonian *H) {
	int i, j, is_sym=1;
	char row_bin[env->N], col_bin[env->N];

	printf("--------------------------------- Hamiltonian ---------------------------------\n");
	printf("%8s%8s%8s%16s%16s%16s\n", "i", "row", "col", "row(2)", "col(2)", "val");
	for(i=0; i<H->nnz; i++) {
		dec2bin(env->N, H->row[i], row_bin);
		dec2bin(env->N, H->col[i], col_bin);

		for(j=0; j<H->nnz; j++) {
			if(((H->row[i] == H->col[j]) && (H->row[j] == H->col[i])) && (fabs(H->val[i] - H->val[j]) > 1e-6)) {
				is_sym = 0;
				break;
			}
		}
		if(is_sym) printf("%8d%8d%8d%16s%16s%16f\n", i, H->row[i], H->col[i], row_bin, col_bin, creal(H->val[i]));
		else {
			printf("\nBroken symmetry:\n");
			printf("%8d%8d%8d%16s%16s%16f\n",   i, H->row[i], H->col[i], row_bin, col_bin, creal(H->val[i]));
			printf("%8d%8d%8d%16s%16s%16f\n\n", j, H->row[j], H->col[j], row_bin, col_bin, creal(H->val[j]));
			break;
		}
	}
	printf("-------------------------------------------------------------------------------\n\n");

	if(is_sym) {
		printf("------- check coo2csc -------\n");
		printf("%8s%8s%8s\n", "col_csc", "", "col_coo");
		for(i=0; i<env->M+1; i++) printf("%8d%8s%8d\n", H->col_csc[i], "", i);
		printf("-----------------------------\n\n");
	}
}

