#include "hub1d.h"

void dec2bin(int len, int dec, char *bin) {
	int i;
	for(i=0; i<len; i++) {
		if((dec >> i) & 1) bin[i] = '1';
		else               bin[i] = '0';
	}
	bin[i] = '\0';
}

void gen_dir_output(params *pm, char *dir_output) {
	sprintf(dir_output, "output/%s_Ni%d_Ne%d", SYSTEM, pm->Ni, pm->Ne);
	mkdir(dir_output, 0755);
}

void print_basis(params *pm, basis *b, char *dir_output, char *method) {
	char fn[1024];
	sprintf(fn, "%s/%s_basis_U%.1f.txt", dir_output, method, pm->U);
	FILE *f = fopen(fn, "w");

	int i;
	char val_bin[pm->N+1];

	if(!pm->verbose) freopen("/dev/null", "w", stdout);

	printf("---------------- basis ----------------\n");
	printf("%8s%8s%20s\n", "idx", "val(10)", "val(2)");
	for(i=0; i<pm->Nb; i++) {
		dec2bin(pm->N, b[i].val, val_bin);
		printf("%8d%8d%20s\n", b[i].idx, b[i].val, val_bin);
		fprintf(f, "%8d%8d%20s\n", b[i].idx, b[i].val, val_bin);
	}
	printf("---------------------------------------\n");
	fclose(f);

	if(!pm->verbose) freopen("/dev/tty", "w", stdout);

	printf("%s: %s\n\n", __func__, fn);
}

void print_H(params *pm, hamiltonian *H, char *dir_output, char *method) {
	char fn[1024];
	sprintf(fn, "%s/%s_H_U%.1f.txt", dir_output, method, pm->U);
	FILE *f = fopen(fn, "w");

	int i, j, is_sym=1;

	if(!pm->verbose) freopen("/dev/null", "w", stdout);

	printf("----------------------- Hamiltonian -----------------------\n");
	printf("%8s%8s%16s%16s\n", "row", "col", "val_real", "val_imag");
	for(i=0; i<H->nnz; i++) {
		for(j=0; j<H->nnz; j++) {
			if(((H->row[i] == H->col[j]) && (H->row[j] == H->col[i])) && (fabs(H->val[i] - H->val[j]) > 1e-6)) {
				is_sym = 0;
				break;
			}
		}
		printf("%8d%8d%16f%16f\n", H->row[i], H->col[i], creal(H->val[i]), cimag(H->val[i]));
		fprintf(f, "%8d%8d%16f%16f\n", H->row[i], H->col[i], creal(H->val[i]), cimag(H->val[i]));
		if(!is_sym) {
			printf("\nBroken symmetry:\n");
			printf("%8d%8d%16f%16f\n\n", H->row[j], H->col[j], creal(H->val[i]), cimag(H->val[j]));
			printf("-------------------------- ERROR --------------------------\n");
			break;
		}
	}
	printf("-----------------------------------------------------------\n");
	fclose(f);

	/*
	if(is_sym) {
		printf("------- check coo2csc -------\n");
		printf("%8s%8s%8s\n", "col_csc", "", "col_coo");
		for(i=0; i<pm->Nb+1; i++) printf("%8d%8s%8d\n", H->col_csc[i], "", i);
		printf("-----------------------------\n\n");
	}
	*/

	if(!pm->verbose) freopen("/dev/tty", "w", stdout);

	printf("%s: %s\n\n", __func__, fn);
}

