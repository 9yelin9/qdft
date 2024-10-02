#include "hub1d.h"

int compare_basis(const void *key, const void *p) {
	int key_val = *(int*)key;
	int p_val = ((basis*)p)->val;
	return (key_val > p_val) - (key_val < p_val);
}

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

void print_basis(params *pm, basis *b, char *dir_output, char *method, int verbose) {
	char fn[1024];
	sprintf(fn, "%s/%s_basis.txt", dir_output, method);
	FILE *f = fopen(fn, "w");

	int i;
	char val_bin[pm->N+1];

	if(!verbose) freopen("/dev/null", "w", stdout);

	printf("---------------- basis ----------------\n");
	printf("%8s%8s%20s\n", "idx", "val(10)", "val(2)");
	for(i=0; i<pm->Nb; i++) {
		dec2bin(pm->N, b[i].val, val_bin);
		printf("%8d%8d%20s\n", b[i].idx, b[i].val, val_bin);
		fprintf(f, "%8d%8d%20s\n", b[i].idx, b[i].val, val_bin);
	}
	printf("---------------------------------------\n");
	fclose(f);

	if(!verbose) freopen("/dev/tty", "w", stdout);

	printf("%s: %s\n\n", __func__, fn);
}

int check_hamiltonian_hermitian(params *pm, hamiltonian *h) {
	int i, j;

	for(i=0; i<h->nnz; i++) {
		for(j=0; j<h->nnz; j++) {
			if(((h->row[i] == h->col[j]) && (h->row[j] == h->col[i])) && (fabs(h->val[i] - h->val[j]) > 1e-6)) {
				printf("\nERROR: hamiltonian is not hermitian (h[%d] != h[%d])\n", i, j);
				return 1;
			}
		}
	}
	return 0;
}

int check_coo2csc(params *pm, hamiltonian *h) {
	int i, col_old=-1, cnt=0;

	for(i=0; i<h->nnz; i++) {
		if(h->col[i] != col_old) {
			//printf("%d\t%d\t%d\n", h->col[i], i, h->col_csc[cnt]);
			if(i != h->col_csc[cnt]) {
				printf("\nERROR: coo2csc fail (%d != h_col_csc[%d])\n", i, cnt);
				return 1;
			}
			else cnt++;
		}
		col_old = h->col[i];
	}
	return 0;
}

void print_hamiltonian(params *pm, hamiltonian *h, char *dir_output, char *method, int verbose) {
	char fn[1024];
	sprintf(fn, "%s/%s_hamiltonian_U%.1f_e%f.txt", dir_output, method, pm->U, h->e_grd);
	FILE *f = fopen(fn, "w");

	int i;

	if(!verbose) freopen("/dev/null", "w", stdout);

	printf("----------------------- hamiltonian -----------------------\n");
	printf("%8s%8s%16s%16s\n", "row", "col", "val_real", "val_imag");
	for(i=0; i<h->nnz; i++) {
		printf("%8d%8d%16f%16f\n", h->row[i], h->col[i], creal(h->val[i]), cimag(h->val[i]));
		fprintf(f, "%8d%8d%16f%16f\n", h->row[i], h->col[i], creal(h->val[i]), cimag(h->val[i]));
	}
	printf("-----------------------------------------------------------\n");
	fclose(f);

	if(!verbose) freopen("/dev/tty", "w", stdout);

	printf("%s: %s\n\n", __func__, fn);
}

