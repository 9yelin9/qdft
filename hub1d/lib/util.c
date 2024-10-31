#include "hub1d.h"

int compare_basis(const void *key, const void *p) {
	int key_val = *(int*)key;
	int p_val = ((basis*)p)->val;
	return (key_val > p_val) - (key_val < p_val);
}

void dec2bin(int len, int dec, char *bin) {
	int i;
	for(i=0; i<len; i++) bin[(len-1)-i] = (dec >> i) & 1 ? '1' : '0';
	bin[i] = '\0';
}

void gen_dir_output(params *pm, char *dir_output) {
	sprintf(dir_output, "output/N%d_Ne%d", pm->N, pm->Ne);
	mkdir(dir_output, 0755);
}

void save_basis(params *pm, basis *b, char *method) {
	char dir_output[1024], fn[1024];
	gen_dir_output(pm, dir_output);
	sprintf(fn, "%s/%s_basis.txt", dir_output, method);

	int i;
	char buf[pm->Nx+1];
	FILE *f = fopen(fn, "w");

	for(i=0; i<pm->Nb; i++) {
		dec2bin(pm->Nx, b[i].val, buf);
		fprintf(f, "%8d%8d%22s\n", b[i].idx, b[i].val, buf);
	}
	fclose(f);

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

void save_hamiltonian(params *pm, hamiltonian *h, char *method) {
	char dir_output[1024], fn[1024];
	gen_dir_output(pm, dir_output);
	sprintf(fn, "%s/%s_hamiltonian_U%.1f_e%f.txt", dir_output, method, pm->U, h->e_grd);

	int i;
	FILE *f = fopen(fn, "w");

	for(i=0; i<h->nnz; i++) fprintf(f, "%8d%8d%16f%16f\n", h->row[i], h->col[i], creal(h->val[i]), cimag(h->val[i]));
	fclose(f);

	printf("%s: %s\n\n", __func__, fn);
}

