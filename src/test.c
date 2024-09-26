#include "hub1d.h"

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("Usage: %s <Ni>\n", argv[0]);
		exit(1);
	}

	params pm = {
		.Ni = atoi(argv[1]),
		.N = 2 * pm.Ni,
		.Nb = 1 << pm.N,
		.verbose = 1,
	};
	printf("Ni=%d\nN=%d\nNb=%d\n\n", pm.Ni, pm.N, pm.Nb);
	char dir_output[128]; sprintf(dir_output, "output/test_Ni%d", pm.Ni);
	mkdir(dir_output, 0755);

	hamiltonian H = {
		.nnz = pm.Nb * pm.Nb,
		.row = (int*)malloc(sizeof(int) * H.nnz),
		.col = (int*)malloc(sizeof(int) * H.nnz),
		.col_csc = (int*)malloc(sizeof(int) * (pm.Nb + 1)),
		.val = (double_complex*)malloc(sizeof(double_complex) * H.nnz),
	};

	int i, j;
	H.nnz = 0;
	for(i=0; i<pm.Nb; i++) {
		H.col_csc[i] = H.nnz;
		for(j=0; j<pm.Nb; j++) {
			if(abs(j-i) < 2) {
				H.row[H.nnz] = j;
				H.col[H.nnz] = i;
				H.val[H.nnz] = j >= i ? j+1 : i+1;
				H.nnz++;	
			}
		}
	}
	H.col_csc[i] = H.nnz;
	print_H(&pm, &H, dir_output, "fci");
	/* H = 
	 * 1, 2, 0, 0,
	 * 2, 2, 3, 0,
	 * 0, 3, 3, 4,
	 * 0, 0, 4, 4
	 */

	printf("\n");
	printf("FCI(LAPACK): %f\n", laeig(&pm, &H));
	printf("FCI(ARPACK): %f\n", areig(&pm, &H));

	free(H.row);
	free(H.col);
	free(H.col_csc);
	free(H.val);

	return 0;
}
