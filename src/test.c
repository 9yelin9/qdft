#include "qdft.h"

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("Usage: %s <Ne>\n", argv[0]);
		exit(1);
	}

	environment env = {
		.Ne = atoi(argv[1]),
		.N = 2 * env.Ne,
		.M = 1 << env.N,
		.U = 0,
	};
	printf("Ne=%d\tN=%d\tM=%d\tU=%.1f\n\n", env.Ne, env.N, env.M, env.U);

	hamiltonian H = {
		.nnz = env.M * env.M,
		.row = (int*)malloc(sizeof(int) * H.nnz),
		.col = (int*)malloc(sizeof(int) * H.nnz),
		.col_csc = (int*)malloc(sizeof(int) * (env.M + 1)),
		.val = (double_complex*)malloc(sizeof(double_complex) * H.nnz),
	};

	int i, j;
	H.nnz = 0;
	for(i=0; i<env.M; i++) {
		H.col_csc[i] = H.nnz;
		for(j=0; j<env.M; j++) {
			if(abs(j-i) < 2) {
				H.row[H.nnz] = j;
				H.col[H.nnz] = i;
				H.val[H.nnz] = j >= i ? j+1 : i+1;
				H.nnz++;	
			}
		}
	}
	H.col_csc[i] = H.nnz;
	print_H(&env, &H);
	/* H = 
	 * 1, 2, 0, 0,
	 * 2, 2, 3, 0,
	 * 0, 3, 3, 4,
	 * 0, 0, 4, 4
	 */

	printf("laeig: %f\n", laeig(&env, &H));
	printf("areig: %f\n", areig(&env, &H));

	return 0;
}
