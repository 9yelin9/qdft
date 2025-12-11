#define QDFT_PATH "/home/yerin/qdft/hchain/lib/QDFT"
#define double_complex double _Complex

#include <stdio.h>
#include <Python.h>
#include <complex.h>

void vqe(int *npw, int *npwx, int *nvec, int *nvecx, int *npol, double_complex *psi, double_complex *hpsi) {
	int i, j, idx;
	PyObject *pPath, *pModule, *pFunc, *pArgs, *pValue;

	Py_Initialize();
	pPath = PySys_GetObject("path");
	PyList_Append(pPath, PyUnicode_FromString(QDFT_PATH));

	pModule = PyImport_ImportModule("vqe"); //PyErr_Print();
	pFunc = PyObject_GetAttrString(pModule, "vqe");

	pArgs = PyTuple_Pack(4,
			PyLong_FromLong(*npw),
			PyLong_FromLong(*npwx),
			PyLong_FromLong(*nvec),
			PyLong_FromLong(*nvecx));
	pValue = PyObject_CallObject(pFunc, pArgs);

	printf("vqe.c: %d\t%d\t%d\n", *npw, *nvec, *npol);
	for(i=0; i<(*npw)*(*npol); i++) {
		printf("%d:", i);
		for(j=0; j<*nvec; j++) {
			idx = *nvec * i + j;
			printf("%d: %f+%fi\n", j, creal(hpsi[idx]), cimag(hpsi[idx]));
		}
		printf("\n");
	}
	printf("\n");

	Py_DECREF(pPath);
	Py_DECREF(pModule);
	Py_DECREF(pFunc);
	Py_DECREF(pArgs);
	Py_DECREF(pValue);
	Py_Finalize();
}
