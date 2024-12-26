#include <stdio.h>
#include <unistd.h>
#include <Python.h>

#define QDFT_PATH "/home/yerin/qdft/hchain/lib/QDFT"

double run_qdft_(double a, double b) {
	double res;
	PyObject *pPath, *pModule, *pFunc, *pArgs, *pValue;

	//PyErr_Print();
	Py_Initialize();

	pPath = PySys_GetObject("path");
	PyList_Append(pPath, PyUnicode_FromString(QDFT_PATH));
	PyErr_Print();

	pModule = PyImport_ImportModule("qdft");
	PyErr_Print();
	pFunc = PyObject_GetAttrString(pModule, "run_qdft");
	PyErr_Print();
	pArgs = PyTuple_Pack(2, PyFloat_FromDouble(a), PyFloat_FromDouble(b));
	pValue = PyObject_CallObject(pFunc, pArgs);

	res = PyFloat_AsDouble(pValue);

	Py_DECREF(pModule);
	Py_DECREF(pFunc);
	Py_DECREF(pArgs);
	Py_DECREF(pValue);
	Py_Finalize();

	printf("%f\t%f\t%f\n", a, b, res);

	return res;
}

/*
int main() {
	double a=3.0, b=5.0;
	run_qdft_(a, b);
	return 0;
}
*/
