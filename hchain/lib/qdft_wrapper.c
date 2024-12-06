#include <Python.h>

double run_qdft(double a, double b) {
	PyObject *pModule, *pFunc, *pArgs, *pValue;
	double res=0;

	Py_Initialize();
	PyRun_SimpleString("import sys; sys.path.append('.')");

	pModule = PyImport_ImportModule("qdft");
	pFunc = PyObject_GetAttrString(pModule, "run_qdft");
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
