import ctypes
import numpy as np
import pennylane as qml

def vqe(na, nev, a, lda, ev, q, ldq, nblk):
	"""
	na          Order of matrix a

	nev         Number of eigenvalues needed.
				The smallest nev eigenvalues/eigenvectors are calculated.

	a(lda,*)    Distributed matrix for which eigenvalues are to be computed.
				Distribution is like in Scalapack.
				The full matrix must be set (not only one half like in scalapack).
				Destroyed on exit (upper and lower half).

	lda         Leading dimension of a

	ev(na)      On output: eigenvalues of a, every processor gets the complete set

	q(ldq,*)    On output: Eigenvectors of a
				Distribution is like in Scalapack.
				Must be always dimensioned to the full size (corresponding to (na,na))
				even if only a part of the eigenvalues is needed.

	ldq         Leading dimension of q

	nblk        blocksize of cyclic distribution, must be the same in both directions!
	"""
	print("vqe.py:", na, nev, lda, ldq, nblk)
	print("a:\n", a)
	print("ev:\n", ev)
	print("q:\n", q)
