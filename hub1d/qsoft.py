#!/usr/bin/env python3

import os
import sys
import scipy
import numpy as np
import pennylane as qml
from ctypes import *

np.set_printoptions(suppress=True)

class c_double_complex(Structure):
	_fields_ = [('creal', c_double),
				('cimag', c_double)]

class params(Structure):
	_fields_ = [('N',    c_int),
				('Nx',   c_int),
				('Ne',   c_int),
				('Nb',   c_int),
				('U',    c_double),
				('beta', c_double)]

class basis(Structure):
	_fields_ = [('idx', c_int),
				('val', c_int)]

class hamiltonian(Structure):
	_fields_ = [('nnz',     c_int),
				('row',     POINTER(c_int)),
				('col',     POINTER(c_int)),
				('col_csc', POINTER(c_int)),
				('e_grd',   c_double),
				('val',     POINTER(c_double_complex))]

class qsoft:
	def __init__(self, N, Ne, U):
		self.method = 'qsoft'
		self.libsoft = cdll.LoadLibrary('lib/libsoft.so')

		self.N = N
		self.M = int(np.log2(self.N))
		self.Nb = self.N
		self.Ne = Ne
		self.Nocc = (self.Ne + 1) // 2
		self.U = U
		self.beta = self.gen_beta(self.U)
		print(f'N = {self.N}\nNe = {self.Ne}\nNb = {self.Nb}\nU = {self.U}\nbeta = {self.beta}', end='\n\n')

		self.dir_output = f'output/N{self.N}_Ne{self.Ne}'
		os.makedirs(self.dir_output, exist_ok=True)

		self.wires = list(range(self.M))
		self.dev = qml.device('default.qubit', wires=self.M)
		self.qnode = qml.QNode(self.circuit, self.dev)

		self.e_grd = 100 # init e_grd

	def gen_beta(self, U, verbose=0):
		gen_beta_c = self.libsoft.gen_beta
		gen_beta_c.argtypes = [c_double, c_int]
		gen_beta_c.restype = c_double
		return gen_beta_c(U, verbose)

	def gen_basis_soft(self, pm, bs):
		gen_basis_soft_c = self.libsoft.gen_basis_soft
		gen_basis_soft_c.argtypes = [POINTER(params), POINTER(basis)]
		gen_basis_soft_c(byref(pm), bs)

	def init_soft(self):
		pm = params(self.N, self.N, self.Ne, self.Nb, self.U, self.beta)
		pm.beta = self.gen_beta(pm.U)

		bs = (basis * pm.Nb)()
		self.gen_basis_soft(pm, bs)

		hs = hamiltonian(pm.Nb**2)
		hs.row = (c_int * hs.nnz)()
		hs.col = (c_int * hs.nnz)()
		hs.val = (c_double_complex * hs.nnz)()

		return pm, bs, hs

	def gen_basis_qsoft(self):
		return np.array([list(np.binary_repr(i, width=self.M)) for i in range(self.Nb)], dtype='int')

	def gen_projector(self, bq):
		X  = np.array([[0,   1], [1,  0]], dtype='complex')	
		Y  = np.array([[0, -1j], [1j, 0]], dtype='complex')
		iY = 1j * Y 

		projector = []
		for i, j in [(i, j) for i in bq for j in bq]:
			p_ij = np.eye(self.Nb, dtype='complex')

			for mu in self.wires:
				Il = np.eye(2 ** mu)
				Ir = np.eye(2 ** (self.M - mu - 1))

				sign = (i[mu] << 1) - 1
				if i[mu] ^ j[mu]: p = (X - sign * iY) / 2
				else:             p = (X - sign * iY) @ (X + sign * iY) / 4
				p_ij @= np.kron(np.kron(Il, p), Ir)
			p_ij = (p_ij + p_ij.T.conj()) / 2
			projector.append(qml.Hermitian(p_ij, wires=self.wires))

		return projector

	"""
	def gen_projector(self, bq):
		operator = [[lambda mu: qml.FermiA(mu) * qml.FermiC(mu), lambda mu: qml.FermiA(mu)],
					[lambda mu: qml.FermiC(mu), lambda mu: qml.FermiC(mu) * qml.FermiA(mu)]]

		projector = []
		for i, j in [(i, j) for i in bq for j in bq]:
			p_ij = qml.matrix(qml.fermi.jordan_wigner(np.prod([operator[i[mu]][j[mu]](mu) for mu in self.wires])))
			p_ij = np.abs(p_ij + p_ij.T.conj()) / 2
			projector.append(qml.Hermitian(p_ij, wires=self.wires))
			#print(i, j, '\n', p_ij)

		return projector
	"""

	def circuit(self, theta, bq, hq, n_layer):
		qml.BasisState(bq, wires=self.wires)
		#qml.templates.StronglyEntanglingLayers(theta.reshape((n_layer, self.M, 3)), wires=self.wires)
		for m in self.wires: qml.RY(theta[m], wires=m)
		for n in range(n_layer):
			for m in self.wires[:-1]: qml.CNOT(wires=[m, m+1])
			for m in self.wires:      qml.RY(theta[self.M * (n+1) + m], wires=m)
		return qml.expval(hq), qml.probs(wires=self.wires)

	def weight(self, k):
		return (1 + self.Nocc - k) / (self.Nocc * (self.Nocc + 1) / 2)

	def energy_weighted(self, theta, bq, hq, n_layer):
		return np.sum([self.qnode(theta, bq[k], hq, n_layer)[0] * self.weight(k) for k in range(self.Nocc)])

	def energy(self, theta, bq, hq, n_layer):
		return np.sum([self.qnode(theta, bq[i//2], hq, n_layer)[0] for i in range(self.Ne)])

	def occupation(self, theta, bq, hq, n_layer):
		return np.sum([self.qnode(theta, bq[i//2], hq, n_layer)[1] for i in range(self.Ne)], axis=0)

	def energy_hartree(self, occ):
		return self.U * (occ / 2) * (occ / 2)
	def energy_hartree_deriv(self, occ):
		return self.U * (occ / 2)
	def energy_bethe(self, occ, beta):
		return -(2 * beta / np.pi) * np.sin(np.pi * occ / beta)
	def energy_bethe_deriv(self, occ, beta):
		return -2 * np.cos(np.pi * occ / beta)
	def energy_xc(self, occ):
		return self.energy_bethe(occ, self.beta) - self.energy_bethe(occ, 2) - self.energy_hartree(occ)
	def energy_xc_deriv(self, occ):
		return self.energy_bethe_deriv(occ, self.beta) - self.energy_bethe_deriv(occ, 2) - self.energy_hartree_deriv(occ)

	def save_basis(self, bq):
		fn = f'{self.dir_output}/{self.method}_basis.txt'
		with open(fn, 'w') as f:
			for idx, val in enumerate(bq):
				val = ''.join(map(str, val))
				f.write('%8d%8d%22s\n' % (idx, int(val, 2), val.zfill(self.M)))
		print(f'{self.save_basis.__name__}: {fn}', end='\n\n')

	def save_hamiltonian(self, hq, e_grd):
		fn = f'{self.dir_output}/{self.method}_hamiltonian_U{self.U:.1f}_e{e_grd:f}.txt'
		with open(fn, 'w') as f:
			for idx, val in np.ndenumerate(hq):
				f.write('%8d%8d%16f%16f\n' % (idx[0], idx[1], val.real, val.imag))
		print(f'{self.save_hamiltonian.__name__}: {fn}', end='\n\n')
	
	def run_qsoft(self, n_layer=4, itr_max=100, occ_mix=0.1):	
		pm, bs, hs = self.init_soft()
		hs_sorted = np.zeros((self.Nb, self.Nb), dtype='complex')
		gen_hamiltonian_soft_c = self.libsoft.gen_hamiltonian_soft
		gen_hamiltonian_soft_c.argtypes = [POINTER(params), POINTER(basis), POINTER(hamiltonian), POINTER(c_double)]

		bq = self.gen_basis_qsoft()
		pq = self.gen_projector(bq)
		self.save_basis(bq)

		print(f'---------------------------------------- kohn-sham self-consistent field + vqe ({n_layer}-layer) ----------------------------------------');
		print('%8s%12s%s' % ('itr', 'e_grd', ''.join(['%10s%02d' % ('occ', i) for i in range(self.N)])))

		#theta = np.full((n_layer, self.M, 3), 0.5) # init theta
		theta = np.array([0.5 for _ in range(self.M * (n_layer+1))]) # init theta
		occ = np.array([0.5 for _ in range(self.N)]) # init occ
		e_grd = 100 # init e_grd
		print('%8d%12s%s' % (0, '-', ''.join(map(lambda x: f'{x:12f}', occ))))

		for itr in range(1, itr_max+1):
			e_grd_old = e_grd
			occ_old = occ

			hs_sorted.fill(0)
			gen_hamiltonian_soft_c(byref(pm), bs, byref(hs), (c_double * self.N)(*occ))
			for i in range(hs.nnz): hs_sorted[hs.row[i], hs.col[i]] = hs.val[i].creal + 1j * hs.val[i].cimag
			hq = qml.Hamiltonian(np.ravel(hs_sorted).real, pq)
			#print(hs_sorted.reshape((self.Nb, self.Nb)), '\n', qml.matrix(hq))

			res = scipy.optimize.minimize(self.energy_weighted, theta, args=(bq, hq, n_layer), method='L-BFGS-B', options={'disp': False})

			theta = res.x
			occ = self.occupation(theta, bq, hq, n_layer)
			e_grd = self.energy(theta, bq, hq, n_layer)
			e_grd += np.sum(self.energy_xc(occ) + self.energy_hartree(occ))
			e_grd -= np.sum((self.energy_xc_deriv(occ_old) + self.energy_hartree_deriv(occ_old)) * occ)
			print('%8d%12f%s' % (itr, e_grd, ''.join(map(lambda x: f'{x:12f}', occ))))

			occ_tot = np.sum(occ)
			if not np.abs(occ_tot - self.Ne) < 1e-6:
				print(f'\nERROR: occ_tot({occ_tot}) != Ne({self.Ne})')
				sys.exit(1)

			if np.abs(e_grd - e_grd_old) < 1e-6: break
			occ = (1 - occ_mix) * occ_old + occ_mix * occ
		print('---------------------------------------------------------------------------------------------------------------------------------', end='\n\n');
		self.save_hamiltonian(qml.matrix(hq), e_grd)

if len(sys.argv) < 2:
	print(f'Usage: {sys.argv[0]} <N> <Ne> <U> [n_layer=4]')
	sys.exit(1)

qs = qsoft(int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]))
qs.run_qsoft(n_layer=4 if len(sys.argv) == 4 else int(sys.argv[4]))
