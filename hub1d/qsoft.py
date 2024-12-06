#!/usr/bin/env python3

import io
import os
import sys
import h5py
import time
import scipy
import argparse
import paramiko
import subprocess
import numpy as np
import pennylane as qml
from ctypes import *

np.set_printoptions(suppress=True)

parser = argparse.ArgumentParser()
parser.add_argument('N', type=int)
parser.add_argument('Ne', type=int)
parser.add_argument('U', type=float)
parser.add_argument('-l', '--n_layer', type=int, default=4)
parser.add_argument('-gpu', '--use_gpu', action='store_true')
args = parser.parse_args()

class c_double_complex(Structure):
	_fields_ = [('creal', c_double),
				('cimag', c_double)]

class c_params(Structure):
	_fields_ = [('N',    c_int),
				('Nx',   c_int),
				('Ne',   c_int),
				('Nb',   c_int),
				('U',    c_double),
				('beta', c_double)]

class c_basis(Structure):
	_fields_ = [('idx', c_int),
				('val', c_int)]

class c_hamiltonian(Structure):
	_fields_ = [('nnz',     c_int),
				('row',     POINTER(c_int)),
				('col',     POINTER(c_int)),
				('col_csc', POINTER(c_int)),
				('e_grd',   c_double),
				('val',     POINTER(c_double_complex))]

class QSOFT:
	def __init__(self, N, Ne, U, n_layer, use_gpu):
		self.method = 'qsoft'
		self.libsoft = cdll.LoadLibrary('lib/libsoft.so')

		self.N = N
		self.M = int(np.log2(self.N))
		self.Nb = self.N
		self.Ne = Ne
		self.Nocc = (self.Ne + 1) // 2
		self.U = U
		self.beta = self.gen_beta(self.U)
		self.dir_output = f'output/N{self.N}_Ne{self.Ne}'
		os.makedirs(self.dir_output, exist_ok=True)

		self.dir_circuit = 'circuit'
		self.path_run_circuit = f'{self.dir_circuit}/run_circuit_{self.method}.py'
		self.wires = list(range(self.M))
		self.n_layer = n_layer
		self.use_gpu = use_gpu
		if self.use_gpu:
			self.ssh = self.init_ssh('172.26.123.11', 'cmat')
			self.run_circuit = self.run_circuit_gpu
		else:
			self.run_circuit = self.run_circuit_cpu

		print(f'N = {self.N}\nNe = {self.Ne}\nNb = {self.Nb}\nU = {self.U}\nbeta = {self.beta:f}\ndir_output = {self.dir_output}', end='\n\n')
		print(f'n_layer = {self.n_layer}\nuse_gpu = {self.use_gpu}', end='\n\n')

	def init_ssh(self, hostname, username):
		ssh = paramiko.SSHClient()
		ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
		ssh.connect(hostname, username=username)
		ssh.exec_command(f'mkdir -p {self.dir_circuit}')

		sftp = ssh.open_sftp()
		sftp.put(f'{self.path_run_circuit}', f'/home/{username}/{self.path_run_circuit}')
		sftp.close()

		return ssh

	def gen_beta(self, U):
		gen_beta_c = self.libsoft.gen_beta
		gen_beta_c.argtypes = [c_double]
		gen_beta_c.restype = c_double
		return gen_beta_c(U)

	def gen_basis_soft(self, params, basis):
		gen_basis_soft_c = self.libsoft.gen_basis_soft
		gen_basis_soft_c.argtypes = [POINTER(c_params), POINTER(c_basis)]
		gen_basis_soft_c(byref(params), basis)

	def init_soft(self):
		params = c_params(self.N, self.N, self.Ne, self.Nb, self.U, self.beta)
		params.beta = self.gen_beta(params.U)

		basis = (c_basis * params.Nb)()
		self.gen_basis_soft(params, basis)

		hamiltonian = c_hamiltonian(params.Nb**2)
		hamiltonian.row = (c_int * hamiltonian.nnz)()
		hamiltonian.col = (c_int * hamiltonian.nnz)()
		hamiltonian.val = (c_double_complex * hamiltonian.nnz)()

		gen_hamiltonian_soft_c = self.libsoft.gen_hamiltonian_soft
		gen_hamiltonian_soft_c.argtypes = [POINTER(c_params), POINTER(c_basis), POINTER(c_hamiltonian), POINTER(c_double)]

		return params, basis, hamiltonian, gen_hamiltonian_soft_c

	def gen_basis_qsoft(self):
		return np.array([list(np.binary_repr(i, width=self.M)) for i in range(self.Nb)], dtype='int')

	def gen_projector(self, basis):
		X  = np.array([[0,   1], [1,  0]], dtype='complex')	
		Y  = np.array([[0, -1j], [1j, 0]], dtype='complex')
		iY = 1j * Y 

		projector = []
		for i, j in [(i, j) for i in basis for j in basis]:
			p_ij = np.eye(self.Nb, dtype='complex')

			for mu in self.wires:
				Il = np.eye(2 ** mu)
				Ir = np.eye(2 ** (self.M - mu - 1))

				sign = (i[mu] << 1) - 1
				if i[mu] ^ j[mu]: p = (X - sign * iY) / 2
				else:             p = (X - sign * iY) @ (X + sign * iY) / 4
				p_ij @= np.kron(np.kron(Il, p), Ir)
			p_ij = (p_ij + p_ij.T.conj()) / 2
			#projector.append(qml.Hermitian(p_ij, wires=self.wires))
			projector.append(p_ij)

		return projector

	def gen_projector_malfunc(self, basis):
		operator = [[lambda mu: qml.FermiA(mu) * qml.FermiC(mu), lambda mu: qml.FermiA(mu)],
					[lambda mu: qml.FermiC(mu), lambda mu: qml.FermiC(mu) * qml.FermiA(mu)]]

		projector = []
		for i, j in [(i, j) for i in basis for j in basis]:
			p_ij = qml.matrix(qml.fermi.jordan_wigner(np.prod([operator[i[mu]][j[mu]](mu) for mu in self.wires])))
			p_ij = np.abs(p_ij + p_ij.T.conj()) / 2
			projector.append(qml.Hermitian(p_ij, wires=self.wires))

		return projector

	def save_basis(self, basis):
		fn = f'{self.dir_output}/{self.method}_basis.txt'
		with open(fn, 'w') as f:
			for idx, val in enumerate(basis):
				val = ''.join(map(str, val))
				f.write('%8d%8d%22s\n' % (idx, int(val, 2), val.zfill(self.M)))
		print(f'{self.save_basis.__name__}: {fn}', end='\n\n')

	def save_hamiltonian(self, hamiltonian, e_grd):
		fn = f'{self.dir_output}/{self.method}_hamiltonian_U{self.U:.1f}_e{e_grd:f}.txt'
		with open(fn, 'w') as f:
			for idx, val in np.ndenumerate(hamiltonian):
				f.write('%8d%8d%16f%16f\n' % (idx[0], idx[1], val.real, val.imag))
		print(f'{self.save_hamiltonian.__name__}: {fn}', end='\n\n')

	def init_circuit(self, theta, basis, hamiltonian, projector):
		input = io.BytesIO()
		with h5py.File(input, 'w') as f:
			f.create_dataset('M', data=self.M)
			f.create_dataset('n_layer', data=self.n_layer)
			f.create_dataset('use_gpu', data=self.use_gpu)
			f.create_dataset('theta', data=theta)
			f.create_dataset('basis', data=basis[:self.Nocc])
			f.create_dataset('hamiltonian', data=hamiltonian)
			f.create_dataset('projector', data=projector)
		input.seek(0)
		return input

	def run_circuit_gpu(self, theta, basis, hamiltonian, projector):
		input = self.init_circuit(theta, basis, hamiltonian, projector)
		stdin, stdout, stderr = self.ssh.exec_command(f'python3 {self.path_run_circuit}')
		stdin.write(input.getvalue())
		stdin.channel.shutdown_write()
		output = io.BytesIO(stdout.read())
		with h5py.File(output, 'r') as f:
			energy = f['energy'][:]
			occupation = f['occupation'][:]
		return energy, occupation

	def run_circuit_cpu(self, theta, basis, hamiltonian, projector):
		input = self.init_circuit(theta, basis, hamiltonian, projector)
		result = subprocess.run(['python3', f'{self.path_run_circuit}'], input=input.read(), capture_output=True)
		output = io.BytesIO(result.stdout)
		with h5py.File(output, 'r') as f:
			energy = f['energy'][:]
			occupation = f['occupation'][:]
		return energy, occupation

	def energy_weighted(self, theta, basis, hamiltonian, projector, weight):
		return (self.run_circuit(theta, basis, hamiltonian, projector)[0] * weight).sum()

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
	
	def run_qsoft(self, itr_max=100, occ_mix=0.1):	
		# SOFT
		pm, bs, hs, gen_hamiltonian_soft_c = self.init_soft()

		# QSOFT
		hq = np.zeros((self.Nb, self.Nb), dtype='complex')
		bq = self.gen_basis_qsoft()
		pq = self.gen_projector(bq)
		self.save_basis(bq)
		weight = [(1 + self.Nocc - k) / (self.Nocc * (self.Nocc + 1) / 2) for k in range(self.Nocc)]

		# init
		theta = np.full((self.M, self.n_layer+1), 0.5).ravel()
		occ = np.full(self.N, 0.5)
		e_grd = 100

		print(f'---------------------------------------- kohn-sham self-consistent field + vqe ({self.n_layer}-layer) ----------------------------------------');
		print('%8s%12s%s' % ('itr', 'e_grd', ''.join(['%10s%02d' % ('occ', i) for i in range(self.N)])))
		print('%8d%12s%s' % (0, '-', ''.join(map(lambda x: f'{x:12f}', occ))))

		t0 = time.time()
		for itr in range(1, itr_max+1):
			e_grd_old = e_grd
			occ_old = occ

			hq.fill(0)
			gen_hamiltonian_soft_c(byref(pm), bs, byref(hs), (c_double * self.N)(*occ))
			for i in range(hs.nnz): hq[hs.row[i], hs.col[i]] = hs.val[i].creal + 1j * hs.val[i].cimag

			res = scipy.optimize.minimize(self.energy_weighted, theta, args=(bq, hq, pq, weight), method='L-BFGS-B', options={'disp': False})

			theta = res.x
			e_grd_single, occ_single = self.run_circuit(theta, bq, hq, pq)
			occ = np.repeat(occ_single, 2, axis=0)[:self.Ne].sum(axis=0)
			e_grd = np.repeat(e_grd_single, 2)[:self.Ne].sum()
			e_grd += np.sum(self.energy_xc(occ) + self.energy_hartree(occ))
			e_grd -= np.sum((self.energy_xc_deriv(occ_old) + self.energy_hartree_deriv(occ_old)) * occ)
			print('%8d%12f%s' % (itr, e_grd, ''.join(map(lambda x: f'{x:12f}', occ))))

			occ_tot = np.sum(occ)
			if not np.abs(occ_tot - self.Ne) < 1e-6:
				print(f'\nERROR: occ_tot({occ_tot}) != Ne({self.Ne})')
				sys.exit(1)

			if np.abs(e_grd - e_grd_old) < 1e-6: break
			occ = (1 - occ_mix) * occ_old + occ_mix * occ
		tm, ts = divmod(int(time.time() - t0), 60)

		print('---------------------------------------------------------------------------------------------------------------------------------');
		print(f'{self.run_qsoft.__name__} Done: {tm}m {ts}s', end='\n\n')
		self.save_hamiltonian(hq, e_grd)
		if self.use_gpu: self.ssh.close()

qsoft = QSOFT(args.N, args.Ne, args.U, args.n_layer, args.use_gpu)
qsoft.run_qsoft()
