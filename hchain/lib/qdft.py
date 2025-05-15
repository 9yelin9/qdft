import os
import time
#import torch
import scipy
import numpy as np
import pennylane as qml
#import multiprocessing as mp

class QDFT:
	def __init__(self, N):
		self.N  = N # num of spin-orbitals
		self.Ne = self.N # num of electrons
		#self.Nocc = self.Ne // 2
		self.Nocc = self.Ne # num of occupied orbitals

		self.M = int(np.log2(self.N)) # num of qubits
		self.Nl = 3 # num of layers
		self.wires = list(range(self.M))

		self.dev = qml.device('default.qubit', wires=self.M)
		self.qnode = qml.QNode(self.circuit, self.dev, diff_method='parameter-shift', interface='autograd')
		#self.qnode = qml.QNode(self.circuit, self.dev, diff_method='best', interface='torch')

		self.bases = np.array([list(np.binary_repr(i, width=self.M)) for i in range(self.N)], dtype='int')
		self.projector = self.gen_projector()
		self.weight = [(self.Nocc - k) / (self.Nocc * (self.Nocc + 1) / 2) for k in range(self.Nocc)]

		np.random.seed(42)
		self.theta = np.random.rand(self.M * (self.Nl+1))

		#self.theta = torch.tensor(self.theta, dtype=torch.float, requires_grad=True)
		#self.opt = torch.optim.Adam([self.theta], lr=0.01)
		#self.n_epoch = 1000
		#self.loss_tol = 1e-6
		#self.patience = 3
		print(f'<QDFT> Parameters\nN = {self.N}  Ne = {self.Ne}  Nocc = {self.Nocc}  M = {self.M}  Nl = {self.Nl}')

	def gen_projector(self):
		X  = np.array([[0,   1], [1,  0]], dtype='complex')	
		Y  = np.array([[0, -1j], [1j, 0]], dtype='complex')
		iY = 1j * Y 

		projector = []
		for i, j in [(i, j) for i in self.bases for j in self.bases]:
			p_ij = np.eye(self.N, dtype='complex')

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

	def circuit(self, theta, basis, hamiltonian):
		qml.BasisState(basis, wires=self.wires)
		for m in self.wires: qml.RY(theta[m], wires=m)
		for n in range(self.Nl):
			for m in self.wires[:-1]: qml.CNOT(wires=[m, m+1])
			for m in self.wires:      qml.RY(theta[self.M * (n+1) + m], wires=m)
		return qml.expval(hamiltonian), qml.state() #qml.probs(wires=self.wires)

	def run_circuit(self, theta, bases, hamiltonian):
		#with mp.Pool(len(bases)) as pool:
		#	res = pool.starmap(self.qnode, [(theta, basis, hamiltonian) for basis in bases])
		#return zip(*res)
		energy, state = zip(*[self.qnode(theta, basis, hamiltonian) for basis in bases])
		return energy, state

	def energy_weighted(self, theta, bases, hamiltonian):
		energy, _ = self.run_circuit(theta, bases, hamiltonian)
		return sum([e * w for e, w in zip(energy, self.weight)])

	def energy_weighted_grad(self, theta, bases, hamiltonian):
		return qml.grad(self.energy_weighted, argnum=0)(theta, bases, hamiltonian)

	def vqe(self, hamiltonian):
		#print(np.linalg.eig(np.reshape(hamiltonian, (self.N, self.N))))
		hamiltonian = qml.Hamiltonian(hamiltonian, self.projector)
		
		t0 = time.time()
		res = scipy.optimize.minimize(self.energy_weighted, self.theta,\
				args=(self.bases[:self.Nocc], hamiltonian),\
				method='L-BFGS-B',\
				options={'disp': True})	
		self.theta = res.x
		tm, ts = divmod(int(time.time() - t0), 60)

		energy, state = self.run_circuit(self.theta, self.bases[:self.Nocc], hamiltonian)
		energy = np.array(energy).tolist()
		state = np.array(state).real.ravel().tolist()
		#print(energy, state)

		"""
		loss_best, cnt = 100, 0
		for epoch in range(self.n_epoch):
			self.opt.zero_grad()
			loss = self.energy_weighted(self.theta, self.bases[:self.Nocc], hamiltonian)
			loss.backward()
			self.opt.step()

			if loss.item() < loss_best:
				loss_best = loss.item()
				cnt = 0
			elif abs(loss_best - loss.item()) < self.loss_tol:
				if cnt >= self.patience: break
				else: cnt += 1

		energy, state = self.run_circuit(self.theta, self.bases[:self.Nocc], hamiltonian)
		energy = [x.item() for x in energy]
		state = np.ravel([x.tolist() for x in state]).real.tolist()
		print(f'<VQE>  num_epoch = {epoch}  energy_weighted = {loss.item()}')
		"""

		print(f'{self.vqe.__name__} time elapsed: {tm}m {ts}s', end='\n\n')
		return energy, state

#QDFT(2).vqe([1, 2, 2, 3])


