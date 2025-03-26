import sys
import scipy
import numpy as np
import pennylane as qml

class QDFT:
	def __init__(self, N):
		self.N  = N
		self.Nl = 4
		self.Nb = self.N
		self.Nocc = self.N // 2
		self.M = int(np.log2(self.N))
		self.wires = list(range(self.M))
		self.dev = qml.device('default.qubit', wires=self.M)
		self.qnode = qml.QNode(self.circuit, self.dev)

		self.bases = np.array([list(np.binary_repr(i, width=self.M)) for i in range(self.N)], dtype='int')
		self.projector = self.gen_projector()
		self.theta = np.full((self.M, self.Nl+1), 0.5).ravel()
		self.weight = [(self.Nocc - k) / (self.Nocc * (self.Nocc + 1) / 2) for k in range(self.Nocc)]

	def gen_projector(self):
		X  = np.array([[0,   1], [1,  0]], dtype='complex')	
		Y  = np.array([[0, -1j], [1j, 0]], dtype='complex')
		iY = 1j * Y 

		projector = []
		for i, j in [(i, j) for i in self.bases for j in self.bases]:
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

	def circuit(self, theta, basis, hamiltonian):
		qml.BasisState(basis, wires=self.wires)
		for m in self.wires: qml.RY(theta[m], wires=m)
		for n in range(self.Nl):
			for m in self.wires[:-1]: qml.CNOT(wires=[m, m+1])
			for m in self.wires:      qml.RY(theta[self.M * (n+1) + m], wires=m)
		return qml.expval(hamiltonian), qml.state() #qml.probs(wires=self.wires)

	def run_circuit(self, theta, bases, hamiltonian):
		energy, state = zip(*[self.qnode(theta, basis, hamiltonian) for basis in bases])
		return np.array(energy), np.array(state)

	def energy_weighted(self, theta, bases, hamiltonian):
		energy, _ = self.run_circuit(theta, bases, hamiltonian)
		print(energy)
		return (energy * self.weight).sum()

	def vqe(self, hamiltonian):
		hamiltonian = qml.Hamiltonian(hamiltonian, self.projector)

		res = scipy.optimize.minimize(self.energy_weighted, self.theta,\
				args=(self.bases[:self.Nocc], hamiltonian),\
				method='L-BFGS-B',\
				options={'disp': True})	
		theta = res.x

		ev, state = self.run_circuit(self.theta, self.bases, hamiltonian)
		print(ev)

		return ev.real, state.real.ravel()

#QDFT(4).vqe([1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6, 4, 5, 6, 7])
