import io
import sys
import h5py
import numpy as np
import pennylane as qml

input = io.BytesIO(sys.stdin.buffer.read())
with h5py.File(input, 'r') as f:
	M = f['M'][()]
	n_layer = f['n_layer'][()]
	use_gpu = f['use_gpu'][()]
	theta = f['theta'][:]
	basis = f['basis'][:]
	hamiltonian = f['hamiltonian'][:]
	projector = f['projector'][:]

wires = list(range(M))
hamiltonian = qml.Hamiltonian(hamiltonian.real.ravel(), [qml.Hermitian(p, wires=wires) for p in projector])
dev = qml.device('lightning.gpu' if use_gpu else 'default.qubit', wires=M)
print(qml.devices)

@qml.qnode(dev)
def circuit(theta, basis, hamiltonian):
	qml.BasisState(basis, wires=wires)
	for m in wires: qml.RY(theta[m], wires=m)
	for n in range(n_layer):
		for m in wires[:-1]: qml.CNOT(wires=[m, m+1])
		for m in wires:      qml.RY(theta[M * (n+1) + m], wires=m)
	return qml.expval(hamiltonian), qml.probs(wires=wires)

energy, occupation = [], []
for basis_k in basis:
	energy_k, occupation_k = circuit(theta, basis_k, hamiltonian)
	energy.append(energy_k)
	occupation.append(occupation_k)

output = io.BytesIO()
with h5py.File(output, 'w') as f:
	f.create_dataset('energy', data=energy)
	f.create_dataset('occupation', data=occupation)
sys.stdout.buffer.write(output.getvalue())
