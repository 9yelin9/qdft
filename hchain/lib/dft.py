import os
import time
import shutil
import numpy as np
from ase import Atoms
from ase.calculators.vasp import Vasp
from ase.calculators.espresso import Espresso, EspressoProfile

class DFT:
	def __init__(self, method, N, R, keep_old):
		self.method = method
		self.N = N
		self.R = R
		self.dir_output = f'{os.getcwd()}/output/N{self.N}/{self.method}_R{self.R:.2f}'
		if not keep_old:
			if os.path.isdir(self.dir_output): shutil.rmtree(self.dir_output)
			os.makedirs(self.dir_output, exist_ok=True)

		self.atom = 'H'
		self.structure = Atoms(
			self.atom * self.N,
			positions=[[self.R * i, 0, 0] for i in range(self.N)],
			cell=[self.R * self.N, 15, 15],
			pbc=True,
		)

		print(f'method={self.method}\nN = {self.N}\nR = {self.R}\ndir_output = {self.dir_output}', end='\n\n')

	def run_dft_vasp(self):
		t0 = time.time()

		calc = Vasp(
			directory=self.dir_output,
			command='mpirun -n 4 vasp_std',
			encut=300,  # Cutoff energy
			ismear=0,   # Gaussian smearing
			sigma=0.20, # Smearing width
			ediff=1e-6, # Convergence criteria for electronic steps
			icharg=2,   # Initial charge density
			nelm=100,   # Maximum iteration
			nsw=0,		# Ionic movement
			lwave='.FALSE.',
			lcharg='.FALSE.',
		)
		self.structure.calc = calc
		self.structure.get_potential_energy()

		tm, ts = divmod(int(time.time() - t0), 60)
		print(f'{self.run_dft_vasp.__name__} Done: {tm}m {ts}s', end='\n\n')

	def run_dft_espresso(self, diagonalization='david', pseudo_dir='/home/yerin/espresso/pseudo_pre'):
		t0 = time.time()
		
		profile = EspressoProfile(
			command='mpirun -n 4 pw.x',
			pseudo_dir=pseudo_dir,
		)
		input_data = {
			'control': {
				'calculation': 'scf',
				'restart_mode': 'from_scratch',
				'pseudo_dir': pseudo_dir,
				'outdir': self.dir_output,
			},
			'system': {
				'nat': self.N,
				'ntyp': 1,
				'ibrav': 0,
				'ecutwfc': 40,
				'occupations': 'smearing',
				'degauss': 0.2,
			},
			'electrons': {
				'diagonalization': diagonalization,
				'mixing_beta': 0.5,
				'conv_thr': 1e-6,
			},
		}
		pseudopotentials = {
			'H': 'H_ONCV_PBE-1.0.oncvpsp.upf',
		}

		calc = Espresso(
			directory=self.dir_output,
			profile=profile,
			input_data=input_data,
			pseudopotentials=pseudopotentials,
		)
		self.structure.calc = calc
		self.structure.get_potential_energy()

		tm, ts = divmod(int(time.time() - t0), 60)
		print(f'{self.run_dft_espresso.__name__} Done: {tm}m {ts}s', end='\n\n')
