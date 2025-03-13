import os
import sys
import time
import shutil
import inspect
import subprocess
import numpy as np

from ase import Atoms
from ase.dft.kpoints import monkhorst_pack
from ase.calculators.vasp import Vasp
from ase.calculators.openmx import OpenMX
from ase.calculators.espresso import Espresso, EspressoProfile

class DFT:
	def __init__(self, method, N, R, np=4, keep_old=False):
		self.method = method
		self.N = N
		self.R = R
		self.np = np
		self.dir_output = f'{os.getcwd()}/output/N{self.N}/{self.method}_R{self.R:.2f}'
		if not keep_old:
			if os.path.isdir(self.dir_output): shutil.rmtree(self.dir_output)
			os.makedirs(self.dir_output, exist_ok=True)

		self.atom = 'H'
		self.atoms = Atoms(
			self.atom * self.N,
			positions=[[self.R * i, 0, 0] for i in range(self.N)],
			cell=[self.R * self.N, 15, 15],
			pbc=True,
		)

		print(f'method={self.method}\nN = {self.N}\nR = {self.R}\nnp = {self.np}\ndir_output = {self.dir_output}', end='\n\n')

	def run_vasp(self):
		self.atoms.calc = Vasp(
			directory=self.dir_output,
			command=f'mpirun -n {self.np} vasp_std',
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

		t0 = time.time()
		self.atoms.get_potential_energy()
		tm, ts = divmod(int(time.time() - t0), 60)

		print(f'{self.run_vasp.__name__} Done: {tm}m {ts}s', end='\n\n')

	def run_espresso(self, calculation='scf', kpts=None, diagonalization='david', pseudo_dir='/home/yerin/espresso/pseudo_pre'):
		args = locals(); del args['self']
		
		profile = EspressoProfile(
			command=f'mpirun -n {self.np} pw.x',
			pseudo_dir=pseudo_dir,
		)
		input_data = {
			'control': {
				'calculation': calculation,
				'restart_mode': 'from_scratch',
				'pseudo_dir': pseudo_dir,
				'outdir': self.dir_output,
			},
			'system': {
				'nat': self.N,
				'nosym': True,
				'noinv': True,
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

		self.atoms.calc = Espresso(
			directory=self.dir_output,
			profile=profile,
			input_data=input_data,
			pseudopotentials=pseudopotentials,
			kpts=kpts,
		)

		t0 = time.time()
		try: self.atoms.get_potential_energy()
		except:
			if calculation == 'nscf': pass
			else: raise
		tm, ts = divmod(int(time.time() - t0), 60)

		print(f'{self.run_espresso.__name__}{args}\nDone: {tm}m {ts}s', end='\n\n')

	def run_wannier(self, kpts=(4, 1, 1), basename='wannier'):
		self.run_espresso()
		self.run_espresso(calculation='nscf', kpts=kpts)

		with open(f'{self.dir_output}/{basename}.win', 'w') as f:
			f.write(f'num_bands = {self.atoms.calc.results['nbands']}\n')
			f.write(f'num_wann = {self.atoms.calc.results['nbands']}\n')
			f.write(f'write_hr = .true.\n\n')

			f.write(f'begin projections\n')
			f.write(f'\trandom\n')
			f.write(f'end projections\n\n')

			f.write(f'mp_grid = {' '.join(map(str, kpts))}\n\n')

			f.write(f'begin unit_cell_cart\n')
			f.write(f'\t{'\n\t'.join(' '.join([f'{i:10f}' for i in row]) for row in self.atoms.cell[:])}\n')
			f.write(f'end unit_cell_cart\n\n')

			f.write(f'begin atoms_cart\n')
			f.write(f'\t{'\n\t'.join(self.atom + ' '.join([f'{i:10f}' for i in row]) for row in self.atoms.get_positions())}\n')
			f.write(f'end atoms_cart\n\n')

			f.write(f'begin kpoints\n')
			f.write(f'\t{'\n\t'.join(' '.join([f'{i:10f}' for i in row]) for row in self.atoms.calc.results['ibz_kpoints'])}\n')
			f.write(f'end kpoints\n\n')

		with open(f'{self.dir_output}/{basename}.p2win', 'w') as f:
			f.write(f'&INPUTPP\n')
			f.write(f'\tprefix = pwscf\n')
			f.write(f'\toutdir = {self.dir_output}\n')
			f.write(f'/\n')

		t0 = time.time()
		subprocess.run(['wannier90.x', '-pp', basename], cwd=self.dir_output)
		with open(f'{self.dir_output}/{basename}.p2win', 'r') as fi, open(f'{self.dir_output}/{basename}.p2wout', 'w') as fo:
			subprocess.run(['mpirun', '-n', self.np, 'pw2wannier90.x'], cwd=self.dir_output, stdin=fi, stdout=fo)
		subprocess.run(['wannier90.x', basename], cwd=self.dir_output)
		tm, ts = divmod(int(time.time() - t0), 60)

		print(f'{self.run_wannier.__name__} Done: {tm}m {ts}s', end='\n\n')

	def run_openmx(self, eigensolver='cluster'):
		self.atoms.calc = OpenMX(
			label=self.dir_output + '/openmx',
			data_path='/home/yerin/openmx3.9/DFT_DATA19',
			command='openmx',
			mpi={'processes': self.np},
			definition_of_atomic_species=[['H', 'H6.0-s1', 'H_CA19']],
			xc='LDA',
			maxiter=1000,
			energy_cutoff=150.,
			smearing=300,
			kpts=(1, 1, 1),
			convergence=1e-6,
			spinpol=None,
			eigensolver=eigensolver,
			mixer='rmm-diis',
		)

		t0 = time.time()
		self.atoms.get_potential_energy()
		tm, ts = divmod(int(time.time() - t0), 60)

		print(f'{self.run_openmx.__name__} Done: {tm}m {ts}s', end='\n\n')
