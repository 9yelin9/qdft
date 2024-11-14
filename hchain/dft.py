#!/usr/bin/env python3

import os
import sys
import shutil
import argparse
import numpy as np
import matplotlib.pyplot as plt

from scipy.constants import physical_constants

from pymatgen.core import Lattice, Structure, Molecule
from pymatgen.io.vasp import Poscar, Potcar, PotcarSingle, Kpoints, Incar
from pymatgen.io.vasp.outputs import Outcar

parser = argparse.ArgumentParser()
parser.add_argument('N', type=int)
parser.add_argument('R', type=float)
parser.add_argument('-i', '--gen_input', const='false', nargs='?', metavar='remove_old=false/true')
parser.add_argument('-o', '--get_output', action='store_true')
args = parser.parse_args()

path_pmgrc = '%s/.pmgrc.yaml' % os.environ['HOME']
if not os.path.exists(path_pmgrc):
	with open(path_pmgrc, 'w') as f:
		f.write('PMG_VASP_PSP_DIR=/APP/enhpc/VASP/POT\n')
PotcarSingle.functional_dir['PBE_54'] = 'PAW_PBE_54'

class DFT:
	def __init__(self, N, R):
		self.method = 'dft'

		self.N = N
		self.R = R
		print(f'N = {self.N}\nR = {self.R}', end='\n\n')

		self.dir_output = f'output/N{self.N}/{self.method}_R{self.R:.1f}'
		os.makedirs(self.dir_output, exist_ok=True)

		fn_list = ['poscar', 'potcar', 'kpoints', 'incar', 'outcar']
		self.fn = {fn: f'{self.dir_output}/{fn.upper()}' for fn in fn_list}

	def gen_input(self, remove_old='false'):
		if remove_old == 'true': 
			print(f'Remove {self.dir_output} ... ', end='')
			shutil.rmtree(self.dir_output)
			os.makedirs(self.dir_output)
			print('Done', end='\n\n')

		a = self.R * self.N
		lattice = Lattice.from_parameters(a=a, b=a, c=a, alpha=90, beta=90, gamma=90, pbc=(True, False, False))
		structure = Structure(lattice, ['H']*self.N, [[i / self.N, 0.5, 0.5] for i in range(self.N)])

		Poscar(structure).write_file(self.fn['poscar'])
		Potcar(symbols=['H'], functional='PBE_54').write_file(self.fn['potcar'])
		#Kpoints.automatic_density(structure, 1000).write_file(self.fn['kpoints'])
		Kpoints.monkhorst_automatic([1, 1, 1]).write_file(self.fn['kpoints'])
		Incar({
			'GGA': 'CA',		  # LDA functional
			'ENCUT': 400,         # Cutoff energy
			'ISMEAR': 0,          # Gaussian smearing for insulators/semiconductors
			'SIGMA': 0.01,        # Smearing width
			'EDIFF': 1E-6,        # Convergence criteria for electronic steps
			'IBRION': -1,         # Atomic relaxation
			'ISIF': 0,            # Stress and relaxation
			'NSW': 0,             # Number of ionic steps
			'ISPIN': 2,			  # Spin-polarized calculation
			'MAGMOM': '8*1.0',	  # Initial magnetic moments
			'LWAVE': '.FALSE.',   # Do not write WAVECAR
			'LCHARG': '.FALSE.',  # Do not write CHGCAR
		}).write_file(self.fn['incar'])

		print(f'{self.gen_input.__name__}:')
		print(*[self.fn[fn] for fn in ['poscar', 'potcar', 'kpoints', 'incar']], sep='\n', end='\n\n')
	
	def get_output(self):
		outcar = Outcar(self.fn['outcar'])
		ev2ht = physical_constants['electron volt-hartree relationship'][0]

		print(f'Final energy: {outcar.final_energy:f} eV = {outcar.final_energy * ev2ht:f} Hartree')
		print()
	
if len(sys.argv) < 2:
	print(f'Usage: {sys.argv[0]} <N> <R>')
	sys.exit(1)

dft = DFT(args.N, args.R)
if args.gen_input: dft.gen_input(args.gen_input)
elif args.get_output: dft.get_output()
else: parser.print_help()
