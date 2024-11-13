#!/usr/bin/env python3

import os
import sys
import shutil
import argparse
import numpy as np

from scipy.constants import physical_constants

from pymatgen.core import Lattice, Structure, Molecule
from pymatgen.io.vasp import Poscar, Potcar, PotcarSingle, Kpoints, Incar
from pymatgen.io.vasp.outputs import Outcar

parser = argparse.ArgumentParser()
parser.add_argument('N', type=int)
parser.add_argument('R', type=int)
parser.add_argument('-i', '--gen_input', const='false', nargs='?', metavar='remove_old=false/true')
parser.add_argument('-o', '--show_output', action='store_true')
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

		self.dir_output = f'output/N{self.N}_R{self.R:.1f}/{self.method}/'
		os.makedirs(self.dir_output, exist_ok=True)

	def gen_input(self, remove_old='false'):
		if remove_old == 'true': 
			print(f'Remove {self.dir_output}... ', end='')
			shutil.rmtree(self.dir_output)
			os.makedirs(self.dir_output)
			print('Done', end='\n\n')

		lattice = Lattice.from_parameters(a=self.R * self.N, b=10, c=10, alpha=90, beta=90, gamma=90)
		structure = Structure(lattice, ['H']*self.N, [[i / self.N, 0.5, 0.5] for i in range(self.N)])

		fn = {fn: f'{self.dir_output}/{fn.upper()}' for fn in ['poscar', 'potcar', 'kpoints', 'incar']}

		Poscar(structure).write_file(fn['poscar'])
		Potcar(symbols=['H'], functional='PBE_54').write_file(fn['potcar'])
		Kpoints.automatic_density(structure, 1000).write_file(fn['kpoints'])
		Incar({
			'ENCUT': 400,         # Cutoff energy
			'ISMEAR': 0,          # Gaussian smearing for insulators/semiconductors
			'SIGMA': 0.01,        # Smearing width
			'EDIFF': 1E-6,        # Convergence criteria for electronic steps
			'IBRION': 2,          # Atomic relaxation
			'ISIF': 3,            # Stress and relaxation
			'NSW': 50,            # Number of ionic steps
			'LWAVE': '.FALSE.',   # Do not write WAVECAR
			'LCHARG': '.FALSE.',  # Do not write CHGCAR
		}).write_file(fn['incar'])

		print(f'{self.gen_input.__name__}:')
		print(*fn.values(), sep='\n')
		
	def show_output(self):
		outcar = Outcar(f'{self.dir_output}/OUTCAR')
		ev2ht = physical_constants['electron volt-hartree relationship'][0]

		print(f'Final energy: {outcar.final_energy} eV = {outcar.final_energy * ev2ht} Hartree')
		print(f'Minimum distance: {outcar.distance}')
	
if len(sys.argv) < 2:
	print(f'Usage: {sys.argv[0]} <N> <R>')
	sys.exit(1)

dft = DFT(args.N, args.R)
if args.gen_input: dft.gen_input(remove_old=args.gen_input)
elif args.show_output: dft.show_output()
else: parser.print_help()
