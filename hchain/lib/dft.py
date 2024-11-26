import os
import shutil
import numpy as np
from pymatgen.core import Lattice, Structure, Molecule
from pymatgen.io.vasp import Poscar, Potcar, PotcarSingle, Kpoints, Incar
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath

path_pmgrc = '%s/.pmgrc.yaml' % os.environ['HOME']
if not os.path.exists(path_pmgrc):
	with open(path_pmgrc, 'w') as f:
		f.write('PMG_VASP_PSP_DIR=/APP/enhpc/VASP/POT\n')
PotcarSingle.functional_dir['PBE_54'] = 'PAW_PBE_54'
PotcarSingle.functional_dir['LDA_54'] = 'PAW_LDA_54'

class DFT:
	def __init__(self, method, N, R):
		self.method = method
		self.N = N
		self.R = R
		print(f'method={self.method}\nN = {self.N}\nR = {self.R}', end='\n\n')

		self.dir_output = f'output/N{self.N}/{self.method}_R{self.R:.2f}'
		os.makedirs(self.dir_output, exist_ok=True)

		fn_list = ['poscar', 'potcar', 'kpoints', 'incar', 'outcar']
		self.fn = {fn: f'{self.dir_output}/{fn.upper()}' for fn in fn_list}

		self.incar = {
			'ENCUT': 300,  # Cutoff energy
			'ISMEAR': 0,   # Gaussian smearing
			'SIGMA': 0.20, # Smearing width
			'EDIFF': 1E-6, # Convergence criteria for electronic steps
			'ICHARG': 2,   # Initial charge density
			'ISPIN': 2,	   # Spin-polarized calculation
			'NUPDOWN': 0,  # Difference between spin up and down
			'NELM': 100,   # Maximum iteration
			'NSW': 0,	   # Ionic movement
			'LWAVE': '.FALSE.',
			'LCHARG': '.FALSE.',
		}
		if self.method == 'qdft':
			self.incar['NELM'] = 0
	
	def remove_old(self, dir_name):
		print(f'Remove {dir_name} ... ', end='')
		shutil.rmtree(dir_name)
		os.makedirs(dir_name)
		print('Done', end='\n\n')

	def gen_input(self, remove_old='true', atom='H'):
		if remove_old == 'true': self.remove_old(self.dir_output)

		lattice = Lattice.from_parameters(a=15, b=15, c=self.R, alpha=90, beta=90, gamma=90)
		structure = Structure(lattice, [atom], [[0, 0, 0]])
		supercell = structure.copy().make_supercell([1, 1, self.N])

		Poscar(supercell).write_file(f'{self.dir_output}/POSCAR')
		Potcar(symbols=[atom], functional='PBE_54').write_file(f'{self.dir_output}/POTCAR')
		Kpoints.gamma_automatic([1, 1, 1]).write_file(f'{self.dir_output}/KPOINTS')
		Incar(self.incar).write_file(f'{self.dir_output}/INCAR')

		print(f'{self.gen_input.__name__}:')
		print(*[f'{self.dir_output}/{fn}' for fn in ['POSCAR', 'POTCAR', 'KPOINTS', 'INCAR']], sep='\n', end='\n\n')

	def gen_band(self, remove_old='true'):
		dir_band = f'{self.dir_output}/band'
		os.makedirs(dir_band, exist_ok=True)
		if remove_old == 'true': self.remove_old(dir_band)
		
		for fn in ['POSCAR', 'POTCAR', 'CHGCAR', 'WAVECAR', 'vasprun.xml']:
			shutil.copyfile(f'{self.dir_output}/{fn}', f'{dir_band}/{fn}')
		structure = Poscar.from_file(f'{dir_band}/POSCAR').structure
		structure = SpacegroupAnalyzer(structure).get_primitive_standard_structure()
		Kpoints.automatic_linemode(divisions=10, ibz=HighSymmKpath(structure)).write_file(f'{dir_band}/KPOINTS')
		Incar(dict(self.incar, ICHARG=11)).write_file(f'{dir_band}/INCAR')

		print(f'{self.gen_band.__name__}:')
		print(*[f'{dir_band}/{fn}' for fn in ['POSCAR', 'POTCAR', 'CHGCAR', 'WAVECAR', 'vasprun.xml', 'KPOINTS', 'INCAR']], sep='\n', end='\n\n')
	
	def get_output(self):
		outcar = Outcar(self.fn['outcar'])

		print(f'Final energy: {outcar.final_energy:f} eV')
		print(*outcar.final_energy_contribs.items(), sep='\n')
		print()
