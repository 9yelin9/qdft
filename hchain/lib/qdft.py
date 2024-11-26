import os
import numpy as np
import pennylane as qml
from vaspwfc import vaspwfc

class QDFT:
	def __init__(self, dir_output):
		self.dir_output = dir_output

		self.method = 'qdft'

	def func(self):
		pswfc = vaspwfc(f'{self.dir_output}/WAVECAR')
		
