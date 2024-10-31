#!/usr/bin/env python3

import os
import ctypes
import numpy as np
import pennylane as qml

class QSoft:
	def __init__(self, N, Ne, U):
		self.N = N
		self.Nx = 2 * self.N
		self.Ne = Ne
		self.Nb = self.N
		self.U = U
		#self.beta = beta

	def func(self):
		pass

qs = QSoft(8, 4, 0)
