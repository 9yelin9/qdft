#!/usr/bin/env python3

import sys
import argparse
import numpy as np
from lib.dft import DFT

np.set_printoptions(suppress=True)

parser = argparse.ArgumentParser()
parser.add_argument('method', type=str, choices=['vasp', 'espresso', 'wannier', 'openmx', 'openmx-vqe'])
parser.add_argument('N', type=int)
parser.add_argument('R', type=float)
parser.add_argument('--keep_old', action='store_true')
args = parser.parse_args()

dft = DFT(args.method, args.N, args.R, np=4, keep_old=args.keep_old)
if   args.method == 'vasp':     dft.run_vasp()
elif args.method == 'espresso': dft.run_espresso()
elif args.method == 'wannier':  dft.run_wannier()
elif args.method == 'openmx':   dft.run_openmx()
elif args.method == 'openmx-vqe':
	dft.np = 1
	dft.run_openmx(eigensolver='VQE')
