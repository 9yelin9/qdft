#!/usr/bin/env python3

import sys
import argparse
import numpy as np
from lib.dft import DFT

np.set_printoptions(suppress=True)

parser = argparse.ArgumentParser()
parser.add_argument('method', type=str)
parser.add_argument('N', type=int)
parser.add_argument('R', type=float)
parser.add_argument('--keep_old', action='store_true')
args = parser.parse_args()

if not args.method in ['vasp', 'espresso', 'wannier', 'qdft', 'qdftw']:
	print(f'ERROR: wrong method ({args.method})')
	sys.exit(1)

dft = DFT(args.method, args.N, args.R, keep_old=args.keep_old)
if   args.method == 'espresso': dft.run_espresso()
elif args.method == 'qdft':     dft.run_espresso(diagonalization='qdft')
elif args.method == 'wannier':  dft.run_wannier()
elif args.method == 'qdftw':    dft.run_wannier(use_vqe=True)
else: dft.run_vasp()
