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

if not args.method in ['vasp', 'espresso', 'qdft']:
	print(f'ERROR: wrong method ({args.method})')
	sys.exit(1)

dft = DFT(args.method, args.N, args.R, keep_old=args.keep_old)
dft.run_dft
