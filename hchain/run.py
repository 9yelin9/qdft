#!/usr/bin/env python3

import sys
import argparse
import numpy as np
from lib.dft import DFT
from lib.qdft import QDFT

np.set_printoptions(suppress=True)

parser = argparse.ArgumentParser()
parser.add_argument('method', type=str)
parser.add_argument('N', type=int)
parser.add_argument('R', type=float)
parser.add_argument('-i', '--gen_input', const='true', nargs='?', metavar='remove_old=true/false')
parser.add_argument('-b', '--gen_band', const='true', nargs='?', metavar='remove_old=true/false')
parser.add_argument('-o', '--get_output', action='store_true')
args = parser.parse_args()

if not args.method in ['dft', 'qdft']:
	print(f'ERROR: wrong method ({args.method})')
	sys.exit(1)

dft = DFT(args.method, args.N, args.R)
qdft = QDFT(dft.dir_output)

if args.gen_input: dft.gen_input(args.gen_input)
elif args.gen_band: dft.gen_band(args.gen_band)
elif args.get_output: dft.get_output()
#elif args.get_output: qdft.func()
else: parser.print_help()
