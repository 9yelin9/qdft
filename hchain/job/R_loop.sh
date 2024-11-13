#!/bin/bash

np=8
N=8

cd output
#for R in `seq 0.5 0.25 3.0`
for R in '2.0'
do
	cd N${N}_R${R}/dft
	mpirun -np $np vasp_std
	cd ../..
done
cd ..
