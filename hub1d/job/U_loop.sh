#!/bin/bash

N=8
Ne=4
verbose=""

make clean; make
for U in `seq 0 0.5 9`
do
	#./fci  $N $Ne $U $verbose
	./soft $N $Ne $U $verbose
	./qsoft.py $N $Ne $U
done
