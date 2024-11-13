#!/bin/bash

N=8
Ne=4

make clean; make
for U in `seq 0 0.5 9`
do
	#./fci  $N $Ne $U
	./soft $N $Ne $U
	./qsoft.py $N $Ne $U
done
