#!/bin/bash

t0=$(date +%s.%N)
t0_string=$(date)

np=4
N=$1
R=$2

./dft.py $N $R -i "true"
cd output/N${N}/dft_R${R}
mpirun -np $np vasp_std
cd -
./dft.py $N $R -o

t1=$(date +%s.%N)
t1_string=$(date)

t=$(echo "$t1 - $t0"|bc)
h=$(echo "($t/3600)"|bc)
m=$(echo "($t%3600)/60"|bc)
s=$(echo "($t%3600)%60"|bc)

echo ""
echo "# Time Start   : $t0_string"
echo "# Time End     : $t1_string"
echo "# Time Elapsed : ${h}h ${m}m ${s}s"
echo ""
