#!/bin/bash
#$ -S /bin/bash
#$ -q all.q
#$ -j y
#$ -cwd
#$ -o log/

#method="openmx-vqe-shot"
method=$1
N=$2
R=$3
#N=4
#R=2.5
#vqe_seed=$1

t0=$(date +%s.%N)
t0_string=$(date)

#./run.py $method $N $R --vqe_seed $vqe_seed
./run.py $method $N $R

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
