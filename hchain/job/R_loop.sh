#!/bin/bash
#$ -S /bin/bash
#$ -q all.q
#$ -j y
#$ -cwd
#$ -o log/

method=$1
N=8

t0=$(date +%s.%N)
t0_string=$(date)

#./run.py $method $N 1
for R in `seq 0.50 0.10 3.00`
do
	./run.py $method $N $R
done	

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
