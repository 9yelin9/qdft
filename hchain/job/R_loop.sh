#!/bin/bash

if [[ "$0" == "$BASH_SOURCE" ]] && [[ "$1" != "bg" ]]; then
	file_job=$(basename "$0")
	file_log="log/${file_job%.*}_$(date +%Y%m%d).log"

	cnt=1
	while [[ -e "$file_log" ]]
	do
		file_log="log/${file_job%.*}_$(date +%Y%m%d)_${cnt}.log"
		((cnt++))
	done
	file_log="$(pwd)/${file_log}"
	echo $file_log

	nohup "$0" bg "$@" > $file_log 2>&1 &
	echo "Run in background (PID: $!)"
	exit 0
fi

t0=$(date +%s.%N)
t0_string=$(date)

np=4
N=8

for R in '1.00'
#for R in `seq 0.70 0.10 3.00`
do
	./dft.py $N $R -i "true"
	cd output/N${N}/dft_R${R}
	mpirun -np $np vasp_std
	cd -
	./dft.py $N $R -o
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
