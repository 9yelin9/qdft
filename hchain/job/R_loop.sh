#!/bin/bash

t0=$(date +%s.%N)
t0_string=$(date)

method=$1
N=8

if [[ "$0" == "$BASH_SOURCE" ]] && [[ "$2" != "bg" ]]; then
	file_job=$(basename "$0")
	file_log_prefix="log/${file_job%.*}_${method}_$(date +%Y%m%d)"

	file_log="${file_log_prefix}.log"
	cnt=1
	while [[ -e "$file_log" ]]
	do
		file_log="${file_log_prefix}_${cnt}.log"
		((cnt++))
	done
	file_log="$(pwd)/${file_log}"
	echo $file_log

	nohup "$0" "$method" bg "$@" > $file_log 2>&1 &
	echo "Run in background (PID: $!)"
	exit 0
fi

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
