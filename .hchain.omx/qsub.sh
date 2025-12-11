#!/bin/bash

for i in `seq 0 1 19`
do
	qsub job/openmx-vqe-shot.sh $i
done
