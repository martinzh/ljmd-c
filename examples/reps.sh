#!/bin/bash

for i in 1 2 3 4 5
do
    for j in 1
    do
        sbatch job.sh $j
    done
done
