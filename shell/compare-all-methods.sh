#!/bin/bash

dynamics=("dw" "SIS" "genereg" "mutualistic")

for i in ${dynamics[@]}; do
    jobname=compare-methods-$i
    sbatch --job-name=${jobname} --output=${jobname}.out compare-methods.sh $i
done
