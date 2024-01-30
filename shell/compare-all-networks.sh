#!/bin/bash

declare -a dynamics=("dw" "SIS" "genereg" "mutualistic")

for i in ${dynamics[@]}; do
    sbatch compare-networks.sh $i
done

