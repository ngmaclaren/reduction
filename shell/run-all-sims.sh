#!/bin/bash

declare -a dynamics=("dw" "SIS" "genereg" "mutualistic")

sbatch dolphin-demo.sh

sbatch degree-sequences.sh

sbatch knnfig.sh

for i in ${dynamics[@]}; do
    sbatch compare-networks.sh $i
done

sbatch transfer-learning.sh


