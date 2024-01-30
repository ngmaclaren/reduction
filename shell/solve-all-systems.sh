#!/bin/bash

declare -a networks=("dolphin" "celegans" "proximity" "euroroad" "email" "er" "gkk" "ba" "hk" "lfr")

for i in "${networks[@]}"; do
    sbatch solve-full-system.sh $i
    sbatch solve-alt-system.sh $i
done
