#!/bin/bash

# module load gcc openmpi r

cd ../analysis

dynamics=("dw" "SIS" "genereg" "mutualistic")

for j in ${dynamics[@]}; do
    Rscript knnfig-plots.R --dynamics=$j --save-plots --network=hk --random-seed=16318
done
