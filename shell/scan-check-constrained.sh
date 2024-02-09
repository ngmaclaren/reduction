#!/bin/bash

cd ../analysis

networks=("dolphin" "celegans" "proximity" "euroroad" "email" "er" "ba" "hk" "gkk" "lfr")
dynamics=("dw" "SIS" "mutualistic" "genereg")

for net in ${networks[@]}; do
    for dyn in ${dynamics[@]}; do
	Rscript check-constrained.R --network=$net --dynamics=$dyn
    done
done
