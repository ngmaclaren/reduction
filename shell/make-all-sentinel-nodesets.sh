#!/bin/bash

networks=("dolphin" "celegans" "proximity" "euroroad" "email" "er" "gkk" "ba" "hk" "lfr")
dynamics=("dw" "mutualistic" "SIS" "genereg")
ntrials=100
# nnodes=5

for network in ${networks[@]}; do
    for dynamic in ${dynamics[@]}; do
	jobname=sn-$network-$dynamic
	sbatch --job-name=${jobname} --output=./out/${jobname}.out request-sentinel-nodesets.sh ${ntrials} ${network} ${dynamic} # ${nnodes} (as $1)
    done
done
