#!/bin/bash

networks=(dolphin celegans proximity euroroad email er gkk ba hk lfr)

for network in ${networks[@]}; do
    jobname=fs-$network

    sbatch --job-name=$jobname --output=./output/${jobname}.out request-full-system.sh $network
done
