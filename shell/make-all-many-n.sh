#!/bin/bash

networks=(dolphin celegans proximity euroroad email er gkk ba hk lfr)

for network in ${networks[@]}; do
    jobname=many-n-${network}

    sbatch --job-name=$jobname --output=./output/${jobname}.out request-many-n.sh $network
done
