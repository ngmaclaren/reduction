#!/bin/bash

networks=(dolphin celegans proximity euroroad email er gkk ba hk lfr)
large_networks=(drosophila powergrid reactome route_views spanish foldoc tree_of_life word_assoc internet_as enron)

for network in ${networks[@]}; do
    jobname=fs-$network

    sbatch --job-name=$jobname --output=./output/${jobname}.out request-full-system.sh $network
done

for network in ${large_networks[@]}; do
    jobname=fs-$network

    sbatch --time=72:00:00 --job-name=$jobname --output=./output/${jobname}.out request-full-system.sh $network
done
