#!/bin/bash

networks=(dolphin celegans proximity euroroad email)

for network in ${networks[@]}; do
    time Rscript ../sims/select-sentinel-nodesets.R --network=$network --dynamics=doublewell --ntrials=9
done
