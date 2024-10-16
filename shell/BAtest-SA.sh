#!/bin/bash

networks=(BAtest50 BAtest75 BAtest100 BAtest125 BAtest150 BAtest175 BAtest200 BAtest225 BAtest250)

for network in ${networks[@]}; do
    time Rscript ../sims/select-sentinel-nodesets.R --network=$network --dynamics=doublewell --ntrials=9
done
