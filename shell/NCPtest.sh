#!/bin/bash

NCPs=(20 30 40 50 60 70 80 90 100)

for NCP in ${NCPs[@]}; do
    time Rscript ../sims/simulate-full-system.R --network=BAtest150 --ncparam=$NCP
done
