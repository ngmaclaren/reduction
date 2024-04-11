#!/bin/bash

networks=(train_bombers foodweb_dry foodweb_wet windsurfers highschool macaques residence_hall flights celegans_neural unicodelang proximity_weighted)

for network in ${networks[@]}; do
    jobname=fs-$network

    sbatch --job-name=$jobname --output=./output/${jobname}.out request-full-system.sh $network
done
