#!/bin/bash

networks=(windsurfers macaques train_terrorists highschool drug residence_hall netsci_weighted proximity_weighted gap_junction_herm intl_trade)
networks=()

for network in ${networks[@]}; do
    jobname=fs-$network
    # echo $jobname
    sbatch --job-name=$jobname --output=./output/${jobname}.out request-full-system.sh $network
done
