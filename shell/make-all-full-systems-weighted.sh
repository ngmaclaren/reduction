#!/bin/bash

networks=(windsurfers train_terrorists netsci_weighted intl_trade gap_junction_herm gap_junction_male drug biblenouns NZcollab proximity_weighted)
# networks=(train_bombers foodweb_dry foodweb_wet windsurfers highschool macaques residence_hall flights celegans_neural unicodelang proximity_weighted)

for network in ${networks[@]}; do
    jobname=fs-$network
    # echo $jobname
    sbatch --job-name=$jobname --output=./output/${jobname}.out request-full-system.sh $network
done
