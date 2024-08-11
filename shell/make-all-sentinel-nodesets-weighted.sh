#!/bin/bash

networks=(windsurfers macaques train_terrorists highschool drug residence_hall netsci_weighted proximity_weighted gap_junction_herm intl_trade)
dynamics=(doublewell SIS genereg mutualistic)
ntrials=100
optweights=(false true)

for network in ${networks[@]}; do
    for dynamic in ${dynamics[@]}; do
	for optweight in ${optweights[@]}; do
	    jobname=ns-${network}-${dynamic}-${optweight}

	    sbatch --job-name=$jobname --output=./output/${jobname}.out request-sentinel-nodesets.sh $network $dynamic $ntrials $optweight
	done
    done
done

	    
