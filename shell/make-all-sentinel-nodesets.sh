#!/bin/bash

networks=(dolphin celegans proximity euroroad email er gkk ba hk lfr)
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

	    
