#!/bin/bash

# networks=(dolphin celegans proximity euroroad email er gkk ba hk lfr drosophila powergrid reactome route_views spanish foldoc tree_of_life word_assoc internet_as enron)

# large_networks=(marker_cafe prosper)
# dynamics=(doublewell SIS genereg mutualistic)

# correcting mutualistic
large_networks=(prosper)
dynamics=(mutualistic)

### This is the simple version and should work for small- to medium-sized networks and/or dynamics that compute fairly quickly.
# for network in ${networks[@]}; do
#     jobname=fs-$network

#     sbatch --job-name=$jobname --output=./output/${jobname}.out request-full-system.sh $network
# done

for network in ${large_networks[@]}; do
    for dynamic in ${dynamics[@]}; do
	splitsims=TRUE
	split50=TRUE
	whichset=(last)

	for which50 in ${whichset[@]}; do
	    jobname=fs-$network-$dynamic-$which50
		
	    sbatch --time=90:00:00 --job-name=$jobname --output=./output/${jobname}.out request-full-system.sh $network $splitsims $dynamic $split50 $which50
	done
    done
done
