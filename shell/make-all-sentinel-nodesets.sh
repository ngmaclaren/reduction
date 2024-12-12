#!/bin/bash

# dropping powergrid and internet_as here to avoid confusion (getting down to 20 networks).

# full list
networks=(dolphin celegans proximity euroroad email er gkk ba hk lfr drosophila reactome route_views spanish foldoc tree_of_life word_assoc enron marker_cafe prosper)
dynamics=(doublewell SIS genereg mutualistic)
ntrials=50 # going down to 50 to split each run in half. Need to adjust file names accordingly
# optweights=(false true)
optweights=(false)
whichhalf=(first last)

for network in ${networks[@]}; do
    for dynamic in ${dynamics[@]}; do
	for optweight in ${optweights[@]}; do
	    for half in ${whichhalf[@]}; do
		jobname=ns-${network}-${dynamic}-${optweight}-${half}

		if [[ $network == "marker_cafe" || $network == "prosper" ]]; then
		    time="48:00:00"
		else
		    time="24:00:00"
		fi

		sbatch --time=$time --job-name=$jobname --output=./output/${jobname}.out request-sentinel-nodesets.sh $network $dynamic $ntrials $optweight $half
	    done
	done
    done
done

	    
