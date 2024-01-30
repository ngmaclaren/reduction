#!/bin/bash

networks=("dolphin" "celegans" "proximity" "euroroad" "email" "er" "gkk" "ba" "hk" "lfr")
dynamics=("dw" "SIS" "mutualistic" "genereg")

for network in ${networks[@]}; do
    for dynamicA in ${dynamics[@]}; do
	for dynamicB in ${dynamics[@]}; do
	    if [ "$dynamicA" == "$dynamicB" ]; then
		continue
	    fi
	    
	    Rscript ../analysis/dynamics-comparison.R --network=${network} --dynamicsA=$dynamicA --dynamicsB=$dynamicB
	done
    done
done

	    
