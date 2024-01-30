#!/bin/bash

networks=(email email er hk gkk lfr)

counter=1
for net in ${networks[@]}; do
    sbatch --job-name=knnfig-${counter} --output=knnfig-${counter}.out knnfig.sh $net
    ((counter++))
done
