#!/bin/bash

networks=(canton ecoli email_company email_uni faa flamingo jung-c physician_trust polblogs protein_binding usair yeast)

for network in ${networks[@]}; do
    jobname=fs-$network

    sbatch --job-name=$jobname --output=./output/${jobname}.out request-full-system.sh $network
done

