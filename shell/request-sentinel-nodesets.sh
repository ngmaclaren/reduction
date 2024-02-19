#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=50000 # check the hardware here. I think it would be better as 100000

# job-name and output go in the top file
#SBATCH --mail-user=neilmacl@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc openmpi r

network=$1
dynamics=$2
ntrials=$3
optweights=$4

Rscript select-sentinel-nodesets.R --network=$network --dynamics=$dynamics --ntrials=$ntrials --optimize-weights=$optweights
