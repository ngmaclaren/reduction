#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=52
#SBATCH --mem=100000

# job-name and output go in the top file
#SBATCH --mail-user=neil.g.maclaren@gmail.com # neilmacl@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

# for timing purposes, constrain the hardware
#SBATCH --constraint=CPU-Gold-6448Y

module load gcc openmpi r

network=$1
dynamics=$2
ntrials=$3
optweights=$4
half=$5

Rscript ../sims/select-sentinel-nodesets.R --network=$network --dynamics=$dynamics --ntrials=$ntrials --optimize-weights=$optweights --whichhalf=$half
