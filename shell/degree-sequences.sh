#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=50000

#SBATCH --job-name=degree-sequences
#SBATCH --output=degree-sequences.out
#SBATCH --mail-user=neilmacl@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc
module load openmpi
module load r

Rscript ../sims/degree-sequences.R --ntrials=100
Rscript ../sims/degree-sequences.R --ntrials=100 --optimize-weights

Rscript ../sims/degree-sequences.R --ntrials=100 --network=ba
Rscript ../sims/degree-sequences.R --ntrials=100 --network=ba --optimize-weights
