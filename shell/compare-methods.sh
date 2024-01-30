#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=50000
#SBATCH --time=01:00:00

# #SBATCH --job-name=compare-methods
# #SBATCH --output=compare-methods.out
#SBATCH --mail-user=neilmacl@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc openmpi r

dynamics=$1

Rscript ../sims/compare-methods.R --dynamics=$dynamics --ntrials=100

