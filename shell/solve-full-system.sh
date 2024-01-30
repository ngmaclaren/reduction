#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=26 # this is the number I care about for parallel
#SBATCH --mem=50000

#SBATCH --job-name=solve-full-systems
#SBATCH --output=solve-full-systems.out
#SBATCH --mail-user=neilmacl@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc
module load openmpi
module load r

network="${1}"

srun Rscript ../sims/solve-full-system.R ${network}
