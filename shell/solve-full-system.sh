#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
### These next two are backed off slightly from the per-node max
#SBATCH --cpus-per-task=51 # this is the number I care about for parallel
#SBATCH --mem=500000

## #SBATCH --job-name="solve-full-systems-newdolphintest"
## #SBATCH --output=solve-full-systems-newdolphintest.out
#SBATCH --mail-user=neilmacl@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc/11.2.0 openmpi/4.1.1 r/4.2.0

network="${1}"

## srun R CMD BATCH solve-full-systems.R dolphintest.Rout dolphin
srun Rscript ../sims/solve-full-system.R ${network}
