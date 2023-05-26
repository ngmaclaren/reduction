#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=50000

#SBATCH --job-name=dolphin-demo
#SBATCH --output=dolphin-demo.out
#SBATCH --mail-user=neilmacl@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc/11.2.0 openmpi/4.1.1 r/4.2.0

Rscript ../sims/dolphin-demo.R
Rscript ../sims/dolphin-demo.R -w
