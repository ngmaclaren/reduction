#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=50000
#SBATCH --time=70:00:00

#SBATCH --job-name=compare-networks
#SBATCH --output=compare-networks.out
#SBATCH --mail-user=neilmacl@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc/11.2.0 openmpi/4.1.1 r/4.2.0

declare -a dynamics=("dw" "SIS" "genereg" "mutualistic")

for i in "${dynamics[@]}"; do
    Rscript ../sims/compare-networks.R --dynamics=$i --ntrials=100 
    Rscript ../sims/compare-networks.R --dynamics=$i --ntrials=100 --optimize-weights
done
