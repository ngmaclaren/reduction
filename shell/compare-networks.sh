#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=50000
#SBATCH --time=24:00:00

#SBATCH --job-name=compare-networks
#SBATCH --output=compare-networks.out
#SBATCH --mail-user=neilmacl@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc
module load openmpi
module load r

# declare -a dynamics=("dw" "SIS" "genereg" "mutualistic")
dynamics="${1}"

# for i in "${dynamics[@]}"; do
Rscript ../sims/compare-networks.R --dynamics=$dynamics --ntrials=100 
Rscript ../sims/compare-networks.R --dynamics=$dynamics --ntrials=100 --optimize-weights
# done
