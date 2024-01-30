#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=50000
#SBATCH --time=00:15:00

# #SBATCH --job-name=knnfig-road
# #SBATCH --output=knnfig-road.out
#SBATCH --mail-user=neilmacl@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc
module load openmpi
module load r

seed=$RANDOM
network=$1

# Default network for this simulation is the email network
declare -a dynamics=("dw" "SIS" "genereg" "mutualistic")

for i in "${dynamics[@]}"; do
    Rscript ../sims/knnfig.R --dynamics=$i --network=$network --random-seed=$seed
    # Rscript ../sims/knnfig.R --dynamics=$i --optimize-weights
done
