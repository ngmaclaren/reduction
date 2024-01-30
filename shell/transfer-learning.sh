#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=50000

#SBATCH --job-name=transfer-learning
#SBATCH --output=transfer-learning.out
#SBATCH --mail-user=neilmacl@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc
module load openmpi
module load r

#declare -a networks=("dolphin" "celegans" "proximity" "euroroad" "email" "er" "gkk" "ba" "hk" "lfr")
declare -a networks=("dolphin")
declare -a dynamics=("dw" "SIS" "genereg" "mutualistic")

for i in "${networks[@]}"; do
    for j in "${dynamics[@]}"; do
	Rscript ../sims/transfer-learning.R --network=$i --dynamics=$j --ntrials=100
    done
done
