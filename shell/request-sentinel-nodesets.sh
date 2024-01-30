#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=100000
#SBATCH --time=00:15:00
## For making very large batches
# #SBATCH --time=04:00:00

#SBATCH --mail-user=neilmacl@buffalo.edu
#SBATCH --mail-type=NONE
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc openmpi r

# nnodes=$1
ntrials=$1
network=$2
dynamics=$3

Rscript ../sims/make-sentinel-nodesets.R --ntrials=${ntrials} --network=${network} --dynamics=${dynamics} # --nnodes=${nnodes} 
