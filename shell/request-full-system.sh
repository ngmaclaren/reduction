#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=52
#SBATCH --mem=100000

# job-name and output go in the top file
#SBATCH --mail-user=neilmacl@buffalo.edu
#SBATCH --mail-type=ALL
# #SBATCH --partition=general-compute
# #SBATCH --qos=general-compute
# #SBATCH --cluster=ub-hpc

# One-time cluster request to finish largest sim
#SBATCH --cluster=faculty
#SBATCH --partition=rabideau
#SBATCH --account=rabideau
#SBATCH --qos=rabideau

module load gcc openmpi r

network=$1

if [ $# -gt 1 ]
then
    splitsims=$2
    dynamics=$3
    split50=$4
    which50=$5
fi

echo $network
if [ $# -eq 1 ]
then
    Rscript ../sims/simulate-full-system.R --network=$network
else
    echo $splitsims
    echo $dynamics
    echo $which50
    
    Rscript ../sims/simulate-full-system.R --network=$network --splitsims=$splitsims --dynamics=$dynamics --split50=$split50 --which50=$which50
fi
