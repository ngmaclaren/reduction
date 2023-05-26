#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=50000

#SBATCH --job-name=compare-dynamics
#SBATCH --output=compare-dynamics.out
#SBATCH --mail-user=neilmacl@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc/11.2.0 openmpi/4.1.1 r/4.2.0

declare -a networks=("dolphin" "celegans" "proximity" "euroroad" "email" "er" "gkk" "ba" "hk" "lfr")

for i in "${networks[@]}"; do
    # Rscript ../sims/compare-dynamics.R -n 25 -g $i # 100
    Rscript ../sims/compare-dynamics.R -n 25 -g $i -w # 100
done
