#!/bin/bash

# establish environment
## module load gcc/11.2.0 openmpi/4.1.1 r/4.2.0
declare -a networks=("dolphin" "celegans" "proximity" "euroroad" "email" "er" "gkk" "ba" "hk" "lfr")
declare -a dynamics=("dw" "SIS" "genereg" "mutualistic")

cd ../analysis/

# dolphin-demo
Rscript dolphin-demo-plots.R --save-plots
Rscript dolphin-demo-plots.R --save-plots --use-weights

# degree-sequences
Rscript degree-sequences-plots.R --save-plots
Rscript degree-sequences-plots.R --save-plots --use-weights

# knnfig
Rscript knnfig-plots.R --save-plots

# compare-networks
for i in "${dynamics[@]}"; do
    Rscript compare-networks-plots.R --dynamics=$i --save-plots
    Rscript compare-networks-plots.R --dynamics=$i --save-plots --use-weights
done

# transfer-learning
for i in "${networks[@]}"; do
    for j in "${dynamics[@]}"; do
	Rscript transfer-learning-plots.R --network=$i --dynamics=$j
    done
done
