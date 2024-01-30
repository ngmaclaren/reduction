#!/bin/bash

# establish environment
module load gcc
module load openmpi
module load r

declare -a networks=("dolphin" "celegans" "proximity" "euroroad" "email" "er" "gkk" "ba" "hk" "lfr")
declare -a dynamics=("dw" "SIS" "genereg" "mutualistic")

cd ../analysis/

# dolphin-demo
Rscript dolphin-demo-plots.R --save-plots
Rscript dolphin-demo-plots.R --save-plots --use-weights

# degree-sequences
Rscript degree-sequences-plots.R --save-plots
Rscript degree-sequences-plots.R --save-plots --use-weights
Rscript degree-sequences-plots.R --network=ba --save-plots
Rscript degree-sequences-plots.R --network=ba --save-plots --use-weights

# knnfig
for j in ${dynamics[@]}; do
    Rscript knnfig-plots.R --dynamics=$j --save-plots
done

# compare-networks
for i in ${dynamics[@]}; do
    Rscript compare-networks-plots.R --dynamics=$i --save-plots
    Rscript compare-networks-plots.R --dynamics=$i --save-plots --use-weights
done

# transfer-learning
# for i in "${networks[@]}"; do
for j in ${dynamics[@]}; do
    Rscript transfer-learning-plots.R --network=dolphin --dynamics=$j
done
# done

# degrees and weights
Rscript analyze-weights.R --network=dolphin
Rscript analyze-weights.R --network=ba
