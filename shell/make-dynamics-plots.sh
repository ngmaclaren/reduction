#!/bin/bash

for i in "dolphin" "celegans" "proximity" "euroroad" "email" "er" "gkk" "ba" "hk" "lfr"; do
    Rscript ../analysis/compare-dynamics-plots.R -s -w -g $i
done
