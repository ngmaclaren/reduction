#!/bin/bash

cd ../analysis

dynamics=("dw" "SIS" "genereg" "mutualistic")

for i in ${dynamics[@]}; do
    Rscript compare-methods-plots.R --dynamics=$i
done
