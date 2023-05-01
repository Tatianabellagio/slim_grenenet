#!/bin/bash
# Basic for loop
optimas=(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33)
for optima_value in "${optimas[@]}"
do
slim -d optima=$optima_value -d replicates=12 -d gen_to_run=3 arabidopsis_evolve.slim

done
echo simulation done