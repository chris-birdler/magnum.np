#!/bin/bash

echo "" > timings_setup.dat
echo "" > timings_eval.dat

# run benchmarks
for e in `seq 5 18`
do
    echo $e
    python run_magnumnp.py $e &> out_e${e}
done

# cleanup summary
cat out_e* | grep NN= | awk 'BEGIN {print "#NN n_eval t_total throughput"} {print $2, $4, $6, $8}' | sort -n  | column -t > timings.dat

# plot results
python plot.py
