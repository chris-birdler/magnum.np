#!/bin/bash

echo "" > timings_setup.dat
echo "" > timings_eval.dat

# run benchmarks
NN="10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250 260 270 280 290 300"
for N in $NN
do
    echo $N
    python run_magnumnp.py $N &> out_N${N}
    cat out_N${N} | grep t_setup | cut -f 2- >> timings_setup.dat
    cat out_N${N} | grep t_eval | cut -f 2- >> timings_eval.dat
done

# cleanup summary
cat out_N10 | grep names | cut -f 2- > tmp
cat timings_setup.dat >> tmp
cat tmp | column -t > timings_setup.dat

cat out_N10 | grep names | cut -f 2- > tmp
cat timings_eval.dat >> tmp
cat tmp | column -t > timings_eval.dat

rm tmp

