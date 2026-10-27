#!/bin/bash
# box 11, PLAN 14b D7: speed result (d7_eval.py) and check G (k_eval.py, D7 amendment) from the runs of d7.sh
W=/root/d7; R=/root/d7res; mkdir -p $R
python /root/s/d7_eval.py $W --json $R/speed.json > $R/speed.txt 2>&1; echo "speed exit $?"
for c in off108 off156 q4_108 q4_156; do
  python /root/s/k_eval.py --g_old $W/${c}_old_r1.json $W/${c}_old_r2.json $W/${c}_old_r3.json \
    --g_new $W/${c}_new_r1.json $W/${c}_new_r2.json $W/${c}_new_r3.json \
    --g_null $W/${c}_null_reverse.json $W/${c}_null_rotate.json --json $R/g_$c.json > $R/g_$c.txt 2>&1; echo "G $c exit $?"
done
for c in q8_108 q8_156; do
  python /root/s/k_eval.py --g_new $W/${c}_new_r1.json $W/${c}_new_r2.json $W/${c}_new_r3.json --json $R/g_$c.json > $R/g_$c.txt 2>&1
  echo "G $c (information) exit $?"
done
