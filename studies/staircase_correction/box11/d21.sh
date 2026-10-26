#!/bin/bash
# box 11, PLAN 16 E2 and L2 (PROTOKOLL 3.46), on the CPU of the box at the same time as d20.sh (GPU): the set-up of
# the staircase correction of AS3 against AS5 (fq_check.py: md5 of the centroid offsets, B, V, S, T, rows, _d, _v and
# the small-budget offsets), float64 as in production. T3 geometry (d/l_ex 96, phi 0.55, N 48) at q 2, 4, 8, 16; N 108
# at q 4, 8, 16; N 156 at q 8. 24 CPU threads, so that the GPU jobs keep CPU time.
export CUDA_DEVICE=-1 OMP_NUM_THREADS=24
P=/root/progress.txt
W=/root/p16/e2; mkdir -p $W
OLD=/root/AS3/magnumnp/field_terms/demag_staircase.py; NEW=/root/AS5/magnumnp/field_terms/demag_staircase.py
C108=/root/w5/run_q8_108_AS3_r1/config.json; C156=/root/w5/run_q8_156_AS3_r1/config.json
fq() {  # NAME CONFIG ARGS...
  n=$1; c=$2; shift 2; [ -f $W/$n.json ] && return
  echo "$(date -u +%T) start P16 E2 $n" >> $P; cd /root/s
  PYTHONPATH=/root/AS5:/root/AS5/studies/beta_bpeak python /root/s/fq_check.py --old $OLD --new $NEW --config $c "$@" --json $W/$n.json > $W/$n.log 2>&1
  echo "$(date -u +%T) end P16 E2 $n exit $?" >> $P
}
for q in 2 4 8 16; do fq t3_q$q $C108 --q $q --d_lex 96 --phi 0.55; done
fq n108_q4 $C108 --q 4
fq n108_q8 $C108 --q 8
fq n156_q8 $C156 --q 8
fq n108_q16 $C108 --q 16
echo "$(date -u +%T) DONE d21" >> $P
