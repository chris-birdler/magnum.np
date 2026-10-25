#!/bin/bash
# box 11, PLAN 14b D7, cause of the G fail in q4_108, part 3 (PROTOKOLL 3.41): double-precision reference. c63e862
# (old) and ABP on q4_108 with the arguments of D7 and --precision_cycles double (ARGS_B as in K1), --trace. Each run
# gets a new kernel folder (the kernel file of the single-precision runs is not used for double precision).
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
W=/root/d7c; mkdir -p $W
run() {  # NAME TREE
  n=$1; t=$2; [ -f $W/$n.json ] && return
  mkdir -p /root/kc/c_d_$n
  echo "$(date -u +%T) start D7c3 $n" >> $P
  cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d7c_$n PYTHONPATH=/root/$t python /root/s/prof_step.py --skip 30 --plain 96 --window 0 --trace --kcache /root/kc/c_d_$n --json $W/$n.json -- $ARGS_B --b_list 0.006 --out $W/run_$n > $W/$n.log 2>&1
  rc=$?; echo "$(date -u +%T) end D7c3 $n exit $rc" >> $P
}
run q4_108_c63_double B_c63; run q4_108_ABP_double ABP
echo "$(date -u +%T) DONE d7c3" >> $P
