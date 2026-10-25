#!/bin/bash
# box 11, PLAN 14b D7, part 4 (PROTOKOLL 3.41): the double-precision reference of part 3 also for q4_156 (G passes
# there with a small margin, and interval 115 has a local excess). c63e862 and ABP with the arguments of q4_156 of D7
# and --precision_cycles double, --trace, a new kernel folder for each run.
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
B="${ARGS_B/--init_from \/root\/states\/Bprep\/init.pt/}"
N156() { x="${1/--d_lex 212.0/--d_lex 300.0}"; x="${x/--phi 0.62/--phi 0.55}"; x="${x/--Leff_lex 30.0/--Leff_lex 18.0}"; echo "${x/--seed 2013/--seed 2043}"; }
C="$(N156 "$B") --init_from /root/states/T_q4_156/init.pt"
W=/root/d7c; mkdir -p $W
run() {  # NAME TREE
  n=$1; t=$2; [ -f $W/$n.json ] && return
  mkdir -p /root/kc/c_d_$n
  echo "$(date -u +%T) start D7c4 $n" >> $P
  cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d7c_$n PYTHONPATH=/root/$t python /root/s/prof_step.py --skip 30 --plain 96 --window 0 --trace --kcache /root/kc/c_d_$n --json $W/$n.json -- $C --b_list 0.006 --out $W/run_$n > $W/$n.log 2>&1
  rc=$?; echo "$(date -u +%T) end D7c4 $n exit $rc" >> $P
}
run q4_156_c63_double B_c63; run q4_156_ABP_double ABP
echo "$(date -u +%T) DONE d7c4" >> $P
