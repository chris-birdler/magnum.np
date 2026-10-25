#!/bin/bash
# box 11, PLAN 15 (PROTOKOLL 3.44), after W2 (d12.sh): (1) W3, repeatability of AS3 in single precision (as D3b):
# step_dump --k1_redo all --k1_thr 0 0 --k1_reps 2, windows k0 8 and 144, one arm with the correction (B, refB.pt)
# and one without (A, refA.pt); (2) the check of Dref: AS3 in double precision in q4_108, q4_156, off108 and off156
# (as d7c3.sh, d7c4.sh and d10a.sh; a new kernel folder for each run).
export CUDA_DEVICE=0
P=/root/progress.txt
while ! grep -q "DONE d12$" $P; do sleep 30; done
eval "$(grep '^ARGS_A=' /root/s/k1.sh)"; eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
A32="${ARGS_A/--precision_cycles double /}"; B32="${ARGS_B/--precision_cycles double /}"
W=/root/w3; mkdir -p $W
w3() {  # GEO REF ARGS
  n=W3_$1; [ -f $W/$n.pt ] && return
  echo "$(date -u +%T) start D13 $n" >> $P
  cd /root/AS3/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d13_$n PYTHONPATH=/root/AS3 python /root/s/step_dump.py --k1_arm $2 --k1_out $W/$n.pt --k1_ck 8 144 --k1_reps 2 --k1_redo all --k1_thr 0 0 -- $3 --b_list 0.006 --out $W/run_$n > $W/$n.log 2>&1
  rc=$?; echo "$(date -u +%T) end D13 $n exit $rc" >> $P
  python /root/s/glitch_count.py --reps $W/$n.pt --term_rel 0 --json $W/${n}_count.json > $W/${n}_count.txt 2>&1
}
w3 B /root/k1/B/refB.pt "$B32"
w3 A /root/k1/A/refA.pt "$A32"
N156() { x="${1/--d_lex 212.0/--d_lex 300.0}"; x="${x/--phi 0.62/--phi 0.55}"; x="${x/--Leff_lex 30.0/--Leff_lex 18.0}"; echo "${x/--seed 2013/--seed 2043}"; }
B2="${ARGS_B/--init_from \/root\/states\/Bprep\/init.pt/}"; A2="${ARGS_A/--init_from \/root\/states\/K1_2_sat\/init.pt/}"
D=/root/d10
run() {  # NAME ARGS
  n=$1; [ -f $D/$n.json ] && return
  mkdir -p /root/kc/c_d_$n
  echo "$(date -u +%T) start D13 $n" >> $P
  cd /root/AS3/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d13_$n PYTHONPATH=/root/AS3 python /root/s/prof_step.py --skip 30 --plain 96 --window 0 --trace --kcache /root/kc/c_d_$n --json $D/$n.json -- $2 --b_list 0.006 --out $D/run_$n > $D/$n.log 2>&1
  rc=$?; echo "$(date -u +%T) end D13 $n exit $rc" >> $P
}
run q4_108_AS3_double "$B2 --init_from /root/states/Bprep/init.pt"
run q4_156_AS3_double "$(N156 "$B2") --init_from /root/states/T_q4_156/init.pt"
run off108_AS3_double "$A2 --init_from /root/states/K1_2_sat/init.pt"
run off156_AS3_double "$(N156 "$A2") --init_from /root/states/K4_2_sat/init.pt"
echo "$(date -u +%T) DONE d13" >> $P
