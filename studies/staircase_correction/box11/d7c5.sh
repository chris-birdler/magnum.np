#!/bin/bash
# box 11, PLAN 14b D7, part 5 (PROTOKOLL 3.41): the new group of part 4 (q4_156) has 1 realisation (ABP is bitwise
# repeatable). 2 runs of AB (7f755a0, CSR product, not repeatable) on q4_156 with the single-precision arguments of D7
# give 2 more realisations of the new code; the distance to the double reference of part 4 is evaluated as in part 4.
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
B="${ARGS_B/--precision_cycles double /}"; B="${B/--init_from \/root\/states\/Bprep\/init.pt/}"
N156() { x="${1/--d_lex 212.0/--d_lex 300.0}"; x="${x/--phi 0.62/--phi 0.55}"; x="${x/--Leff_lex 30.0/--Leff_lex 18.0}"; echo "${x/--seed 2013/--seed 2043}"; }
C="$(N156 "$B") --init_from /root/states/T_q4_156/init.pt"
W=/root/d7c; mkdir -p $W
[ -d /root/kc/c_ABP_156 ] || cp -a /root/kc/ABP_156 /root/kc/c_ABP_156
run() {  # NAME TREE
  n=$1; t=$2; [ -f $W/$n.json ] && return
  echo "$(date -u +%T) start D7c5 $n" >> $P
  cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d7c_$n PYTHONPATH=/root/$t python /root/s/prof_step.py --skip 30 --plain 96 --window 0 --trace --kcache /root/kc/c_ABP_156 --json $W/$n.json -- $C --b_list 0.006 --out $W/run_$n > $W/$n.log 2>&1
  rc=$?; echo "$(date -u +%T) end D7c5 $n exit $rc" >> $P
}
run q4_156_AB_r1 AB; run q4_156_AB_r2 AB
echo "$(date -u +%T) DONE d7c5" >> $P
