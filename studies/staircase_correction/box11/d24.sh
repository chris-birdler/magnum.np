#!/bin/bash
# box 11, PLAN 16 (PROTOKOLL 3.46), after d22.sh: the final speed measurement (request of Chris 2026-10-26), as D7 and
# W5: prof_step.py --skip 30 --plain 96 --window 0 --trace on the production path (single precision), one run on the
# GPU at a time, the codes in turn; 3 runs for each code and configuration. Codes: the start code without the speed
# patches (M_9e0: main 9e0ea7c, correction off; B_c63: c63e862, correction on), ABP (released), AS5 (final). The start
# code runs in off108, off156, q4_108, q4_156 (it does not fit at q 8); ABP and AS5 also in q8_108 and q8_156. Then
# q 16 with AS5 only (information): N 108 3 runs, N 156 2 runs (WQ of d22.sh is the third), if the transfer state
# exists. The kernel caches of D7 and W5 are used again (one folder for each tree and N; AS5 a copy of AS3).
export CUDA_DEVICE=0
P=/root/progress.txt
while ! grep -q "DONE d22$" $P; do sleep 30; done
eval "$(grep '^ARGS_A=' /root/s/k1.sh)"; eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
A="${ARGS_A/--precision_cycles double /}"; A="${A/--init_from \/root\/states\/K1_2_sat\/init.pt/}"
B="${ARGS_B/--precision_cycles double /}"; B="${B/--init_from \/root\/states\/Bprep\/init.pt/}"
N156() { x="${1/--d_lex 212.0/--d_lex 300.0}"; x="${x/--phi 0.62/--phi 0.55}"; x="${x/--Leff_lex 30.0/--Leff_lex 18.0}"; echo "${x/--seed 2013/--seed 2043}"; }
declare -A CONF NN
CONF[off108]="$A --init_from /root/states/K1_2_sat/init.pt"; NN[off108]=108
CONF[off156]="$(N156 "$A") --init_from /root/states/K4_2_sat/init.pt"; NN[off156]=156
CONF[q4_108]="$B --init_from /root/states/Bprep/init.pt"; NN[q4_108]=108
CONF[q4_156]="$(N156 "$B") --init_from /root/states/T_q4_156/init.pt"; NN[q4_156]=156
CONF[q8_108]="${B/--staircase_q 4/--staircase_q 8} --init_from /root/states/T_q8_108/init.pt"; NN[q8_108]=108
CONF[q8_156]="$(N156 "${B/--staircase_q 4/--staircase_q 8}") --init_from /root/states/T_q8_156/init.pt"; NN[q8_156]=156
CONF[q16_108]="${B/--staircase_q 4/--staircase_q 16} --init_from /root/states/T16_q16_108/init.pt"; NN[q16_108]=108
CONF[q16_156]="$(N156 "${B/--staircase_q 4/--staircase_q 16}") --init_from /root/states/T17_q16_156/init.pt"; NN[q16_156]=156
W=/root/p16/speed; mkdir -p $W
run() {  # NAME TREE CONF
  n=$1; t=$2; c=$3; [ -f $W/$n.json ] && return
  echo "$(date -u +%T) start P16 speed $n" >> $P
  cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/p16s_$n PYTHONPATH=/root/$t python /root/s/prof_step.py --skip 30 --plain 96 --window 0 --trace --kcache /root/kc/${t}_${NN[$c]} --json $W/$n.json -- ${CONF[$c]} --b_list 0.006 --out $W/run_$n > $W/$n.log 2>&1
  rc=$?; echo "$(date -u +%T) end P16 speed $n exit $rc" >> $P
}
for r in 1 2 3; do
  for c in off108 off156 q4_108 q4_156 q8_108 q8_156; do
    case $c in off*) old=M_9e0 ;; q4_*) old=B_c63 ;; *) old= ;; esac
    [ -n "$old" ] && run ${c}_${old}_r$r $old $c
    run ${c}_ABP_r$r ABP $c
    run ${c}_AS5_r$r AS5 $c
  done
done
[ -f /root/states/T16_q16_108/init.pt ] && for r in 1 2 3; do run q16_108_AS5_r$r AS5 q16_108; done
[ -f /root/states/T17_q16_156/init.pt ] && for r in 1 2; do run q16_156_AS5_r$r AS5 q16_156; done
echo "$(date -u +%T) DONE d24" >> $P
