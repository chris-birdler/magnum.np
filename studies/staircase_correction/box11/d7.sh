#!/bin/bash
# box 11, PLAN 14b D7 / PLAN 14 "Speed measurement": prof_step.py on the production path (single precision), one run on
# the GPU at a time, old and new in turn; 30 sample intervals of warm-up, 96 measured; 3 runs for each arm and
# configuration. Old: M_9e0 (main 9e0ea7c, correction off), B_c63 (c63e862, correction on). New: ABP (AB with the
# patch). Kernel cache: one folder for each tree and N; the first run of each folder (kernel build) is a separate run.
# --trace: t, tries, <Ms m_x> and accepted steps after each sample interval, for check G; null runs of the old code at the end.
export CUDA_DEVICE=0
P=/root/progress.txt
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
W=/root/d7; mkdir -p $W /root/kc
run() {  # NAME TREE CONF SKIP PLAIN [PROF_STEP_OPTIONS]
  n=$1; t=$2; c=$3; [ -f $W/$n.json ] && return
  echo "$(date -u +%T) start D7 $n" >> $P
  cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d7_$n PYTHONPATH=/root/$t python /root/s/prof_step.py --skip $4 --plain $5 --window 0 --trace $6 --kcache /root/kc/${t}_${NN[$c]} --json $W/$n.json -- ${CONF[$c]} --b_list 0.006 --out $W/run_$n > $W/$n.log 2>&1
  rc=$?; echo "$(date -u +%T) end D7 $n exit $rc" >> $P
}
# kernel builds (one for each tree and N)
for x in M_9e0:off108 ABP:off108 B_c63:q4_108 M_9e0:off156 ABP:off156 B_c63:q4_156; do
  t=${x%%:*}; c=${x##*:}; mkdir -p /root/kc/${t}_${NN[$c]}; run kb_${t}_${NN[$c]} $t $c 0 1
done
# measurement: old and new in turn
for r in 1 2 3; do
  for c in off108 off156 q4_108 q4_156 q8_108 q8_156; do
    case $c in off*) old=M_9e0 ;; q4_*) old=B_c63 ;; *) old= ;; esac
    [ -n "$old" ] && run ${c}_old_r$r $old $c 30 96
    run ${c}_new_r$r ABP $c 30 96
  done
done
# null runs of the old code (scale of G if the old runs have no scatter; PLAN 14b D7 amendment)
for c in off108 off156 q4_108 q4_156; do
  case $c in off*) old=M_9e0 ;; *) old=B_c63 ;; esac
  for nl in reverse rotate; do run ${c}_null_${nl} $old $c 30 96 "--null $nl"; done
done
echo "$(date -u +%T) DONE d7" >> $P
