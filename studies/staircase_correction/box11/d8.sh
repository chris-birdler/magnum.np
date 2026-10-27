#!/bin/bash
# box 11, PROTOKOLL 3.42: (1) the old code M_9e0 with DemagFieldPBC in the drive cycles (--strayfield_cycles fd), the
# arguments of off108 and off156 of D7 otherwise, n 3 for each N, then 1 control run with Newell (the D7 arm off_old).
# (2) Profile of ABP in q4_108, q4_156, off108 and off156 (--window 10): GPU busy share and the split into the terms.
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
W=/root/d8; mkdir -p $W
run() {  # NAME TREE CONF PROF_STEP_OPTIONS [RUN_LOOPS_OPTIONS]
  n=$1; t=$2; c=$3; [ -f $W/$n.json ] && return
  k=/root/kc/c8_${t}_${NN[$c]}; [ -d $k ] || cp -a /root/kc/${t}_${NN[$c]} $k
  echo "$(date -u +%T) start D8 $n" >> $P
  cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d8_$n PYTHONPATH=/root/$t python /root/s/prof_step.py $4 --kcache $k --json $W/$n.json -- ${CONF[$c]} $5 --b_list 0.006 --out $W/run_$n > $W/$n.log 2>&1
  rc=$?; echo "$(date -u +%T) end D8 $n exit $rc" >> $P
}
M="--skip 30 --plain 96 --window 0"
for c in off108 off156; do
  for r in 1 2 3; do run ${c}_fd_r$r M_9e0 $c "$M" "--strayfield_cycles fd"; done
  run ${c}_newell_ctrl M_9e0 $c "$M"
done
for c in q4_108 q4_156 off108 off156; do run prof_ABP_$c ABP $c "--skip 30 --plain 20 --window 10"; done
echo "$(date -u +%T) DONE d8" >> $P
