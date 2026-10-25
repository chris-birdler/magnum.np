#!/bin/bash
# box 11, PROTOKOLL 3.43: feasibility of q 16. For N 108 and N 156: transfer state as in d7t.sh with --staircase_q 16
# (tree B_det), then a speed run of ABP from it as the q 8 arms of D7 (n 1). nvidia-smi records the GPU memory.
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
B="${ARGS_B/--precision_cycles double /}"; B="${B/--init_from \/root\/states\/Bprep\/init.pt/}"
B16="${B/--staircase_q 4/--staircase_q 16}"
N156() { x="${1/--d_lex 212.0/--d_lex 300.0}"; x="${x/--phi 0.62/--phi 0.55}"; x="${x/--Leff_lex 30.0/--Leff_lex 18.0}"; echo "${x/--seed 2013/--seed 2043}"; }
W=/root/d9; mkdir -p $W
nvidia-smi --query-gpu=timestamp,memory.used --format=csv,noheader -l 5 > $W/smi.csv 2>&1 &
SMI=$!
xfer() {  # NAME ARGS SRC
  n=$1; a=$2; src=$3; S=/root/states/$n; [ -f $S/init.pt ] && return; mkdir -p $S
  echo "$(date -u +%T) start D9 transfer $n" >> $P
  cd /root/B_det/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d9t_$n PYTHONPATH=/root/B_det python /root/s/step_dump.py --at 1 --save $S/prep_m.pt -- $a --b_list 0.006 --init_transfer $src --out $S > $S.log 2>&1
  rc=$?; echo "$(date -u +%T) end D9 transfer $n exit $rc ($(ls $S/init.pt 2>/dev/null))" >> $P
}
run() {  # NAME N ARGS
  n=$1; [ -f $W/$n.json ] && return
  k=/root/kc/c9_ABP_$2; [ -d $k ] || cp -a /root/kc/ABP_$2 $k
  echo "$(date -u +%T) start D9 $n" >> $P
  cd /root/ABP/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d9_$n PYTHONPATH=/root/ABP python /root/s/prof_step.py --skip 30 --plain 96 --window 0 --trace --kcache $k --json $W/$n.json -- $3 --b_list 0.006 --out $W/run_$n > $W/$n.log 2>&1
  rc=$?; echo "$(date -u +%T) end D9 $n exit $rc" >> $P
}
xfer T_q16_108 "$B16" /root/states/K1_2_sat/init.pt
[ -f /root/states/T_q16_108/init.pt ] && run q16_108 108 "$B16 --init_from /root/states/T_q16_108/init.pt"
xfer T_q16_156 "$(N156 "$B16")" /root/states/K4_2_sat/init.pt
[ -f /root/states/T_q16_156/init.pt ] && run q16_156 156 "$(N156 "$B16") --init_from /root/states/T_q16_156/init.pt"
kill $SMI
echo "$(date -u +%T) DONE d9" >> $P
