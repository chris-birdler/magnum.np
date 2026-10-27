#!/bin/bash
# box 11, PLAN 14b D7 (speed measurement, start states with the correction): 3 transfer runs as for bs3/prep (box8/b_s3.sh):
# K*_2_sat transferred onto the fraction surface with the correction, at most 6000 relaxation iterations, single
# precision; tree B_det (as bs3/prep with B_df0). q 4 at N 156, q 8 at N 108, q 8 at N 156. One at a time.
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
B="${ARGS_B/--precision_cycles double /}"; B="${B/--init_from \/root\/states\/Bprep\/init.pt/}"
N156() { x="${1/--d_lex 212.0/--d_lex 300.0}"; x="${x/--phi 0.62/--phi 0.55}"; x="${x/--Leff_lex 30.0/--Leff_lex 18.0}"; echo "${x/--seed 2013/--seed 2043}"; }
xfer() {  # NAME ARGS SRC
  n=$1; a=$2; src=$3; W=/root/states/$n; [ -f $W/init.pt ] && return; mkdir -p $W
  echo "$(date -u +%T) start D7 transfer $n" >> $P
  cd /root/B_det/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d7t_$n PYTHONPATH=/root/B_det python /root/s/step_dump.py --at 1 --save $W/prep_m.pt -- $a --b_list 0.006 --init_transfer $src --out $W > $W.log 2>&1
  rc=$?; echo "$(date -u +%T) end D7 transfer $n exit $rc ($(ls $W/init.pt 2>/dev/null))" >> $P
}
xfer T_q4_156 "$(N156 "$B")" /root/states/K4_2_sat/init.pt
xfer T_q8_108 "${B/--staircase_q 4/--staircase_q 8}" /root/states/K1_2_sat/init.pt
xfer T_q8_156 "$(N156 "${B/--staircase_q 4/--staircase_q 8}")" /root/states/K4_2_sat/init.pt
echo "$(date -u +%T) DONE d7t" >> $P
