#!/bin/bash
# box 11, PLAN 15 (PROTOKOLL 3.44), after d15.sh: (1) the halves rule of W5 in off108 (AS3 r1 14.7 %, r2 11.5 % > 10 %):
# the configuration repeats with 100 sample intervals of warm-up, ABP and AS3 in turn, n 3 for each code (both codes,
# so that the gain compares the same intervals; for the speed only). (2) WQ again: d15.sh stopped at the transfer check
# of the run_loops of AS3 (the source states have sat_static True, h_sat 2.0; the q 4 and q 8 transfer states of W5 came
# from step_dump --at 1 with the tree B_det, which has no such check). The transfer now is that of d9.sh / d7t.sh
# (step_dump --at 1, the study folder of B_det) with the magnumnp package of AS3 (FQ): the tree ASW is a copy of AS3
# with the study folder of B_det. The speed run of AS3 from it as in d15.sh. Time limit of WQ: 6 h.
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_A=' /root/s/k1.sh)"; eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
A="${ARGS_A/--precision_cycles double /}"; A="${A/--init_from \/root\/states\/K1_2_sat\/init.pt/}"
B="${ARGS_B/--precision_cycles double /}"; B="${B/--init_from \/root\/states\/Bprep\/init.pt/}"
B16="${B/--staircase_q 4/--staircase_q 16}"
N156() { x="${1/--d_lex 212.0/--d_lex 300.0}"; x="${x/--phi 0.62/--phi 0.55}"; x="${x/--Leff_lex 30.0/--Leff_lex 18.0}"; echo "${x/--seed 2013/--seed 2043}"; }
W=/root/w5
run() {  # NAME TREE
  n=$1; t=$2; [ -f $W/$n.json ] && return
  echo "$(date -u +%T) start D16 $n" >> $P
  cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d16_$n PYTHONPATH=/root/$t python /root/s/prof_step.py --skip 100 --plain 96 --window 0 --trace --kcache /root/kc/${t}_108 --json $W/$n.json -- $A --init_from /root/states/K1_2_sat/init.pt --b_list 0.006 --out $W/run_$n > $W/$n.log 2>&1
  rc=$?; echo "$(date -u +%T) end D16 $n exit $rc" >> $P
}
for r in 1 2 3; do for t in ABP AS3; do run off108_${t}_s100_r$r $t; done; done
# (2) WQ
if [ ! -d /root/ASW ]; then
  cp -a /root/AS3 /root/ASW
  rsync -a --delete /root/B_det/studies/beta_bpeak/ /root/ASW/studies/beta_bpeak/
fi
md5sum /root/ASW/magnumnp/field_terms/demag_staircase.py /root/AS3/magnumnp/field_terms/demag_staircase.py \
  /root/ASW/studies/beta_bpeak/run_loops.py /root/B_det/studies/beta_bpeak/run_loops.py > /root/wq/ASW_md5.txt
Q=/root/wq
T0=$(date +%s); LIM=$((6 * 3600))
left() { echo $(( LIM - ($(date +%s) - T0) )); }
nvidia-smi --query-gpu=timestamp,memory.used --format=csv,noheader -l 5 > $Q/smi2.csv 2>&1 &
SMI=$!
xfer() {  # NAME ARGS SRC
  n=$1; S=/root/states/$n; [ -f $S/init.pt ] && return; [ $(left) -le 60 ] && return
  rm -rf $S; mkdir -p $S
  echo "$(date -u +%T) start D16 WQ transfer $n (limit $(left) s)" >> $P
  cd /root/ASW/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d16t_$n PYTHONPATH=/root/ASW timeout $(left) python /root/s/step_dump.py --at 1 --save $S/prep_m.pt -- $2 --b_list 0.006 --init_transfer $3 --out $S > $S.log 2>&1
  rc=$?; echo "$(date -u +%T) end D16 WQ transfer $n exit $rc ($(ls $S/init.pt 2>/dev/null))" >> $P
}
wq() {  # NAME N ARGS
  n=$1; [ -f $Q/$n.json ] && return; [ $(left) -le 60 ] && return
  echo "$(date -u +%T) start D16 WQ $n (limit $(left) s)" >> $P
  cd /root/AS3/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d16_$n PYTHONPATH=/root/AS3 timeout $(left) python /root/s/prof_step.py --skip 30 --plain 96 --window 0 --trace --kcache /root/kc/AS3_$2 --json $Q/$n.json -- $3 --b_list 0.006 --out $Q/run_$n > $Q/$n.log 2>&1
  rc=$?; echo "$(date -u +%T) end D16 WQ $n exit $rc" >> $P
}
xfer T16_q16_108 "$B16" /root/states/K1_2_sat/init.pt
[ -f /root/states/T16_q16_108/init.pt ] && wq q16_108 108 "$B16 --init_from /root/states/T16_q16_108/init.pt"
xfer T16_q16_156 "$(N156 "$B16")" /root/states/K4_2_sat/init.pt
[ -f /root/states/T16_q16_156/init.pt ] && wq q16_156 156 "$(N156 "$B16") --init_from /root/states/T16_q16_156/init.pt"
kill $SMI
echo "$(date -u +%T) DONE d16" >> $P
