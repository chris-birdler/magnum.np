#!/bin/bash
# box 11, PLAN 15 WQ (PROTOKOLL 3.44), after d14.sh: q 16 with AS3 (FQ). For N 108 and N 156 as the q 8 arms of D7, only
# --staircase_q 16: (1) the transfer state (run_loops --prep_only --init_transfer, the preparation of production with
# the tree AS3; d9.sh used step_dump --at 1 with B_det, which has no FQ); (2) a speed run of AS3 from it (prof_step.py
# as W5, n 1). nvidia-smi records the GPU memory every 5 s. Time limit: 6 h from the start of WQ.
export CUDA_DEVICE=0
P=/root/progress.txt
while ! grep -q "DONE d14$" $P; do sleep 30; done
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
B="${ARGS_B/--precision_cycles double /}"; B="${B/--init_from \/root\/states\/Bprep\/init.pt/}"
B16="${B/--staircase_q 4/--staircase_q 16}"
N156() { x="${1/--d_lex 212.0/--d_lex 300.0}"; x="${x/--phi 0.62/--phi 0.55}"; x="${x/--Leff_lex 30.0/--Leff_lex 18.0}"; echo "${x/--seed 2013/--seed 2043}"; }
W=/root/wq; mkdir -p $W
T0=$(date +%s); LIM=$((6 * 3600))
left() { echo $(( LIM - ($(date +%s) - T0) )); }
nvidia-smi --query-gpu=timestamp,memory.used --format=csv,noheader -l 5 > $W/smi.csv 2>&1 &
SMI=$!
xfer() {  # NAME ARGS SRC
  n=$1; S=/root/states/$n; [ -f $S/init.pt ] && return; [ $(left) -le 60 ] && return
  echo "$(date -u +%T) start WQ transfer $n (limit $(left) s)" >> $P
  cd /root/AS3/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/wq_$n PYTHONPATH=/root/AS3 timeout $(left) python run_loops.py $2 --b_list 0.006 --init_transfer $3 --prep_only --out $S > $S.log 2>&1
  rc=$?; echo "$(date -u +%T) end WQ transfer $n exit $rc ($(ls $S/init.pt 2>/dev/null))" >> $P
}
run() {  # NAME N ARGS
  n=$1; [ -f $W/$n.json ] && return; [ $(left) -le 60 ] && return
  k=/root/kc/AS3_$2
  echo "$(date -u +%T) start WQ $n (limit $(left) s)" >> $P
  cd /root/AS3/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/wq_$n PYTHONPATH=/root/AS3 timeout $(left) python /root/s/prof_step.py --skip 30 --plain 96 --window 0 --trace --kcache $k --json $W/$n.json -- $3 --b_list 0.006 --out $W/run_$n > $W/$n.log 2>&1
  rc=$?; echo "$(date -u +%T) end WQ $n exit $rc" >> $P
}
xfer T16_q16_108 "$B16" /root/states/K1_2_sat/init.pt
[ -f /root/states/T16_q16_108/init.pt ] && run q16_108 108 "$B16 --init_from /root/states/T16_q16_108/init.pt"
xfer T16_q16_156 "$(N156 "$B16")" /root/states/K4_2_sat/init.pt
[ -f /root/states/T16_q16_156/init.pt ] && run q16_156 156 "$(N156 "$B16") --init_from /root/states/T16_q16_156/init.pt"
kill $SMI
echo "$(date -u +%T) DONE d15" >> $P
