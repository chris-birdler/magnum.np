#!/bin/bash
# box 11, PLAN 16 WQ (PROTOKOLL 3.46), after d20.sh and d23.sh: q 16 at N 156 with AS5 (information, n 1). (1) The
# transfer state (run_loops --prep_only --init_transfer from K4_2_sat, as d15.sh, tree AS5); (2) a speed run from it
# (prof_step.py as W5); (3) the pair symmetry (sym_check.py). The GPU memory (nvidia-smi) and the memory of the
# container (cgroup: usage and its maximum) are sampled every 5 s. Time limit: 6 h from the start; then timeout stops
# the running process.
export CUDA_DEVICE=0
P=/root/progress.txt
while ! grep -q "DONE d23$" $P; do sleep 30; done
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
B="${ARGS_B/--precision_cycles double /}"; B="${B/--init_from \/root\/states\/Bprep\/init.pt/}"
B16="${B/--staircase_q 4/--staircase_q 16}"
N156() { x="${1/--d_lex 212.0/--d_lex 300.0}"; x="${x/--phi 0.62/--phi 0.55}"; x="${x/--Leff_lex 30.0/--Leff_lex 18.0}"; echo "${x/--seed 2013/--seed 2043}"; }
W=/root/p16/wq; mkdir -p $W
T0=$(date +%s); LIM=$((6 * 3600))
left() { echo $(( LIM - ($(date +%s) - T0) )); }
nvidia-smi --query-gpu=timestamp,memory.used --format=csv,noheader -l 5 > $W/smi.csv 2>&1 &
SMI=$!
CG=/sys/fs/cgroup/memory
( while true; do echo "$(date -u +%T) $(cat $CG/memory.usage_in_bytes 2>/dev/null) $(cat $CG/memory.max_usage_in_bytes 2>/dev/null)"; sleep 5; done ) > $W/cgroup.txt 2>&1 &
CGP=$!
S=/root/states/T17_q16_156; mkdir -p $S
if [ ! -f $S/init.pt ]; then
  echo "$(date -u +%T) start P16 WQ transfer (limit $(left) s)" >> $P
  cd /root/AS5/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/p16_wq_xfer PYTHONPATH=/root/AS5 timeout $(left) python run_loops.py $(N156 "$B16") --b_list 0.006 --init_transfer /root/states/K4_2_sat/init.pt --prep_only --out $S > $S.log 2>&1
  echo "$(date -u +%T) end P16 WQ transfer exit $? ($(ls $S/init.pt 2>/dev/null))" >> $P
fi
if [ -f $S/init.pt ] && [ ! -f $W/q16_156.json ] && [ $(left) -gt 60 ]; then
  echo "$(date -u +%T) start P16 WQ speed q16_156 (limit $(left) s)" >> $P
  cd /root/AS5/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/p16_wq_q16_156 PYTHONPATH=/root/AS5 timeout $(left) python /root/s/prof_step.py --skip 30 --plain 96 --window 0 --trace --kcache /root/kc/AS5_156 --json $W/q16_156.json -- $(N156 "$B16") --init_from $S/init.pt --b_list 0.006 --out $W/run_q16_156 > $W/q16_156.log 2>&1
  echo "$(date -u +%T) end P16 WQ speed q16_156 exit $?" >> $P
fi
if [ -f $W/run_q16_156/config.json ] && [ ! -f $W/sym_q16_156.json ] && [ $(left) -gt 60 ]; then
  echo "$(date -u +%T) start P16 WQ sym q16_156 (limit $(left) s)" >> $P
  cd /root/s
  PYTHONPATH=/root/AS5:/root/AS5/studies/beta_bpeak timeout $(left) python /root/s/sym_check.py --config $W/run_q16_156/config.json --json $W/sym_q16_156.json > $W/sym_q16_156.log 2>&1
  echo "$(date -u +%T) end P16 WQ sym q16_156 exit $?" >> $P
fi
kill $SMI $CGP
echo "$(date -u +%T) DONE d22" >> $P
