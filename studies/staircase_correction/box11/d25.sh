#!/bin/bash
# box 11, PLAN 16 WQ, second try (PROTOKOLL 3.46), after d24.sh. d22.sh stopped at once with "K4_2_sat/init.pt was
# relaxed with other arguments: sat_static True -> False, h_sat 2.0 -> None": the same error of my script as the first
# try of WQ in PLAN 15 (3.44): the run_loops of the AS trees compares sat_static and h_sat at a transfer. The fix of
# 3.44, not a code change: the transfer as d7t.sh / d9.sh / d16.sh (step_dump --at 1 with the study folder of B_det),
# with the magnumnp package of AS5: the tree ASW5 is a copy of AS5 with the study folder of B_det. Then, with AS5:
# the speed run q16_156 (prof_step.py as W5, 3 runs; the first is the WQ run), the pair symmetry (sym_check.py). The
# GPU memory and the memory of the container are sampled every 5 s. Time limit 6 h from the start of this part.
export CUDA_DEVICE=0
P=/root/progress.txt
while ! grep -q "DONE d24$" $P; do sleep 30; done
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
B="${ARGS_B/--precision_cycles double /}"; B="${B/--init_from \/root\/states\/Bprep\/init.pt/}"
B16="${B/--staircase_q 4/--staircase_q 16}"
N156() { x="${1/--d_lex 212.0/--d_lex 300.0}"; x="${x/--phi 0.62/--phi 0.55}"; x="${x/--Leff_lex 30.0/--Leff_lex 18.0}"; echo "${x/--seed 2013/--seed 2043}"; }
if [ ! -d /root/ASW5 ]; then
  cp -a /root/AS5 /root/ASW5
  rsync -a --delete /root/B_det/studies/beta_bpeak/ /root/ASW5/studies/beta_bpeak/
fi
W=/root/p16/wq; mkdir -p $W
( cd /root/ASW5 && find magnumnp -name "*.py" | sort | xargs md5sum ) > $W/ASW5_pkg.md5
( cd /root/AS5 && md5sum -c --quiet $W/ASW5_pkg.md5 ) > $W/ASW5_pkg_check.txt 2>&1; echo "$(date -u +%T) P16 WQ ASW5 package = AS5 package: exit $?" >> $P
md5sum /root/ASW5/studies/beta_bpeak/run_loops.py /root/B_det/studies/beta_bpeak/run_loops.py > $W/ASW5_md5.txt
T0=$(date +%s); LIM=$((6 * 3600))
left() { echo $(( LIM - ($(date +%s) - T0) )); }
nvidia-smi --query-gpu=timestamp,memory.used --format=csv,noheader -l 5 > $W/smi2.csv 2>&1 &
SMI=$!
CG=/sys/fs/cgroup/memory
( while true; do echo "$(date -u +%T) $(cat $CG/memory.usage_in_bytes 2>/dev/null) $(cat $CG/memory.max_usage_in_bytes 2>/dev/null)"; sleep 5; done ) > $W/cgroup2.txt 2>&1 &
CGP=$!
S=/root/states/T17_q16_156
if [ ! -f $S/init.pt ]; then
  rm -rf $S; mkdir -p $S
  echo "$(date -u +%T) start P16 WQ transfer2 (limit $(left) s)" >> $P
  cd /root/ASW5/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/p16_wq_xfer2 PYTHONPATH=/root/ASW5 timeout $(left) python /root/s/step_dump.py --at 1 --save $S/prep_m.pt -- $(N156 "$B16") --b_list 0.006 --init_transfer /root/states/K4_2_sat/init.pt --out $S > $S.log 2>&1
  echo "$(date -u +%T) end P16 WQ transfer2 exit $? ($(ls $S/init.pt 2>/dev/null))" >> $P
fi
run() {  # NAME
  n=$1; [ -f $W/$n.json ] && return; [ $(left) -le 60 ] && return
  echo "$(date -u +%T) start P16 WQ speed $n (limit $(left) s)" >> $P
  cd /root/AS5/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/p16_wq_$n PYTHONPATH=/root/AS5 timeout $(left) python /root/s/prof_step.py --skip 30 --plain 96 --window 0 --trace --kcache /root/kc/AS5_156 --json $W/$n.json -- $(N156 "$B16") --init_from $S/init.pt --b_list 0.006 --out $W/run_$n > $W/$n.log 2>&1
  echo "$(date -u +%T) end P16 WQ speed $n exit $?" >> $P
}
if [ -f $S/init.pt ]; then
  run q16_156_AS5_r1
  if [ -f $W/run_q16_156_AS5_r1/config.json ] && [ ! -f $W/sym_q16_156.json ] && [ $(left) -gt 60 ]; then
    echo "$(date -u +%T) start P16 WQ sym q16_156 (limit $(left) s)" >> $P
    cd /root/s
    PYTHONPATH=/root/AS5:/root/AS5/studies/beta_bpeak timeout $(left) python /root/s/sym_check.py --config $W/run_q16_156_AS5_r1/config.json --json $W/sym_q16_156.json > $W/sym_q16_156.log 2>&1
    echo "$(date -u +%T) end P16 WQ sym q16_156 exit $?" >> $P
  fi
  run q16_156_AS5_r2
  run q16_156_AS5_r3
fi
kill $SMI $CGP
echo "$(date -u +%T) DONE d25" >> $P
