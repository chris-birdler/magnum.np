#!/bin/bash
# box 11, PLAN 15 F0, second run (PROTOKOLL 3.44): the first run (d10b.sh) gave no stacks (torch 2.7 needs the verbose
# profiler config); prof_step.py with that fix. q4_108, q4_156 and off156, 10 profiled sample intervals. Waits for
# the end of d10c.sh; then the fixed prof_step.py replaces the old file (same behaviour without --stack).
export CUDA_DEVICE=0
P=/root/progress.txt
while ! grep -q "DONE d10c$" $P; do sleep 30; done
mv /root/s/prof_step_new.py /root/s/prof_step.py; echo "$(date -u +%T) D10d prof_step.py $(md5sum /root/s/prof_step.py | cut -c1-32)" >> $P
eval "$(grep '^ARGS_A=' /root/s/k1.sh)"; eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
A="${ARGS_A/--precision_cycles double /}"; A="${A/--init_from \/root\/states\/K1_2_sat\/init.pt/}"
B="${ARGS_B/--precision_cycles double /}"; B="${B/--init_from \/root\/states\/Bprep\/init.pt/}"
N156() { x="${1/--d_lex 212.0/--d_lex 300.0}"; x="${x/--phi 0.62/--phi 0.55}"; x="${x/--Leff_lex 30.0/--Leff_lex 18.0}"; echo "${x/--seed 2013/--seed 2043}"; }
W=/root/d10; mkdir -p $W
run() {  # NAME N ARGS
  n=$1; [ -f $W/$n.json ] && return
  echo "$(date -u +%T) start D10d $n" >> $P
  cd /root/ABP/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d10_$n PYTHONPATH=/root/ABP python /root/s/prof_step.py --skip 30 --plain 20 --window 10 --stack --kcache /root/kc/c10_ABP_$2 --json $W/$n.json -- $3 --b_list 0.006 --out $W/run_$n > $W/$n.log 2>&1
  rc=$?; echo "$(date -u +%T) end D10d $n exit $rc" >> $P
}
run f0b_q4_108 108 "$B --init_from /root/states/Bprep/init.pt"
run f0b_q4_156 156 "$(N156 "$B") --init_from /root/states/T_q4_156/init.pt"
run f0b_off156 156 "$(N156 "$A") --init_from /root/states/K4_2_sat/init.pt"
echo "$(date -u +%T) DONE d10d" >> $P
