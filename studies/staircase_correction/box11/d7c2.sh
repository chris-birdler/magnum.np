#!/bin/bash
# box 11, PLAN 14b D7, cause of the G fail in q4_108, part 2 (PROTOKOLL 3.41). d7c.sh: B_d53P, B_det and BA4P pass G;
# the excess enters at BA4P -> ABP (A1 ... A3 with the correction on). (1) AB (7f755a0: A1 ... A3 and the B patches
# with the old CSR product, no pair sum), 2 runs: the CSR product is not repeatable, so the 2 runs are 2 realisations
# of the rounding. (2) Null runs of ABP (field terms summed in the reverse and in the rotated order): ABP is bitwise
# repeatable, so its 3 runs of D7 are one realisation. Each run is judged by G as in d7c.sh.
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
B="${ARGS_B/--precision_cycles double /}"; B="${B/--init_from \/root\/states\/Bprep\/init.pt/}"
W=/root/d7c; R=/root/d7cres; D7=/root/d7; mkdir -p $W $R
run() {  # NAME TREE [PROF_STEP_OPTIONS]
  n=$1; t=$2; [ -f $W/$n.json ] && return
  echo "$(date -u +%T) start D7c2 $n" >> $P
  cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d7c_$n PYTHONPATH=/root/$t python /root/s/prof_step.py --skip 30 --plain 96 --window 0 --trace $3 --kcache /root/kc/c_ABP_108 --json $W/$n.json -- $B --init_from /root/states/Bprep/init.pt --b_list 0.006 --out $W/run_$n > $W/$n.log 2>&1
  rc=$?; echo "$(date -u +%T) end D7c2 $n exit $rc" >> $P
}
run q4_108_AB_r1 AB; run q4_108_AB_r2 AB
run q4_108_ABP_null_reverse ABP "--null reverse"; run q4_108_ABP_null_rotate ABP "--null rotate"
for n in AB_r1 AB_r2 ABP_null_reverse ABP_null_rotate; do
  python /root/s/k_eval.py --g_old $D7/q4_108_old_r1.json $D7/q4_108_old_r2.json $D7/q4_108_old_r3.json \
    --g_new $W/q4_108_$n.json --g_null $D7/q4_108_null_reverse.json $D7/q4_108_null_rotate.json \
    --json $R/g_q4_108_$n.json > $R/g_q4_108_$n.txt 2>&1
  rc=$?; echo "$(date -u +%T) D7c2 G q4_108 $n exit $rc" >> $P
done
echo "$(date -u +%T) DONE d7c2" >> $P
