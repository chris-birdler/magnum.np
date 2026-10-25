#!/bin/bash
# box 11, PLAN 14b D1b: T1 of B_d53 (CSR product) against B_d53P (PairBlocks) on the 16 states of refB.pt, double
# precision; the two set-ups at the same time.
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
T1=/root/t1/BP; RES=/root/tres_p; mkdir -p $T1 $RES
r64() {  # NAME TREE
  n=$1; t=$2; [ -f $T1/$n/state_240.pt ] && return
  echo "$(date -u +%T) start T1 BP $n" >> $P; cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/t1_BP_$n PYTHONPATH=/root/$t python /root/s/s1_states.py --ref /root/k1/B/refB.pt --out $T1/$n -- $ARGS_B --b_list 0.006 --out /root/t1/run_BP_$n > $T1/$n.log 2>&1
  rc=$?; echo "$(date -u +%T) end T1 BP $n exit $rc" >> $P; }
r64 B_d53 B_d53 & r64 B_d53P B_d53P & wait
python /root/s/k_eval.py --t1 D1b $T1/B_d53 $T1/B_d53P --json $RES/t1_D1b.json > $RES/t1_D1b.txt 2>&1
rc=$?; echo "$(date -u +%T) T1 D1b exit $rc" >> $P
rm -rf $T1/B_d53
echo "$(date -u +%T) DONE d1b" >> $P
