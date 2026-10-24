#!/bin/bash
# box 11, PLAN 14b D3 (2 processes at the same time). D3a: detector control on B_d53P, a fault of 1e-3 (relative) at the
# largest entry of the PairBlocks product in 4 chosen calls inside the tries (k1_inject.py --spmv); D3b: ABP in the
# production precision (the cycle in single precision), 2 repetitions. Both: windows k0 8 and 144, --k1_redo all,
# --k1_thr 0 0 (each difference recorded).
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
B32="${ARGS_B/--precision_cycles double /}"
W=/root/k1/D3; mkdir -p $W
d3a() {
  echo "$(date -u +%T) start D3a" >> $P
  cd /root/B_d53P/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d3a PYTHONPATH=/root/B_d53P python /root/s/k1_inject.py --spmv 1e-3 --calls 100 700 1300 2000 -- --k1_arm /root/k1/B/refB.pt --k1_out $W/D3a.pt --k1_ck 8 144 --k1_reps 1 --k1_redo all --k1_thr 0 0 -- $ARGS_B --b_list 0.006 --out $W/run_D3a > $W/D3a.log 2>&1
  echo "$(date -u +%T) end D3a exit $?" >> $P
  python /root/s/glitch_count.py --reps $W/D3a.pt --term_rel 0 --json $W/D3a_count.json > $W/D3a_count.txt 2>&1
}
d3b() {
  echo "$(date -u +%T) start D3b" >> $P
  cd /root/ABP/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d3b PYTHONPATH=/root/ABP python /root/s/step_dump.py --k1_arm /root/k1/B/refB.pt --k1_out $W/D3b.pt --k1_ck 8 144 --k1_reps 2 --k1_redo all --k1_thr 0 0 -- $B32 --b_list 0.006 --out $W/run_D3b > $W/D3b.log 2>&1
  echo "$(date -u +%T) end D3b exit $?" >> $P
  python /root/s/glitch_count.py --reps $W/D3b.pt --term_rel 0 --json $W/D3b_count.json > $W/D3b_count.txt 2>&1
}
d3a & d3b & wait
echo "$(date -u +%T) DONE d3" >> $P
