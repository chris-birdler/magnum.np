#!/bin/bash
# box 11, PLAN 14b D3: detector control and GPU smoke on B_d53P: a fault of 1e-3 (relative) at the largest entry of the
# PairBlocks product in 4 chosen calls inside the tries (k1_inject.py --spmv); windows k0 8 and 144; --k1_redo all
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
W=/root/k1/D3; mkdir -p $W
echo "$(date -u +%T) start D3 D3_inj" >> $P
cd /root/B_d53P/studies/beta_bpeak
TORCHINDUCTOR_CACHE_DIR=/root/ic/d3 PYTHONPATH=/root/B_d53P python /root/s/k1_inject.py --spmv 1e-3 --calls 100 700 1300 2000 -- --k1_arm /root/k1/B/refB.pt --k1_out $W/D3_inj.pt --k1_ck 8 144 --k1_reps 1 --k1_redo all -- $ARGS_B --b_list 0.006 --out $W/run_D3_inj > $W/D3_inj.log 2>&1
echo "$(date -u +%T) end D3 D3_inj exit $?" >> $P
python /root/s/glitch_count.py --reps $W/D3_inj.pt --json $W/D3_inj_count.json > $W/D3_inj_count.txt 2>&1
echo "$(date -u +%T) D3 count exit $?" >> $P
