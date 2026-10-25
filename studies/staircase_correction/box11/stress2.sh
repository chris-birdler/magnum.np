#!/bin/bash
# box 11: stress test with two states in turn (stale data shows), production try first; d538bde only (PLAN 14a outlier)
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
mkdir -p /root/stress2
t=B_d53
echo "$(date -u +%T) start stress2 $t" >> $P; cd /root/$t/studies/beta_bpeak
TORCHINDUCTOR_CACHE_DIR=/root/ic/st2_$t PYTHONPATH=/root/$t python /root/s/glitch_stress.py --ref /root/k1/B/refB.pt --k0 8 144 --parts try:3600 corr:1200 spmv:600 --idx 1513981 876667 --json /root/stress2/$t.json -- $ARGS_B --b_list 0.006 --out /root/stress2/run_$t > /root/stress2/$t.log 2>&1
rc=$?; echo "$(date -u +%T) end stress2 $t exit $rc" >> $P
echo "$(date -u +%T) DONE stress2" >> $P
