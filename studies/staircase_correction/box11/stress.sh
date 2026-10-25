#!/bin/bash
# box 11: stress test of the field parts on one fixed state (transient outlier of the replay, PLAN 14a); 2 runs at a time
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
mkdir -p /root/stress
one() { t=$1; shift
  echo "$(date -u +%T) start stress $t" >> $P; cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/st_$t PYTHONPATH=/root/$t python /root/s/glitch_stress.py --ref /root/k1/B/refB.pt --k0 8 --parts "$@" --idx 1513981 876667 --json /root/stress/$t.json -- $ARGS_B --b_list 0.006 --out /root/stress/run_$t > /root/stress/$t.log 2>&1
  rc=$?; echo "$(date -u +%T) end stress $t exit $rc" >> $P; }
export -f one; export P ARGS_B CUDA_DEVICE
printf "B_d53 corr:1800 spmv:600 centroid:600 demag:300 dm:900\nB_c63 corr:1800 centroid:600 demag:300 dm:900\n" | xargs -P 2 -L1 bash -c 'one "$@"' _
echo "$(date -u +%T) DONE stress" >> $P
