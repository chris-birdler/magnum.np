#!/bin/bash
# box 11: stress test 3 of the sparse product (outlier test 3 showed the fault in it); two processes with the same job, as in all observed outliers; 3 parts of 1 h; --idx gives row and volume fraction of the event cells
# usage: stress3.sh a|b
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
mkdir -p /root/stress3
t=B_d53; x=$1
echo "$(date -u +%T) start stress3 $t$x" >> $P; cd /root/$t/studies/beta_bpeak
TORCHINDUCTOR_CACHE_DIR=/root/ic/st3_$t$x PYTHONPATH=/root/$t python /root/s/glitch_stress.py --ref /root/k1/B/refB.pt --k0 8 144 --parts spmv:3600 spmv:3600 spmv:3600 --idx 1513981 876667 1634169 3060786 3615150 3632484 2605404 --json /root/stress3/$t$x.json -- $ARGS_B --b_list 0.006 --out /root/stress3/run_$t$x > /root/stress3/$t$x.log 2>&1
echo "$(date -u +%T) end stress3 $t$x exit $?" >> $P
