#!/bin/bash
# box 11, PLAN 14b D1: bench_det.py on B_det (e7ce6ca), 2 processes at the same time, each with its own compile cache
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
mkdir -p /root/d1
one() {
  x=$1
  echo "$(date -u +%T) start D1 bench $x" >> $P
  cd /root/B_det/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/d1_$x PYTHONPATH=/root/B_det python /root/s/bench_det.py --ref /root/k1/B/refB.pt --k0 8 144 --reps 20000 --json /root/d1/bench_$x.json -- $ARGS_B --b_list 0.006 --out /root/d1/run_$x > /root/d1/bench_$x.log 2>&1
  echo "$(date -u +%T) end D1 bench $x exit $?" >> $P
}
one a & one b & wait
echo "$(date -u +%T) DONE d1" >> $P
