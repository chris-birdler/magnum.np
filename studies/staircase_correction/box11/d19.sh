#!/bin/bash
# box 11, PLAN 15 WQ (PROTOKOLL 3.44): method check of the pair asymmetry at q 16 (sym_check 4.4e-11 > 1e-12): the same
# check on the same geometry (N 108, config of q16_108) at q 2, 4 and 8, tree AS3, on the GPU. If the asymmetry grows
# with q, it is the rounding of the coefficient sums (q^6 sub-cell pairs per coefficient).
export CUDA_DEVICE=0
P=/root/progress.txt; Q=/root/wq
for q in 2 4 8; do
  [ -f $Q/sym_q${q}_108.json ] && continue
  echo "$(date -u +%T) start D19 sym q$q 108" >> $P
  cd /root/s
  PYTHONPATH=/root/AS3:/root/AS3/studies/beta_bpeak timeout 1800 python /root/s/sym_check.py --config $Q/run_q16_108/config.json --q $q --json $Q/sym_q${q}_108.json > $Q/sym_q${q}_108.log 2>&1
  rc=$?; echo "$(date -u +%T) end D19 sym q$q 108 exit $rc" >> $P
done
echo "$(date -u +%T) DONE d19" >> $P
