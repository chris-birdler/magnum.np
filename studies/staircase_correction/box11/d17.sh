#!/bin/bash
# box 11, PLAN 15 WQ (PROTOKOLL 3.44), after d16.sh: the v-weighted symmetry of the pair operator at q 16 (information,
# as test_correction_is_symmetric; 2 random fields; float64) on the geometry of the WQ runs q16_108 and q16_156 (if
# that run has a config.json), with the tree AS3.
export CUDA_DEVICE=0
P=/root/progress.txt
while ! grep -q "DONE d16$" $P; do sleep 30; done
Q=/root/wq
for c in q16_108 q16_156; do
  f=$Q/run_$c/config.json; [ -f $f ] || continue; [ -f $Q/sym_$c.json ] && continue
  echo "$(date -u +%T) start D17 sym $c" >> $P
  cd /root/s
  PYTHONPATH=/root/AS3:/root/AS3/studies/beta_bpeak timeout 7200 python /root/s/sym_check.py --config $f --json $Q/sym_$c.json > $Q/sym_$c.log 2>&1
  rc=$?; echo "$(date -u +%T) end D17 sym $c exit $rc" >> $P
done
echo "$(date -u +%T) DONE d17" >> $P
