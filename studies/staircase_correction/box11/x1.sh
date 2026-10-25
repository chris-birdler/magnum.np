#!/bin/bash
# box 11, check X1 (PROTOKOLL 3.41): why do two processes of B_d53P compute slightly different results for the same
# input? Hypothesis: the inductor autotuning of the reduction kernels of block_sum chooses its configuration by time at
# the first call (H3_d53P_a: XBLOCK 32, 4 warps; H3_d53P_b: XBLOCK 128, 8 warps for the sum over the 16 blocks of a
# chunk), and the configuration sets the order of the sum. Each process then keeps its configuration.
# 1. T1 of B_d53P in a new process with a copy of the compile cache of the T1 of D1b: all field terms must be bitwise
#    equal to D1b, else there is a second cause.
# 2. bench_det.py (the product of production) with copies of the caches of H3_d53P_a (two processes) and H3_d53P_b,
#    and one process with an empty cache: the md5 of the products follows the configuration, if the hypothesis holds.
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
T1=/root/t1/BP; X=/root/x1; mkdir -p $X
rm -rf /root/ic/x1_*
cp -a /root/ic/t1_BP_B_d53P /root/ic/x1_t1
cp -a /root/ic/H/H3_d53P_a /root/ic/x1_a1; cp -a /root/ic/H/H3_d53P_a /root/ic/x1_a2; cp -a /root/ic/H/H3_d53P_b /root/ic/x1_b1
echo "$(date -u +%T) start X1 T1" >> $P
cd /root/B_d53P/studies/beta_bpeak
TORCHINDUCTOR_CACHE_DIR=/root/ic/x1_t1 PYTHONPATH=/root/B_d53P python /root/s/s1_states.py --ref /root/k1/B/refB.pt --out $T1/B_d53P_r -- $ARGS_B --b_list 0.006 --out /root/t1/run_BP_B_d53P_r > $X/t1_r.log 2>&1
echo "$(date -u +%T) end X1 T1 exit $?" >> $P
python /root/s/k_eval.py --t1 X1 $T1/B_d53P $T1/B_d53P_r --json $X/t1_X1.json > $X/t1_X1.txt 2>&1
echo "$(date -u +%T) X1 T1 eval exit $?" >> $P
one() {
  x=$1
  echo "$(date -u +%T) start X1 bench $x" >> $P
  cd /root/B_d53P/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/x1_$x PYTHONPATH=/root/B_d53P python /root/s/bench_det.py --ref /root/k1/B/refB.pt --k0 8 144 --reps 20 --time 20 --modes eager --dtypes float64 --json $X/bench_$x.json -- $ARGS_B --b_list 0.006 --out $X/run_$x > $X/bench_$x.log 2>&1
  echo "$(date -u +%T) end X1 bench $x exit $?" >> $P
}
one a1 & one b1 & wait
one a2 & one n & wait
for d in /root/ic/x1_*; do for f in $(find $d -name "*.best_config" | sort); do echo "$(basename $d) $(basename $(dirname $f)) $(cat $f)"; done; done > $X/best_configs.txt
echo "$(date -u +%T) DONE x1" >> $P
