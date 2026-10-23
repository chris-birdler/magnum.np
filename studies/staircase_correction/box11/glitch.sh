#!/bin/bash
# box 11: test for non-deterministic faults of the field terms (finding of the K1 pilot); 2 runs at a time
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"; B32="${ARGS_B/--precision_cycles double /}"
mkdir -p /root/glitch
one() { t=$1; pr=$2; [ $pr = 64 ] && A="$ARGS_B" || A="$B32"
  echo "$(date -u +%T) start glitch $t $pr" >> $P; cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/gl_${t}_$pr PYTHONPATH=/root/$t python /root/s/glitch_test.py --ref /root/k1/B/refB.pt --k0 0 48 144 --reps 300 --json /root/glitch/${t}_$pr.json -- $A --b_list 0.006 --out /root/glitch/run_${t}_$pr > /root/glitch/${t}_$pr.log 2>&1
  echo "$(date -u +%T) end glitch $t $pr exit $?" >> $P; }
export -f one; export P ARGS_B B32 CUDA_DEVICE
printf "B_c63 64\nB_01e 64\nB_d53 64\nB_df0 64\nB_c63 32\nB_01e 32\nB_d53 32\nB_df0 32\n" | xargs -P 2 -L1 bash -c 'one "$@"' _
echo "$(date -u +%T) DONE glitch" >> $P
