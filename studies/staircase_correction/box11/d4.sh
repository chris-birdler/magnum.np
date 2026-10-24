#!/bin/bash
# box 11, PLAN 14b D4: outlier test 3 on the patched tree. H3_d53P (B_d53P) and the control H3_d53 (B_d53, without the
# patch), each the job of H2_d53 (refB.pt, 17 windows, 18 repetitions, --k1_redo all); 2 processes on the GPU all the
# time. An arm whose main run ends first runs fillers of 2 repetitions (new output files; a completed filler counts for
# its arm) until the main run of the other arm has ended; the filler that runs then is stopped by its PID and is not
# counted.
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
W=/root/k1/H; mkdir -p $W
start() {  # NAME TREE REPS: start one step_dump run in the background; its PID in $W/NAME.pid
  echo "$(date -u +%T) start D4 $1" >> $P
  cd /root/$2/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/H/$1 PYTHONPATH=/root/$2 python /root/s/step_dump.py --k1_arm /root/k1/B/refB.pt --k1_out $W/$1.pt --k1_reps $3 --k1_redo all -- $ARGS_B --b_list 0.006 --out $W/run_$1 > $W/$1.log 2>&1 &
  echo $! > $W/$1.pid
}
arm() {  # NAME TREE OTHER
  start $1 $2 18; pid=$(cat $W/$1.pid); wait $pid; echo "$(date -u +%T) end D4 $1 exit $?" >> $P
  touch $W/$1.MAIN_DONE
  i=1
  while [ ! -f $W/$3.MAIN_DONE ]; do
    n=${1}_f$i; start $n $2 2; pid=$(cat $W/$n.pid)
    while kill -0 $pid 2>/dev/null; do
      if [ -f $W/$3.MAIN_DONE ]; then
        kill $pid; echo "$(date -u +%T) D4 filler $n stopped by PID $pid (main run of $3 ended; not counted)" >> $P; break
      fi
      sleep 20
    done
    wait $pid; echo "$(date -u +%T) end D4 $n exit $?" >> $P
    i=$((i + 1))
  done
}
arm H3_d53P B_d53P H3_d53 &
arm H3_d53 B_d53 H3_d53P &
wait
cd $W && python /root/s/glitch_count.py --reps H3_d53P*.pt --json H3_d53P_count.json > H3_d53P_count.txt 2>&1
echo "$(date -u +%T) D4 count H3_d53P exit $?" >> $P
cd $W && python /root/s/glitch_count.py --reps $(ls H3_d53.pt H3_d53_f*.pt 2>/dev/null) --json H3_d53_count.json > H3_d53_count.txt 2>&1
echo "$(date -u +%T) D4 count H3_d53 exit $?" >> $P
echo "$(date -u +%T) DONE d4" >> $P
