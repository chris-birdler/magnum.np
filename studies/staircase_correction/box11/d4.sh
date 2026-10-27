#!/bin/bash
# box 11, PLAN 14b D4: outlier test 3 on the patched tree. P arm (B_d53P, --k1_thr 0 0: each difference recorded): the
# job of H2_d53 twice, one run after the other (H3_d53P_a, H3_d53P_b). Rate control (B_d53, without the patch, the
# thresholds of outlier test 3): the job of H2_d53 once (H3_d53). Job: refB.pt, 17 windows, 18 repetitions,
# --k1_redo all. 2 processes on the GPU all the time: an arm whose main runs end first runs fillers of 2 repetitions
# (new output files) until the main runs of the other arm have ended; the filler that runs then is stopped by its PID
# (its log is read for events, its exposure is not counted). A main run with an exit code other than 0 ends its arm
# (file NAME.FAILED); the other arm then stops its fillers.
export CUDA_DEVICE=0
P=/root/progress.txt
eval "$(grep '^ARGS_B=' /root/s/k1.sh)"
W=/root/k1/H; mkdir -p $W
start() {  # NAME TREE REPS THR...: one step_dump run in the background; its PID in $W/NAME.pid
  n=$1; t=$2; r=$3; shift 3
  echo "$(date -u +%T) start D4 $n" >> $P
  cd /root/$t/studies/beta_bpeak
  TORCHINDUCTOR_CACHE_DIR=/root/ic/H/$n PYTHONPATH=/root/$t python /root/s/step_dump.py --k1_arm /root/k1/B/refB.pt --k1_out $W/$n.pt --k1_reps $r --k1_redo all "$@" -- $ARGS_B --b_list 0.006 --out $W/run_$n > $W/$n.log 2>&1 &
  echo $! > $W/$n.pid
}
arm() {  # ARM TREE OTHER THR MAIN_RUN_NAMES...
  A=$1; t=$2; o=$3; thr=$4; shift 4
  for n in "$@"; do
    start $n $t 18 $thr; pid=$(cat $W/$n.pid); wait $pid; e=$?; echo "$(date -u +%T) end D4 $n exit $e" >> $P
    [ $e -ne 0 ] && { touch $W/$A.FAILED $W/$A.MAIN_DONE; return; }
  done
  touch $W/$A.MAIN_DONE
  i=1
  while [ ! -f $W/$o.MAIN_DONE ]; do
    n=${A}_f$i; start $n $t 2 $thr; pid=$(cat $W/$n.pid)
    while kill -0 $pid 2>/dev/null; do
      if [ -f $W/$o.MAIN_DONE ]; then
        kill $pid; echo "$(date -u +%T) D4 filler $n stopped by PID $pid (main runs of $o ended; exposure not counted)" >> $P; break
      fi
      sleep 20
    done
    wait $pid; rc=$?; echo "$(date -u +%T) end D4 $n exit $rc" >> $P
    i=$((i + 1))
  done
}
arm H3_d53P B_d53P H3_d53 "--k1_thr 0 0" H3_d53P_a H3_d53P_b &
arm H3_d53 B_d53 H3_d53P "" H3_d53 &
wait
cd $W
python /root/s/glitch_count.py --reps $(ls H3_d53P_*.pt) --term_rel 0 --json H3_d53P_count.json > H3_d53P_count.txt 2>&1
rc=$?; echo "$(date -u +%T) D4 count P arm exit $rc" >> $P
python /root/s/glitch_count.py --reps $(ls H3_d53.pt H3_d53_f*.pt 2>/dev/null) --json H3_d53_count.json > H3_d53_count.txt 2>&1
rc=$?; echo "$(date -u +%T) D4 count rate control exit $rc" >> $P
grep -c "redo event" H3_d53P*.log > H3_d53P_events_in_logs.txt 2>&1
echo "$(date -u +%T) DONE d4" >> $P
