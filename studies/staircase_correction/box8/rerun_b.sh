#!/bin/bash
# box 8: S3 of the B branch again for the runs without a result. 3 runs at a time with the correction do not fit into
# the 16 GB of the V100 (about 5.6 GB for each run, PROTOKOLL 3.39): at most 2 runs on the GPU, a new run only with
# at least 7 GB free; a run that stops with CUDA out of memory starts again (at most 3 tries). Same arguments, codes
# and start states as b_s3.sh. Old outputs of a failed try go to $W/failed.
export CUDA_DEVICE=0
P=/root/prof/progress.txt
W=/root/bs3
eval "$(grep '^BARGS=' /root/b_s3.sh)"
I0=$W/prep/init.pt
RUNNING_PIDS="26104 26782"          # b5 null3 and b4 old of b_s3.sh, still on the GPU when this script starts
done_ok() { [ -f $W/$1/$2/m.pt ] && [ -f $W/$1/$2/summary.json ]; }
gpu_free() { nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1; }
gpu_apps() { nvidia-smi --query-compute-apps=pid --format=csv,noheader | grep -c .; }
one() { s=$1; name=$2; code=$3; init=$4; D=$W/$s/$name
  for try in 1 2 3; do
    if [ -e $D ] || [ -e $D.log ]; then mkdir -p $W/failed; ts=$(date +%s); mv $D $W/failed/${s}_${name}_$ts 2>/dev/null; mv $D.log $W/failed/${s}_${name}_$ts.log 2>/dev/null; fi
    mkdir -p $D; echo "$(date -u +%T) start BS3R $s $name try $try" >> $P
    cd /root/$code/studies/beta_bpeak
    PYTHONPATH=/root/$code python /root/step_dump.py --at 256 --save $D/m.pt --go_on -- $BARGS --init_from $init --out $D > $D.log 2>&1
    rc=$?
    echo "$(date -u +%T) end BS3R $s $name try $try exit $rc" >> $P
    if [ $rc -eq 0 ] && done_ok $s $name; then return 0; fi
    grep -q "OutOfMemoryError" $D.log || { echo "$(date -u +%T) FAIL BS3R $s $name (not out of memory)" >> $P; return 1; }
    sleep 120
  done
  echo "$(date -u +%T) FAIL BS3R $s $name (out of memory 3 times)" >> $P; return 1; }
export -f one done_ok; export W P BARGS CUDA_DEVICE
launch_missing() {
  while read s name code init; do
    done_ok $s $name && continue
    case " $SKIP " in *" $s/$name "*) continue;; esac
    while [ $(gpu_apps) -ge 2 ] || [ $(gpu_free) -lt 7000 ]; do sleep 20; done
    one $s $name $code $init < /dev/null &
    sleep 150                                     # the set-ups of the correction do not overlap
  done < $W/jobs.txt
  wait
}
echo "$(date -u +%T) start BS3R (at most 2 runs on the GPU)" >> $P
SKIP="b5/null3 b4/old" launch_missing
for p in $RUNNING_PIDS; do while kill -0 $p 2>/dev/null; do sleep 30; done; done
SKIP="" launch_missing                            # second pass: b5 null3 and b4 old if they failed
for s in p10 b5 b4 b1; do python /root/s3_eval.py --old $W/$s/old --null $W/$s/null1 $W/$s/null2 $W/$s/null3 --new $W/$s/new --init_old $I0 --init_null $W/null1/init.pt $W/null2/init.pt $W/null3/init.pt --json $W/$s.json > $W/$s.txt 2>&1; echo "$(date -u +%T) eval BS3R $s exit $?" >> $P; done
touch $W/DONE
