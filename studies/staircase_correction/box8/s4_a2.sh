#!/bin/bash
# box 8: stage 4 of PLAN 11 for the A bundle (old main 9e0ea7c, new perf/timestep 3d9de71), 8 phase-1 states, 4 cycles.
# Same runs as s4_a.sh; changes: the exit code is read directly after the run, a new run starts only with at most
# 2 others on the GPU and at least 5 GB free, a run that stops with CUDA out of memory starts again (at most 3 tries).
export CUDA_DEVICE=0
P=/root/prof/progress.txt
until [ -f /root/bs3/DONE ]; do sleep 30; done
W=/root/s4a; mkdir -p $W
while read st args; do
  mkdir -p $W/$st
  PYTHONPATH=/root/A_old python /root/perturb_init.py /root/states/$st/init.pt $W/$st/null/init.pt --seed 1 > $W/$st/perturb.txt 2>&1
done < /root/s4_jobs.txt
gpu_free() { nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1; }
gpu_apps() { nvidia-smi --query-compute-apps=pid --format=csv,noheader | grep -c .; }
one() { st=$1; arm=$2; code=$3; init=$4; shift 4; D=$W/$st/$arm
  for try in 1 2 3; do
    if [ -e $D ] || [ -e $D.log ]; then mkdir -p $W/failed; ts=$(date +%s); mv $D $W/failed/${st}_${arm}_$ts 2>/dev/null; mv $D.log $W/failed/${st}_${arm}_$ts.log 2>/dev/null; fi
    mkdir -p $D; echo "$(date -u +%T) start S4 $st $arm try $try" >> $P
    cd /root/$code/studies/beta_bpeak
    PYTHONPATH=/root/$code python /root/step_dump.py --at 1024 --save $D/m.pt --go_on -- "$@" --init_from $init --out $D > $D.log 2>&1
    rc=$?
    echo "$(date -u +%T) end S4 $st $arm try $try exit $rc" >> $P
    if [ $rc -eq 0 ] && [ -f $D/m.pt ] && [ -f $D/summary.json ]; then return 0; fi
    grep -q "OutOfMemoryError" $D.log || { echo "$(date -u +%T) FAIL S4 $st $arm (not out of memory)" >> $P; return 1; }
    sleep 120
  done
  echo "$(date -u +%T) FAIL S4 $st $arm (out of memory 3 times)" >> $P; return 1; }
while read st args; do
  for arm in old null new; do
    case $arm in old) code=A_old; init=/root/states/$st/init.pt;; null) code=A_old; init=$W/$st/null/init.pt;; new) code=A_a3f; init=/root/states/$st/init.pt;; esac
    while [ $(gpu_apps) -ge 3 ] || [ $(gpu_free) -lt 5000 ]; do sleep 15; done
    one $st $arm $code $init $args < /dev/null &
    sleep 60
  done
done < /root/s4_jobs.txt
wait
python /root/s4_eval.py --root $W --states $(awk '{print $1}' /root/s4_jobs.txt) --init_dir /root/states --json $W/s4.json > $W/s4.txt 2>&1
echo "$(date -u +%T) eval S4 exit $?" >> $P
touch $W/DONE
