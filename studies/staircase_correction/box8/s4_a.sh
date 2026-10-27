#!/bin/bash
# box 8: stage 4 of PLAN 11 for the A bundle (old main 9e0ea7c, new perf/timestep 3d9de71), 8 phase-1 states, 4 cycles
export CUDA_DEVICE=0
P=/root/prof/progress.txt
until [ -f /root/bs3/DONE ]; do sleep 30; done
W=/root/s4a; mkdir -p $W
while read st args; do
  mkdir -p $W/$st
  PYTHONPATH=/root/A_old python /root/perturb_init.py /root/states/$st/init.pt $W/$st/null/init.pt --seed 1 > $W/$st/perturb.txt 2>&1
done < /root/s4_jobs.txt
one() { st=$1; arm=$2; code=$3; init=$4; shift 4; D=$W/$st/$arm; mkdir -p $D; echo "$(date -u +%T) start S4 $st $arm" >> $P; cd /root/$code/studies/beta_bpeak; PYTHONPATH=/root/$code python /root/step_dump.py --at 1024 --save $D/m.pt --go_on -- "$@" --init_from $init --out $D > $D.log 2>&1; echo "$(date -u +%T) end S4 $st $arm exit $?" >> $P; }
while read st args; do
  for arm in old null new; do
    case $arm in old) code=A_old; init=/root/states/$st/init.pt;; null) code=A_old; init=$W/$st/null/init.pt;; new) code=A_a3f; init=/root/states/$st/init.pt;; esac
    while [ $(jobs -rp | wc -l) -ge 3 ]; do sleep 15; done
    one $st $arm $code $init $args &
    sleep 60
  done
done < /root/s4_jobs.txt
wait
python /root/s4_eval.py --root $W --states $(awk '{print $1}' /root/s4_jobs.txt) --init_dir /root/states --json $W/s4.json > $W/s4.txt 2>&1
touch $W/DONE
