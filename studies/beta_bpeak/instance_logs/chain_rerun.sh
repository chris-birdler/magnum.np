#!/bin/bash
B=/root/magnum.np/studies/beta_bpeak
while [ ! -f $B/runs/ALL5_DONE ]; do sleep 60; done
cd /root/magnum.np && git fetch -q --depth 1 origin study/beta-bpeak-pinning && git checkout -q -f -B study/beta-bpeak-pinning FETCH_HEAD
cd $B && git log --oneline -1 > queue_rerun.log && ls jobs/baseline.txt jobs/mesh300.txt >> queue_rerun.log
python run_queue.py jobs/baseline.txt --gpus 0,1,2,3 >> queue_rerun.log 2>&1
python run_queue.py jobs/mesh300.txt --gpus 0,1,2,3 >> queue_rerun.log 2>&1
touch runs/ALL6_DONE
