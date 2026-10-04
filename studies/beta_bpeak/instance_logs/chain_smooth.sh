#!/bin/bash
B=/root/magnum.np/studies/beta_bpeak
while [ ! -f $B/runs/ALL_DONE ]; do sleep 60; done
cd /root/magnum.np && git fetch -q --depth 1 origin study/beta-bpeak-pinning && git checkout -q -f -B study/beta-bpeak-pinning FETCH_HEAD
cd $B && git log --oneline -1 > queue_smooth.log && ls jobs/smoothtest.txt >> queue_smooth.log
python run_queue.py jobs/smoothtest.txt --gpus 0,1,2,3 >> queue_smooth.log 2>&1
touch runs/ALL2_DONE
