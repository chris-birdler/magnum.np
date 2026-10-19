#!/bin/bash
B=/root/magnum.np/studies/beta_bpeak
while [ ! -f $B/runs/ALL4_DONE ]; do sleep 60; done
cd /root/magnum.np && git fetch -q --depth 1 origin study/beta-bpeak-pinning && git checkout -q -f -B study/beta-bpeak-pinning FETCH_HEAD
cd $B && git log --oneline -1 > queue_snap.log && ls jobs/snaptest.txt >> queue_snap.log
python run_queue.py jobs/snaptest.txt --gpus 0,1 >> queue_snap.log 2>&1
touch runs/ALL5_DONE
