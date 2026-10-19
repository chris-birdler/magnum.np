#!/bin/bash
B=/root/magnum.np/studies/beta_bpeak
while [ ! -f $B/runs/ALL3_DONE ]; do sleep 60; done
cd /root/magnum.np && git fetch -q --depth 1 origin study/beta-bpeak-pinning && git checkout -q -f -B study/beta-bpeak-pinning FETCH_HEAD
cd $B && git log --oneline -1 > queue_alpha2.log && ls jobs/alphatest2.txt >> queue_alpha2.log
python run_queue.py jobs/alphatest2.txt --gpus 0,1,2 >> queue_alpha2.log 2>&1
touch runs/ALL4_DONE
