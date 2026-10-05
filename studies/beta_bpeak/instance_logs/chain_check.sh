#!/bin/bash
# checks before the Sobol design (PROTOKOLL 7.24): freq9 + phase16, after mesh300
B=/root/magnum.np/studies/beta_bpeak
while [ ! -f $B/runs/ALL6_DONE ]; do sleep 60; done
cd /root/magnum.np && git fetch -q --depth 1 origin study/beta-bpeak-pinning && git checkout -q -f -B study/beta-bpeak-pinning FETCH_HEAD
cd $B && git log --oneline -1 > queue_check.log && cat jobs/phase16.txt jobs/freq9.txt > jobs/check_combined.txt && ls jobs/check_combined.txt >> queue_check.log
python run_queue.py jobs/check_combined.txt --gpus 0,1,2,3 >> queue_check.log 2>&1
touch runs/ALL7_DONE
