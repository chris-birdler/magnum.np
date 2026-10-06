#!/bin/bash
# Sobol design v2 (PLAN 4.0) on the vast.ai instance. Start only after Chris's go and enough credit.
# main (96 + 8) and corners (12) in random order, then the alpha 0.01 subset (32 + 8, from the init.pt of
# the alpha 0.1 runs); post-processing of the snapshots in parallel (instance_postproc.py, standing go).
B=/root/magnum.np/studies/beta_bpeak
while [ ! -f $B/runs/ALL9_DONE ]; do sleep 60; done
cd /root/magnum.np && git fetch -q --depth 1 origin study/beta-bpeak-pinning && git checkout -q -f -B study/beta-bpeak-pinning FETCH_HEAD
cd $B && git log --oneline -1 > queue_sobol.log
python -c "import random; L=open('jobs/sobol_main.txt').readlines()+open('jobs/sobol_corners.txt').readlines(); random.Random(20261006).shuffle(L); open('jobs/sobol_a01_shuffled.txt','w').writelines(L)"
nohup python instance_postproc.py --runs runs --stop runs/ALL10_DONE --workers 4 >> postproc.log 2>&1 &
python run_queue.py jobs/sobol_a01_shuffled.txt --gpus 0,1,2,3 --min_free_gb 8 >> queue_sobol.log 2>&1
python run_queue.py jobs/sobol_a001.txt --gpus 0,1,2,3 --min_free_gb 8 >> queue_sobol.log 2>&1
touch runs/ALL10_DONE
wait
