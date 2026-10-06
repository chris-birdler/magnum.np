#!/bin/bash
# Sobol design v2 (PLAN 4.0) on the vast.ai instance. Start only after Chris's go and enough credit.
# main (96 + 8) and corners (12) in random order, then the 8 centre replicates at alpha 0.01 (from the init.pt
# of the alpha 0.1 centre runs, with snapshots); post-processing of the snapshots in parallel
# (instance_postproc.py, standing go). The mesh check at the corners K2/K3 is NOT part of this chain
# (after the main run, only if needed, separate go).
#
# Restart safety (audit 3, item 2):
# * the code is pinned to the commit in /root/sobol_commit (written once at installation). A restart checks out
#   exactly this commit and never a newer one: run_loops resumes an interrupted run only with the same commit.
#   Without /root/sobol_commit the chain stops.
# * a lock (flock) keeps a second chain from starting; /root/onstart.sh starts this script after an instance
#   restart (finished runs are skipped by run_queue, interrupted runs resume from their checkpoint).
# * all logs are appended, never overwritten; the shuffled job file is made once and then kept.
# Unattended operation (audit 3, item 4): retry once, GPU paused / disabled after failure bursts, job timeout 12 h,
# status_sobol.json every 5 min, alerts in runs/ALERT, post-processing restarted if it dies.
#
#   install: echo <commit> > /root/sobol_commit; rm -f <B>/jobs/sobol_a01_shuffled.txt; cp chain_sobol.sh /root/;
#            append "nohup bash /root/chain_sobol.sh >> /root/chain_sobol.out 2>&1 &" to /root/onstart.sh;
#            nohup bash /root/chain_sobol.sh >> /root/chain_sobol.out 2>&1 &
B=/root/magnum.np/studies/beta_bpeak
exec 9> /root/sobol_chain.lock
flock -n 9 || { echo "$(date '+%F %T') chain already running"; exit 0; }
[ -f /root/sobol_commit ] || { echo "$(date '+%F %T') no /root/sobol_commit: stop"; exit 1; }
C=$(cat /root/sobol_commit)
cd /root/magnum.np || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git fetch -q --depth 50 origin study/beta-bpeak-pinning && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
cd $B || exit 1
echo "$(date '+%F %T') chain (re)start at $(git log --oneline -1)" >> queue_sobol.log
[ -f jobs/sobol_a01_shuffled.txt ] || python -c "import random; L=open('jobs/sobol_main.txt').readlines()+open('jobs/sobol_corners.txt').readlines(); random.Random(20261006).shuffle(L); open('jobs/sobol_a01_shuffled.txt','w').writelines(L)"
# post-processing under a supervisor loop (audit 3, item 4): restarted if it dies, until ALL10_DONE and no run left
( while true; do
      python instance_postproc.py --runs runs --stop runs/ALL10_DONE --workers 4 >> postproc.log 2>&1 && break
      echo "$(date '+%F %T') instance_postproc exited with an error: restart in 60 s" >> postproc.log
      echo "$(date '+%F %T') instance_postproc restarted" >> runs/ALERT
      sleep 60
  done ) &
Q="--gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 12 --status status_sobol.json"
python run_queue.py jobs/sobol_a01_shuffled.txt $Q >> queue_sobol.log 2>&1
rc=$?; echo "$(date '+%F %T') main queue exit code $rc" >> queue_sobol.log
python run_queue.py jobs/sobol_a001.txt $Q >> queue_sobol.log 2>&1
rc=$?; echo "$(date '+%F %T') alpha 0.01 queue exit code $rc" >> queue_sobol.log
touch runs/ALL10_DONE
echo "$(date '+%F %T') queues finished" >> queue_sobol.log
wait
