#!/bin/bash
# Option 1 of Chris with patch (2026-10-14, PROTOKOLL 7.101): tolerance twins of the main run at 30 MHz (jobs/atol30.txt,
# 3 jobs: the production line of K1_k with atol 1e-6 and an own relaxation). Waits until the atol60 queue has ended
# (runs/ATOL60_QUEUE_DONE) and the atol60 chain has released its lock, then pins its own commit (/root/atol30_commit,
# read after the wait), runs the queue, post-processes and writes the pre-registered report runs/atol30_report.md
# (eval_atol30.py). A stop file (runs/STOP_NEW_JOBS) stops the queue before any job.
#   install: echo <commit> > /root/atol30_commit; cp chain_atol30.sh /root/; add the onstart line; start with nohup.
B=/root/magnum.np/studies/beta_bpeak
exec 3> /root/atol30_chain.lock
flock -n 3 || { echo "$(date '+%F %T') atol30 chain already running"; exit 0; }
[ -f /root/atol30_commit ] || { echo "$(date '+%F %T') no /root/atol30_commit: stop"; exit 1; }
[ -f $B/runs/ATOL30_DONE ] && { echo "$(date '+%F %T') atol30 already done"; exit 0; }
sleep 180                                   # at a restart the atol60 chain takes its lock first
echo "$(date '+%F %T') atol30 chain waits for the end of the atol60 queue and chain"
while ! { [ -f $B/runs/ATOL60_QUEUE_DONE ] && flock -n /root/atol60_chain.lock true; }; do
    sleep 300                               # atol60 not done or its chain still running (the lock is only tested)
done
C=$(cat /root/atol30_commit)
cd /root/magnum.np || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git fetch -q --depth 50 origin study/beta-bpeak-pinning && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
cd $B || exit 1
L=queue_atol30.log
echo "$(date '+%F %T') atol30 chain start at $(git log --oneline -1)" >> $L
rm -f runs/ATOL30_QUEUE_DONE
python run_queue.py jobs/atol30.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 10 --gate_maxiter 0 \
    --status status_atol30.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rc=$?; echo "$(date '+%F %T') atol30 queue exit code $rc" >> $L
touch runs/ATOL30_QUEUE_DONE
for i in 1 2 3; do
    python instance_postproc.py --runs runs --stop runs/ATOL30_QUEUE_DONE --workers 4 --pattern '^K1_\d_L6$' >> $L 2>&1
    p=$?; echo "$(date '+%F %T') atol30 post-processing exit code $p" >> $L
    [ $p -eq 0 ] && break
done
python eval_atol30.py > runs/atol30_report.md 2>> $L
[ $rc -eq 0 ] && touch runs/ATOL30_DONE
