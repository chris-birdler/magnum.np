#!/bin/bash
# Option 2 of Chris (2026-10-14, PROTOKOLL 7.100): long tolerance test at 60 MHz (jobs/atol60.txt, 6 jobs, 40 cycles,
# atol 1e-5 and 1e-6). Waits until B0b has its report (runs/conv_b0_report.md) and the sweep/B0b chain has ended,
# then pins its own commit (/root/atol60_commit), runs the queue, post-processes and writes the pre-registered report
# runs/atol60_report.md (eval_atol60.py). A stop file (runs/STOP_NEW_JOBS) stops the queue before any job.
#   install: echo <commit> > /root/atol60_commit; cp chain_atol60.sh /root/; add the onstart line; start with nohup.
B=/root/magnum.np/studies/beta_bpeak
exec 3> /root/atol60_chain.lock
flock -n 3 || { echo "$(date '+%F %T') atol60 chain already running"; exit 0; }
[ -f /root/atol60_commit ] || { echo "$(date '+%F %T') no /root/atol60_commit: stop"; exit 1; }
C=$(cat /root/atol60_commit)
[ -f $B/runs/ATOL60_DONE ] && { echo "$(date '+%F %T') atol60 already done"; exit 0; }
sleep 120                                   # at a restart the sweep/B0b chain takes its lock first
while ! { [ -f $B/runs/conv_b0_report.md ] && flock -n /root/sweep_chain.lock true; }; do
    sleep 300                               # B0b not done or its chain still running (the lock is only tested)
done
cd /root/magnum.np || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git fetch -q --depth 50 origin study/beta-bpeak-pinning && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
cd $B || exit 1
L=queue_atol60.log
echo "$(date '+%F %T') atol60 chain start at $(git log --oneline -1)" >> $L
python run_queue.py jobs/atol60.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 10 --gate_maxiter 0 \
    --status status_atol60.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rc=$?; echo "$(date '+%F %T') atol60 queue exit code $rc" >> $L
touch runs/ATOL60_QUEUE_DONE
for i in 1 2 3; do
    python instance_postproc.py --runs runs --stop runs/ATOL60_QUEUE_DONE --workers 4 --pattern '^K1_\d_sw060L[56]$' >> $L 2>&1
    p=$?; echo "$(date '+%F %T') atol60 post-processing exit code $p" >> $L
    [ $p -eq 0 ] && break
done
python eval_atol60.py > runs/atol60_report.md 2>> $L
[ $rc -eq 0 ] && touch runs/ATOL60_DONE
