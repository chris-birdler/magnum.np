#!/bin/bash
# Option B of Chris (2026-10-14, PROTOKOLL 7.89): first the frequency sweep at 9 mT (jobs/sweep9.txt, 32 short jobs),
# then the B0 rerun with atol 1e-6 (jobs/conv_b0.txt, names B0b_*). Waits until the precision test queue (prec48) has
# ended, pins its commit (/root/sweep_commit, checked at the start), runs both queues and writes both pre-registered
# reports (runs/sweep9_report.md, runs/conv_b0_report.md). A stop file (runs/STOP_NEW_JOBS) stops a queue before any job.
#   install: echo <commit> > /root/sweep_commit; cp chain_sweep_b0.sh /root/; add the onstart line; start with nohup.
B=/root/magnum.np/studies/beta_bpeak
exec 3> /root/sweep_chain.lock
flock -n 3 || { echo "$(date '+%F %T') sweep chain already running"; exit 0; }
[ -f /root/sweep_commit ] || { echo "$(date '+%F %T') no /root/sweep_commit: stop"; exit 1; }
C=$(cat /root/sweep_commit)
sleep 60
while pgrep -f "jobs/prec48.txt" > /dev/null; do sleep 120; done
cd /root/magnum.np || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git fetch -q --depth 50 origin study/beta-bpeak-pinning && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
cd $B || exit 1
L=queue_sweep.log
echo "$(date '+%F %T') sweep chain (re)start at $(git log --oneline -1)" >> $L
python run_queue.py jobs/sweep9.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 10 --gate_maxiter 0 \
    --status status_sweep.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rc=$?; echo "$(date '+%F %T') sweep queue exit code $rc" >> $L
touch runs/SWEEP_QUEUE_DONE
for i in 1 2 3; do
    python instance_postproc.py --runs runs --stop runs/SWEEP_QUEUE_DONE --workers 4 --pattern '^K1(A10)?_\d_sw\d{3}$' >> $L 2>&1
    p=$?; echo "$(date '+%F %T') sweep post-processing exit code $p" >> $L
    [ $p -eq 0 ] && break
done
python eval_sweep.py > runs/sweep9_report.md 2>> $L
[ $rc -eq 0 ] && touch runs/SWEEP_DONE
L2=queue_conv.log
echo "$(date '+%F %T') B0b (atol 1e-6) start at $(git log --oneline -1)" >> $L2
python run_queue.py jobs/conv_b0.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 48 --gate_maxiter 0 \
    --status status_conv.json --stop_file runs/STOP_NEW_JOBS >> $L2 2>&1
rc2=$?; echo "$(date '+%F %T') B0b queue exit code $rc2" >> $L2
touch runs/CONV_QUEUE_DONE
python eval_conv.py > runs/conv_b0_report.md 2>> $L2
[ $rc2 -eq 0 ] && touch runs/CONV_DONE
