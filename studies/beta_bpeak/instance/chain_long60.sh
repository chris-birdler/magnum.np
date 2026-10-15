#!/bin/bash
# Point 5 of the open questions (decision Chris 2026-10-15, PROTOKOLL 7.103): long runs at 60 MHz after the C′ history
# (jobs/long60.txt, 9 jobs: large particle at 9 mT, small particle at 20 and 35 mT, 40 cycles each). No other chain
# runs; the B0b chain was removed from onstart (it has no done marker check and would check out another commit). Pins
# its own commit (/root/long60_commit), runs the queue, post-processes and writes the pre-registered report
# runs/long60_report.md (eval_long60.py). A stop file (runs/STOP_NEW_JOBS) stops the queue before any job.
#   install: echo <commit> > /root/long60_commit; cp chain_long60.sh /root/; add the onstart line; start with nohup.
B=/root/magnum.np/studies/beta_bpeak
exec 3> /root/long60_chain.lock
flock -n 3 || { echo "$(date '+%F %T') long60 chain already running"; exit 0; }
[ -f /root/long60_commit ] || { echo "$(date '+%F %T') no /root/long60_commit: stop"; exit 1; }
[ -f $B/runs/LONG60_DONE ] && { echo "$(date '+%F %T') long60 already done"; exit 0; }
sleep 60                                    # at a restart: give the system time to settle
free=$(df --output=avail -BG /root | tail -1 | tr -dc 0-9)
[ "$free" -ge 6 ] || { echo "$(date '+%F %T') only ${free} GB free: stop"; exit 1; }
C=$(cat /root/long60_commit)
cd /root/magnum.np || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git fetch -q --depth 50 origin study/beta-bpeak-pinning && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
cd $B || exit 1
L=queue_long60.log
echo "$(date '+%F %T') long60 chain start at $(git log --oneline -1)" >> $L
rm -f runs/LONG60_QUEUE_DONE
python run_queue.py jobs/long60.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 6 --gate_maxiter 0 \
    --status status_long60.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rc=$?; echo "$(date '+%F %T') long60 queue exit code $rc" >> $L
touch runs/LONG60_QUEUE_DONE
for i in 1 2 3; do
    python instance_postproc.py --runs runs --stop runs/LONG60_QUEUE_DONE --workers 4 --pattern '^cond_K1(A10)?_\d_f2_b\d+$' >> $L 2>&1
    p=$?; echo "$(date '+%F %T') long60 post-processing exit code $p" >> $L
    [ $p -eq 0 ] && break
done
python eval_long60.py > runs/long60_report.md 2>> $L
e=$?; echo "$(date '+%F %T') long60 evaluation exit code $e" >> $L
[ $rc -eq 0 ] && [ $e -eq 0 ] && touch runs/LONG60_DONE
