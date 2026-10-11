#!/bin/bash
# PLAN 4.5 after the gate stop of 2026-10-11 (PROTOKOLL 7.59, go Chris): reruns with the start path of their original
# (a rerun relaxes in the process when its original did; only the reruns of jobs that started from init.pt load it).
# Waits for the first mc2 chain (its running jobs and post-processing), pins /root/mc2b_commit, renames the init.pt
# variants X_mc2 -> X_mc2init (rerun_jobs.py --retire; kept as a sensitivity record), removes runs/STOP_NEW_JOBS,
# then runs the queue with the determinism gate, then the post-processing of the accumulators. Restart safe.
#   install: echo <commit> > /root/mc2b_commit; cp chain_mc2b.sh /root/; add the onstart line; start with nohup.
B=/root/magnum.np/studies/beta_bpeak
exec 5> /root/mc2b_chain.lock
flock -n 5 || { echo "$(date '+%F %T') mc2b chain already running"; exit 0; }
exec 6> /root/mc2_chain.lock
flock 6                                  # blocks until the first mc2 chain has ended
flock -u 6
sed -i "/chain_mc2.sh/d" /root/onstart.sh
[ -f /root/mc2b_commit ] || { echo "$(date '+%F %T') no /root/mc2b_commit: stop"; exit 1; }
C=$(cat /root/mc2b_commit)
cd /root/magnum.np || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git fetch -q --depth 50 origin study/beta-bpeak-pinning && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
cd $B || exit 1
L=queue_mc2b.log
echo "$(date '+%F %T') mc2b chain (re)start at $(git log --oneline -1)" >> $L
python rerun_jobs.py --reruns jobs/meshK3fill.txt:12 --cap_vh 15 --out jobs/rerun_mc2_fill.txt >> $L 2>&1
cat jobs/rerun_mc2.txt jobs/rerun_mc2_fill.txt jobs/cond.txt > jobs/queue_mc2b.txt
python rerun_jobs.py --retire jobs/queue_mc2b.txt >> $L 2>&1 || { echo "$(date '+%F %T') retire failed: stop" >> $L; exit 1; }
if [ ! -f runs/MC2B_STARTED ]; then      # the stop file of the first chain; a new stop file after this start stays
    rm -f runs/STOP_NEW_JOBS
    touch runs/MC2B_STARTED
fi
echo "$(date '+%F %T') queue: $(wc -l < jobs/queue_mc2b.txt) jobs" >> $L
python check_reruns.py --watch 600 --stop_file runs/STOP_NEW_JOBS --until runs/MC2B_QUEUE_DONE >> $L 2>&1 &
W=$!
python run_queue.py jobs/queue_mc2b.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 30 \
    --status status_mc2b.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rc=$?; echo "$(date '+%F %T') mc2b queue exit code $rc" >> $L
touch runs/MC2B_QUEUE_DONE
wait $W; echo "$(date '+%F %T') determinism gate exit code $?" >> $L
for i in 1 2 3; do                       # a killed worker (BrokenProcessPool) ends the script: restart it
    python instance_postproc.py --runs runs --stop runs/MC2B_QUEUE_DONE --workers 4 \
        --pattern '^(.+_mc2|cond_C\d_b\d+)$' >> $L 2>&1
    rc=$?; echo "$(date '+%F %T') mc2b post-processing exit code $rc" >> $L
    [ $rc -eq 0 ] && break
done
touch runs/MC2B_DONE
