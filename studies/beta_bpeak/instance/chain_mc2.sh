#!/bin/bash
# PLAN 4.5 (go Chris 2026-10-11): reruns with --min_calm 2 and the conditioning test. Starts only after the K3 fill-in
# chain has ended (it waits for the lock of chain_k3fill.sh), then pins its own commit (/root/mc2_commit), removes the
# fill-in chain from onstart, makes the reruns of the fill-in jobs with a calm tail < 2 (cap 15 V100-h), runs one queue
# (jobs/rerun_mc2.txt, the fill-in reruns, jobs/cond.txt) with the determinism gate (check_reruns.py: a mismatch
# creates runs/STOP_NEW_JOBS, no new job starts), then the post-processing of the accumulators. Restart safe like
# chain_k3.sh.
#   install: echo <commit> > /root/mc2_commit; cp chain_mc2.sh /root/; add the onstart line; start with nohup.
B=/root/magnum.np/studies/beta_bpeak
exec 6> /root/mc2_chain.lock
flock -n 6 || { echo "$(date '+%F %T') mc2 chain already running"; exit 0; }
exec 7> /root/k3fill_chain.lock
flock 7                                  # blocks until the fill-in chain (with its post-processing) has ended
flock -u 7
sed -i "/chain_k3fill.sh/d" /root/onstart.sh
[ -f /root/mc2_commit ] || { echo "$(date '+%F %T') no /root/mc2_commit: stop"; exit 1; }
C=$(cat /root/mc2_commit)
cd /root/magnum.np || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git fetch -q --depth 50 origin study/beta-bpeak-pinning && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
cd $B || exit 1
L=queue_mc2.log
echo "$(date '+%F %T') mc2 chain (re)start at $(git log --oneline -1)" >> $L
python rerun_jobs.py --reruns jobs/meshK3fill.txt:12 --cap_vh 15 --out jobs/rerun_mc2_fill.txt >> $L 2>&1
cat jobs/rerun_mc2.txt jobs/rerun_mc2_fill.txt jobs/cond.txt > jobs/queue_mc2.txt
echo "$(date '+%F %T') queue: $(wc -l < jobs/queue_mc2.txt) jobs" >> $L
python check_reruns.py --watch 600 --stop_file runs/STOP_NEW_JOBS --until runs/MC2_QUEUE_DONE >> $L 2>&1 &
W=$!
python run_queue.py jobs/queue_mc2.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 30 \
    --status status_mc2.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rc=$?; echo "$(date '+%F %T') mc2 queue exit code $rc" >> $L
touch runs/MC2_QUEUE_DONE
wait $W; echo "$(date '+%F %T') determinism gate exit code $?" >> $L
for i in 1 2 3; do                       # a killed worker (BrokenProcessPool) ends the script: restart it
    python instance_postproc.py --runs runs --stop runs/MC2_QUEUE_DONE --workers 4 \
        --pattern '^(.+_mc2|cond_C\d_b\d+)$' >> $L 2>&1
    rc=$?; echo "$(date '+%F %T') mc2 post-processing exit code $rc" >> $L
    [ $rc -eq 0 ] && break
done
touch runs/MC2_DONE
