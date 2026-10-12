#!/bin/bash
# Frequency test C′ (decision Chris 2026-10-12, PLAN 4.8): jobs/freq_K1pair.txt (15 jobs: the pair K1 (A 20 pJ/m) and
# K1A10 (A 10 pJ/m) at f/2 and 2 f, and K1 at f (f1), 3 seeds, accumulators with harmonics (7.72, 7.73), no snapshots). Starts only when the pair chain has finished
# (polls for runs/PAIR_DONE; chain_pair.sh writes it only after a queue exit code 0, else this script waits). Then it
# removes the pair chain from onstart, takes and keeps the pair lock, pins its own commit (/root/freq_commit), runs the
# queue and the post-processing of the accumulators. A stop file (runs/STOP_NEW_JOBS) stops the queue before any job.
# Pause (Chris 2026-10-12, PROTOKOLL 7.70): it also waits for runs/FREQ_GO, written by hand after the extended
# accumulators (option O2) are tested and /root/freq_commit points to that commit.
#   install: echo <commit> > /root/freq_commit; cp chain_freq.sh /root/; add the onstart line; start with nohup.
B=/root/magnum.np/studies/beta_bpeak
exec 3> /root/freq_chain.lock
flock -n 3 || { echo "$(date '+%F %T') freq chain already running"; exit 0; }
until [ -f $B/runs/PAIR_DONE ]; do sleep 300; done
until [ -f $B/runs/FREQ_GO ]; do sleep 300; done
sed -i "/chain_pair.sh/d" /root/onstart.sh
exec 4> /root/pair_chain.lock
flock 4                                  # waits for a pair chain that a reboot restarted; held until the end
[ -f /root/freq_commit ] || { echo "$(date '+%F %T') no /root/freq_commit: stop"; exit 1; }
C=$(cat /root/freq_commit)
cd /root/magnum.np || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git fetch -q --depth 50 origin study/beta-bpeak-pinning && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
cd $B || exit 1
L=queue_freq.log
echo "$(date '+%F %T') freq chain (re)start at $(git log --oneline -1)" >> $L
python run_queue.py jobs/freq_K1pair.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 30 \
    --status status_freq.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rc=$?; echo "$(date '+%F %T') freq queue exit code $rc" >> $L
touch runs/FREQ_QUEUE_DONE
for i in 1 2 3; do                       # a killed worker (BrokenProcessPool) ends the script: restart it
    python instance_postproc.py --runs runs --stop runs/FREQ_QUEUE_DONE --workers 4 \
        --pattern '^K1(A10)?_\d_f(half|1|2)$' >> $L 2>&1
    p=$?; echo "$(date '+%F %T') freq post-processing exit code $p" >> $L
    [ $p -eq 0 ] && break
done
[ $rc -eq 0 ] && touch runs/FREQ_DONE
