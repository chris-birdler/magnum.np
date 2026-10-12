#!/bin/bash
# Damping test at small B (decision Chris 2026-10-12, PROTOKOLL 7.66): jobs/cond_a001.txt (16 jobs: C0_a001 … C7_a001
# with the design history, then 30 cycles at 9 mT or at 20 mT; earlier stages capped at the design length of 9 cycles;
# start from C<k>/init.pt as the design runs). Starts only when the size chain has finished: it polls for
# runs/SIZE_DONE (chain_size.sh writes it only after a queue exit code 0; else this script waits). Then it removes the
# size chain from onstart, takes and keeps the size lock (it waits for a size chain that a reboot restarted), pins its
# own commit (/root/alpha_commit), runs the queue with the determinism gate (check_reruns.py: the first cycles must
# repeat C<k>_a001 bit for bit) and the post-processing of the accumulators. A stop file (runs/STOP_NEW_JOBS) stops
# the queue before any job starts; after its removal the chain needs a manual restart. Restart safe.
#   install: echo <commit> > /root/alpha_commit; cp chain_alpha.sh /root/; add the onstart line; start with nohup.
B=/root/magnum.np/studies/beta_bpeak
exec 5> /root/alpha_chain.lock
flock -n 5 || { echo "$(date '+%F %T') alpha chain already running"; exit 0; }
until [ -f $B/runs/SIZE_DONE ]; do sleep 300; done
sed -i "/chain_size.sh/d" /root/onstart.sh
exec 6> /root/size_chain.lock
flock 6                                  # waits for a size chain that a reboot restarted; held until the end
[ -f /root/alpha_commit ] || { echo "$(date '+%F %T') no /root/alpha_commit: stop"; exit 1; }
C=$(cat /root/alpha_commit)
cd /root/magnum.np || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git fetch -q --depth 50 origin study/beta-bpeak-pinning && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
cd $B || exit 1
L=queue_alpha.log
echo "$(date '+%F %T') alpha chain (re)start at $(git log --oneline -1)" >> $L
python check_reruns.py --watch 600 --stop_file runs/STOP_NEW_JOBS --until runs/ALPHA_QUEUE_DONE >> $L 2>&1 &
W=$!
python run_queue.py jobs/cond_a001.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 30 \
    --status status_alpha.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rc=$?; echo "$(date '+%F %T') alpha queue exit code $rc" >> $L
touch runs/ALPHA_QUEUE_DONE
wait $W; echo "$(date '+%F %T') determinism gate exit code $?" >> $L
for i in 1 2 3; do                       # a killed worker (BrokenProcessPool) ends the script: restart it
    python instance_postproc.py --runs runs --stop runs/ALPHA_QUEUE_DONE --workers 4 \
        --pattern '^cond_C\d_a001_b\d+$' >> $L 2>&1
    p=$?; echo "$(date '+%F %T') alpha post-processing exit code $p" >> $L
    [ $p -eq 0 ] && break
done
[ $rc -eq 0 ] && touch runs/ALPHA_DONE
