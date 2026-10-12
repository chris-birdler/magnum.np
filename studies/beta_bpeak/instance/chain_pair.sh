#!/bin/bash
# Paired field images for the A effect (decision Chris 2026-10-12): jobs/pair_K1A10.txt (3 jobs: corner K1 with
# A 10 pJ/m instead of 20, same seeds, snapshots at 9 / 50 / 150 mT). Starts only when the alpha chain has finished
# (polls for runs/ALPHA_DONE; chain_alpha.sh writes it only after a queue exit code 0, else this script waits). Then it
# removes the alpha chain from onstart, takes and keeps the alpha lock, pins its own commit (/root/pair_commit), runs
# the queue and the post-processing of the snapshots and accumulators (domains, phase maps, cuts; the VTI files are
# deleted after the analysis, standing rule). A stop file (runs/STOP_NEW_JOBS) stops the queue before any job starts.
#   install: echo <commit> > /root/pair_commit; cp chain_pair.sh /root/; add the onstart line; start with nohup.
B=/root/magnum.np/studies/beta_bpeak
exec 3> /root/pair_chain.lock
flock -n 3 || { echo "$(date '+%F %T') pair chain already running"; exit 0; }
until [ -f $B/runs/ALPHA_DONE ]; do sleep 300; done
sed -i "/chain_alpha.sh/d" /root/onstart.sh
exec 4> /root/alpha_chain.lock
flock 4                                  # waits for an alpha chain that a reboot restarted; held until the end
[ -f /root/pair_commit ] || { echo "$(date '+%F %T') no /root/pair_commit: stop"; exit 1; }
C=$(cat /root/pair_commit)
cd /root/magnum.np || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git fetch -q --depth 50 origin study/beta-bpeak-pinning && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
cd $B || exit 1
L=queue_pair.log
echo "$(date '+%F %T') pair chain (re)start at $(git log --oneline -1)" >> $L
python run_queue.py jobs/pair_K1A10.txt --gpus 0,1,2 --min_free_gb 3 --retries 1 --job_timeout_h 30 \
    --status status_pair.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rc=$?; echo "$(date '+%F %T') pair queue exit code $rc" >> $L
touch runs/PAIR_QUEUE_DONE
for i in 1 2 3; do                       # a killed worker (BrokenProcessPool) ends the script: restart it
    python instance_postproc.py --runs runs --stop runs/PAIR_QUEUE_DONE --workers 3 --pattern '^K1A10_\d$' >> $L 2>&1
    p=$?; echo "$(date '+%F %T') pair post-processing exit code $p" >> $L
    [ $p -eq 0 ] && break
done
[ $rc -eq 0 ] && touch runs/PAIR_DONE
