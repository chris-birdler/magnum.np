#!/bin/bash
# Mesh check at d/l_ex 212 (PLAN 4.6, option 2 of Chris 2026-10-11): jobs/mesh212.txt (24 jobs: base microstructure,
# dx 1.5 and dx 3, 9 / 50 / 100 mT, no accumulators, no snapshots). Starts only when the mc2b chain has finished:
# it polls for runs/MC2B_DONE (the last line of chain_mc2b.sh; it does not touch the mc2b lock before, so that a
# chain_mc2b.sh started by onstart after a reboot always gets its lock). An early exit of chain_mc2b.sh (no commit
# file, checkout or retire failed) writes no MC2B_DONE: then this script waits. After MC2B_DONE it removes the mc2b
# chain from onstart (else an instance restart would check out the old pinned commit), then takes and keeps the mc2b
# lock (it waits for an mc2b chain that a reboot restarted; a later mc2b start exits with "already running", so the
# two chains never run on two commits at the same time), pins its own commit (/root/size_commit) and runs the queue
# with a 30 h job limit. A stop file (runs/STOP_NEW_JOBS, e.g. from the mc2b determinism gate) stops this queue
# before any job starts; after its removal the chain needs a manual restart. Restart safe.
#   install: echo <commit> > /root/size_commit; cp chain_size.sh /root/; add the onstart line; start with nohup.
B=/root/magnum.np/studies/beta_bpeak
exec 7> /root/size_chain.lock
flock -n 7 || { echo "$(date '+%F %T') size chain already running"; exit 0; }
until [ -f $B/runs/MC2B_DONE ]; do sleep 300; done
sed -i "/chain_mc2b.sh/d" /root/onstart.sh
exec 8> /root/mc2b_chain.lock
flock 8                                  # waits for an mc2b chain that a reboot restarted; held until the end
[ -f /root/size_commit ] || { echo "$(date '+%F %T') no /root/size_commit: stop"; exit 1; }
C=$(cat /root/size_commit)
cd /root/magnum.np || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git fetch -q --depth 50 origin study/beta-bpeak-pinning && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
cd $B || exit 1
L=queue_size.log
echo "$(date '+%F %T') size chain (re)start at $(git log --oneline -1)" >> $L
python run_queue.py jobs/mesh212.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 30 \
    --status status_size.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rc=$?; echo "$(date '+%F %T') size queue exit code $rc" >> $L
[ $rc -eq 0 ] && touch runs/SIZE_DONE
