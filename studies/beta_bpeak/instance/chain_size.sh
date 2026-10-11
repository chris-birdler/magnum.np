#!/bin/bash
# Mesh check at d/l_ex 212 (PLAN 4.6, option 2 of Chris 2026-10-11): jobs/mesh212.txt (24 jobs: base microstructure,
# dx 1.5 and dx 3, 9 / 50 / 100 mT, no accumulators, no snapshots). Starts only when the mc2b chain has finished
# (runs/MC2B_DONE exists): it waits for the lock of chain_mc2b.sh and, if MC2B_DONE is missing (mc2b ended with an
# error, or this script got the lock first after a reboot), waits again. Then it removes the mc2b chain from onstart
# (else an instance restart would check out the old pinned commit), pins its own commit (/root/size_commit) and runs
# the queue with a 30 h job limit. A stop file (runs/STOP_NEW_JOBS, e.g. from the mc2b determinism gate) stops this
# queue before any job starts; after its removal the chain needs a manual restart. Restart safe.
#   install: echo <commit> > /root/size_commit; cp chain_size.sh /root/; add the onstart line; start with nohup.
B=/root/magnum.np/studies/beta_bpeak
exec 7> /root/size_chain.lock
flock -n 7 || { echo "$(date '+%F %T') size chain already running"; exit 0; }
exec 8> /root/mc2b_chain.lock
until [ -f $B/runs/MC2B_DONE ]; do
    flock 8                              # blocks while the mc2b chain runs
    flock -u 8
    [ -f $B/runs/MC2B_DONE ] || sleep 300
done
sed -i "/chain_mc2b.sh/d" /root/onstart.sh
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
