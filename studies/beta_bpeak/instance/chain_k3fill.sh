#!/bin/bash
# Fine-mesh fill-in at K3 (go Chris 2026-10-10, PROTOKOLL 7.51, option 2): jobs/meshK3fill.txt (36 jobs: 35 / 20 / 70 mT,
# dx 1.5 first, then dx 3, 6 seeds, start from the init.pt of the K3 9 mT jobs, with --acc). Starts only after the K3
# chain has ended (it waits for the lock of chain_k3.sh), then pins its own commit (/root/k3fill_commit), removes the
# K3 chain from onstart (else an instance restart would check out the old pinned commit), runs the queue with a 30 h
# job limit and then the post-processing of the accumulators (acc_reduce, md5, delete acc_*.npz). Restart safe like
# chain_k3.sh.
#   install: echo <commit> > /root/k3fill_commit; cp chain_k3fill.sh /root/; add the onstart line; start with nohup.
B=/root/magnum.np/studies/beta_bpeak
exec 7> /root/k3fill_chain.lock
flock -n 7 || { echo "$(date '+%F %T') k3fill chain already running"; exit 0; }
exec 8> /root/k3_chain.lock
flock 8                                  # blocks until the K3 chain has ended
flock -u 8
sed -i "/chain_k3.sh/d" /root/onstart.sh
[ -f /root/k3fill_commit ] || { echo "$(date '+%F %T') no /root/k3fill_commit: stop"; exit 1; }
C=$(cat /root/k3fill_commit)
cd /root/magnum.np || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git fetch -q --depth 50 origin study/beta-bpeak-pinning && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
cd $B || exit 1
echo "$(date '+%F %T') k3fill chain (re)start at $(git log --oneline -1)" >> queue_k3fill.log
python run_queue.py jobs/meshK3fill.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 30 \
    --status status_k3fill.json >> queue_k3fill.log 2>&1
rc=$?; echo "$(date '+%F %T') k3fill queue exit code $rc" >> queue_k3fill.log
touch runs/K3FILL_QUEUE_DONE
for i in 1 2 3; do                       # a killed worker (BrokenProcessPool) ends the script: restart it
    python instance_postproc.py --runs runs --stop runs/K3FILL_QUEUE_DONE --workers 4 \
        --pattern '^MK3_dx(15|3)_s\d+_b(20|35|70)$' >> queue_k3fill.log 2>&1
    rc=$?; echo "$(date '+%F %T') k3fill post-processing exit code $rc" >> queue_k3fill.log
    [ $rc -eq 0 ] && break
done
touch runs/K3FILL_DONE
