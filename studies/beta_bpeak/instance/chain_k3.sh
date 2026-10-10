#!/bin/bash
# Mesh check at the design corner K3 (PLAN 4.0 rule A, go Chris 2026-10-10): jobs/meshK3.txt (36 jobs, mesh300
# protocol, dx 1.5 and dx 3, 6 seeds). Starts only after the Sobol chain has ended completely (it waits for the lock
# of chain_sobol.sh), then pins its own commit (/root/k3_commit), removes the Sobol chain from onstart (else an
# instance restart would check out the old pinned commit) and runs the queue with a 30 h job limit (a dx 1.5 job at
# K3 needs several hours; mesh300 at d/l_ex 300 needed up to ~21 h). Restart safe like chain_sobol.sh.
#   install: echo <commit> > /root/k3_commit; cp chain_k3.sh /root/; add the onstart line; start with nohup.
B=/root/magnum.np/studies/beta_bpeak
exec 8> /root/k3_chain.lock
flock -n 8 || { echo "$(date '+%F %T') k3 chain already running"; exit 0; }
exec 9> /root/sobol_chain.lock
flock 9                                  # blocks until the Sobol chain (and its post-processing) has ended
flock -u 9
sed -i "/chain_sobol/d" /root/onstart.sh
[ -f /root/k3_commit ] || { echo "$(date '+%F %T') no /root/k3_commit: stop"; exit 1; }
C=$(cat /root/k3_commit)
cd /root/magnum.np || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git fetch -q --depth 50 origin study/beta-bpeak-pinning && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
cd $B || exit 1
echo "$(date '+%F %T') k3 chain (re)start at $(git log --oneline -1)" >> queue_k3.log
python run_queue.py jobs/meshK3.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 30 \
    --status status_k3.json >> queue_k3.log 2>&1
rc=$?; echo "$(date '+%F %T') k3 queue exit code $rc" >> queue_k3.log
touch runs/K3MESH_DONE
