#!/bin/bash
# Frequency test at the corner K3 (decision Chris 2026-10-15, PROTOKOLL 7.106): jobs/k3freq_$1.txt (pair_jobs.py
# --k3freq [newell]; 18 jobs: K3 at 30 … 170 MHz, 3 seeds; $1 = fd or newell). Starts only after the go of Chris for this
# variant (runs/K3FREQ_$1_GO, written by hand) and after no other chain uses the GPUs (it waits while chain_newell.sh
# runs). Pins its own commit on feature/newell-pbc (/root/k3freq_commit). The queue keeps the relaxation gate: a run
# whose virgin relaxation does not converge stops (exit code 3) and runs once more with relax_maxiter 200000. Logs every
# exit code; K3FREQ_$1_DONE only if the queue and the report end with code 0. Stop file runs/STOP_NEW_JOBS.
#   install: echo <commit> > /root/k3freq_commit; cp chain_k3freq.sh /root/; add "nohup bash /root/chain_k3freq.sh fd
#   >> /root/k3freq_chain.out 2>&1 &" to onstart; start it the same way; touch runs/K3FREQ_fd_GO.
V=$1
[ "$V" = fd ] || [ "$V" = newell ] || { echo "$(date '+%F %T') variant must be fd or newell"; exit 1; }
P=K3; [ "$V" = newell ] && P=K3N
B=/root/magnum.np/studies/beta_bpeak
exec 3> /root/k3freq_chain_$V.lock
flock -n 3 || { echo "$(date '+%F %T') k3freq $V chain already running"; exit 0; }
[ -f /root/k3freq_commit ] || { echo "$(date '+%F %T') no /root/k3freq_commit: stop"; exit 1; }
[ -f $B/runs/K3FREQ_${V}_DONE ] && { echo "$(date '+%F %T') k3freq $V already done"; exit 0; }
until [ -f $B/runs/K3FREQ_${V}_GO ]; do sleep 300; done
sleep 60                                    # at a restart: give the system time to settle
while pgrep -f "[c]hain_newell.sh" > /dev/null; do sleep 300; done
free=$(df --output=avail -BG /root | tail -1 | tr -dc 0-9)
[ "$free" -ge 6 ] || { echo "$(date '+%F %T') only ${free} GB free: stop"; exit 1; }
C=$(cat /root/k3freq_commit)
cd /root/magnum.np || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git fetch -q --depth 50 origin feature/newell-pbc && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
cd $B || exit 1
L=queue_k3freq_$V.log
echo "$(date '+%F %T') k3freq $V chain start at $(git log --oneline -1)" >> $L
python run_queue.py jobs/k3freq_$V.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 12 \
    --status status_k3freq_$V.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rc=$?; echo "$(date '+%F %T') k3freq $V queue exit code $rc" >> $L
python eval_k3freq.py --prefix $P > runs/k3freq_${V}_report.md 2>> $L
e=$?; echo "$(date '+%F %T') k3freq $V report exit code $e" >> $L
[ $rc -eq 0 ] && [ $e -eq 0 ] && touch runs/K3FREQ_${V}_DONE
