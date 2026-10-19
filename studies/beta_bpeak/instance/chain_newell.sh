#!/bin/bash
# B0b with the reference stray-field solver (PLAN 4.11, PROTOKOLL 7.105): stage $1 (1: I3 and I15; 2: F, T15, T3) of
# jobs/newell_b0_s$1.txt (conv_jobs.py --newell $1; prefix N0b, --strayfield newell --strayfield_shells 4). Starts only
# after the go of Chris for this stage (runs/NEWELL_S$1_GO, written by hand). Pins its own commit on feature/newell-pbc
# (/root/newell_commit), stops below 6 GB free disk space, runs the queue, writes the pre-registered report (eval_conv.py
# --prefix N0b) and logs every exit code; NEWELL_S$1_DONE only if the queue and the report end with code 0. A stop file
# (runs/STOP_NEW_JOBS) stops the queue before any job. Job timeout 48 h and --gate_maxiter 0 as the B0b chain.
#   install: echo <commit> > /root/newell_commit; cp chain_newell.sh /root/; add "nohup bash /root/chain_newell.sh 1
#   >> /root/newell_chain.out 2>&1 &" to onstart; start it the same way; touch runs/NEWELL_S1_GO.
S=$1
[ "$S" = 1 ] || [ "$S" = 2 ] || { echo "$(date '+%F %T') stage must be 1 or 2"; exit 1; }
B=/root/magnum.np/studies/beta_bpeak
exec 3> /root/newell_chain.lock
flock -n 3 || { echo "$(date '+%F %T') newell chain already running"; exit 0; }
[ -f /root/newell_commit ] || { echo "$(date '+%F %T') no /root/newell_commit: stop"; exit 1; }
[ -f $B/runs/NEWELL_S${S}_DONE ] && { echo "$(date '+%F %T') newell stage $S already done"; exit 0; }
until [ -f $B/runs/NEWELL_S${S}_GO ]; do sleep 300; done
sleep 60                                    # at a restart: give the system time to settle
free=$(df --output=avail -BG /root | tail -1 | tr -dc 0-9)
[ "$free" -ge 6 ] || { echo "$(date '+%F %T') only ${free} GB free: stop"; exit 1; }
C=$(cat /root/newell_commit)
cd /root/magnum.np || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git fetch -q --depth 50 origin feature/newell-pbc && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
cd $B || exit 1
L=queue_newell_s$S.log
echo "$(date '+%F %T') newell stage $S chain start at $(git log --oneline -1)" >> $L
python run_queue.py jobs/newell_b0_s$S.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 48 \
    --gate_maxiter 0 --status status_newell_s$S.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rc=$?; echo "$(date '+%F %T') newell stage $S queue exit code $rc" >> $L
python eval_conv.py --prefix N0b > runs/newell_b0_s${S}_report.md 2>> $L
e=$?; echo "$(date '+%F %T') newell stage $S report exit code $e" >> $L
[ $rc -eq 0 ] && [ $e -eq 0 ] && touch runs/NEWELL_S${S}_DONE
