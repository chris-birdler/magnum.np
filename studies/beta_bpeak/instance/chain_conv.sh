#!/bin/bash
# B0, step 1 of the mesh convergence series (PLAN 4.10): jobs/conv_b0.txt (20 jobs: d/l_ex 96, 9 mT, seeds 941 … 944;
# per seed dx 0.75, dx 1.5 and dx 3 from the transferred dx 0.75 state, dx 1.5 and dx 3 independent; no snapshots, no
# accumulators). Starts only after the go of Chris (runs/CONV_GO, written by hand) and after the frequency chain C′ has
# ended: it polls for runs/FREQ_QUEUE_DONE (written after the C′ queue, also after a failed job) and then takes the
# freq lock (held by chain_freq.sh until its post-processing ends). Then it removes the freq chain from onstart, pins
# its own commit (/root/conv_commit) and runs the queue. A stop file (runs/STOP_NEW_JOBS) stops the queue before any
# job. Job timeout 48 h: a dx 0.75 run needs ≤ 80 000 relaxation iterations × 1.07 s + 12 cycles × 1.28 h ≈ 39 h in
# the worst case (measured times, PROTOKOLL 7.80). At the end it writes the pre-registered report (eval_conv.py).
#   install (after the go): echo <commit> > /root/conv_commit; cp chain_conv.sh /root/; add the onstart line;
#   start with nohup; touch runs/CONV_GO.
B=/root/magnum.np/studies/beta_bpeak
exec 3> /root/conv_chain.lock
flock -n 3 || { echo "$(date '+%F %T') conv chain already running"; exit 0; }
until [ -f $B/runs/CONV_GO ]; do sleep 300; done
until [ -f $B/runs/FREQ_QUEUE_DONE ]; do sleep 300; done
exec 4> /root/freq_chain.lock
flock 4                                  # waits for the end of the freq post-processing; held until the end
sed -i "/chain_freq.sh/d" /root/onstart.sh
[ -f /root/conv_commit ] || { echo "$(date '+%F %T') no /root/conv_commit: stop"; exit 1; }
C=$(cat /root/conv_commit)
cd /root/magnum.np || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git fetch -q --depth 50 origin study/beta-bpeak-pinning && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
cd $B || exit 1
L=queue_conv.log
echo "$(date '+%F %T') conv chain (re)start at $(git log --oneline -1)" >> $L
python run_queue.py jobs/conv_b0.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 48 \
    --status status_conv.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rc=$?; echo "$(date '+%F %T') conv queue exit code $rc" >> $L
touch runs/CONV_QUEUE_DONE
python eval_conv.py > runs/conv_b0_report.md 2>> $L
echo "$(date '+%F %T') conv report exit code $?" >> $L
[ $rc -eq 0 ] && touch runs/CONV_DONE
