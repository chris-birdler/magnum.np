#!/bin/bash
# PLAN 4.11, after stage 2 (decision Chris 2026-10-15, PROTOKOLL 7.107): the two stage 1 runs with a drive amplitude
# 5 … 6 % below the target (N0b_I15_s941, N0b_I3_s944) again with --b_tol 0.01 (jobs/newell_b0_s1b.txt), the same drive
# rule as stage 2. Waits for runs/NEWELL_S2_DONE (the stage 2 chain must have ended), pins its own commit on
# feature/newell-pbc (/root/newell_s1b_commit), renames the old run folders to <name>_nobtol (not deleted; local copies
# with md5 manifest exist), runs the queue and writes the report over all arms (eval_conv.py --prefix N0b) to
# runs/newell_b0_s1b_report.md. Logs every exit code; NEWELL_S1B_DONE only if the queue and the report end with code 0.
# Its own lock and script file: the stage 2 chain script must not be changed while it runs.
#   install: echo <commit> > /root/newell_s1b_commit; cp chain_newell_s1b.sh /root/; add "nohup bash
#   /root/chain_newell_s1b.sh >> /root/newell_s1b_chain.out 2>&1 &" to onstart; start it the same way.
B=/root/magnum.np/studies/beta_bpeak
exec 3> /root/newell_s1b_chain.lock
flock -n 3 || { echo "$(date '+%F %T') newell s1b chain already running"; exit 0; }
[ -f /root/newell_s1b_commit ] || { echo "$(date '+%F %T') no /root/newell_s1b_commit: stop"; exit 1; }
[ -f $B/runs/NEWELL_S1B_DONE ] && { echo "$(date '+%F %T') newell s1b already done"; exit 0; }
until [ -f $B/runs/NEWELL_S2_DONE ]; do sleep 300; done
sleep 60
while pgrep -f "[c]hain_newell.sh" > /dev/null; do sleep 120; done
free=$(df --output=avail -BG /root | tail -1 | tr -dc 0-9)
[ "$free" -ge 6 ] || { echo "$(date '+%F %T') only ${free} GB free: stop"; exit 1; }
C=$(cat /root/newell_s1b_commit)
cd /root/magnum.np || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git fetch -q --depth 50 origin feature/newell-pbc && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
cd $B || exit 1
L=queue_newell_s1b.log
echo "$(date '+%F %T') newell s1b chain start at $(git log --oneline -1)" >> $L
for n in N0b_I15_s941 N0b_I3_s944; do
    if [ -d runs/$n ] && ! grep -q "b_tol" runs/$n/config.json 2>/dev/null; then
        [ -e runs/${n}_nobtol ] || mv runs/$n runs/${n}_nobtol
        echo "$(date '+%F %T') renamed runs/$n to runs/${n}_nobtol" >> $L
    elif [ -d runs/$n ] && ! python -c "import json,sys; sys.exit(0 if json.load(open('runs/$n/config.json')).get('b_tol', 0) > 0 else 1)"; then
        [ -e runs/${n}_nobtol ] || mv runs/$n runs/${n}_nobtol
        echo "$(date '+%F %T') renamed runs/$n (b_tol 0) to runs/${n}_nobtol" >> $L
    fi
done
python run_queue.py jobs/newell_b0_s1b.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 12 \
    --gate_maxiter 0 --status status_newell_s1b.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rc=$?; echo "$(date '+%F %T') newell s1b queue exit code $rc" >> $L
python eval_conv.py --prefix N0b > runs/newell_b0_s1b_report.md 2>> $L
e=$?; echo "$(date '+%F %T') newell s1b report exit code $e" >> $L
[ $rc -eq 0 ] && [ $e -eq 0 ] && touch runs/NEWELL_S1B_DONE
