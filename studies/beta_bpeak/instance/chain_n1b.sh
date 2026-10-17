#!/bin/bash
# PLAN 4.12 V3 (PROTOKOLL 7.109; Chris delegated the decision after V1, 2026-10-17): N0b stage 1 and 2 again with the
# final stray-field setting (prefix N1b, jobs/final_n1b.txt: F first, then I3, I15, T15, T3; the T arms wait for
# N1b_F/init.pt) on the 4 GPUs, then the report runs/v3_n1b_report.md (eval_v3.py). Own worktree /root/magnum.np-n1b at
# /root/n1b_commit, runs/ is a link to the shared run folder. Logs every exit code; N1B_DONE only if all are 0.
#   start: echo <commit> > /root/n1b_commit; cp chain_n1b.sh /root/; add "nohup bash /root/chain_n1b.sh >>
#   /root/n1b_chain.out 2>&1 &" to onstart; start it the same way.
B=/root/magnum.np/studies/beta_bpeak
W=/root/magnum.np-n1b
T=$W/studies/beta_bpeak
exec 3> /root/n1b_chain.lock
flock -n 3 || { echo "$(date '+%F %T') n1b chain already running"; exit 0; }
[ -f /root/n1b_commit ] || { echo "$(date '+%F %T') no /root/n1b_commit: stop"; exit 1; }
[ -f $B/runs/N1B_DONE ] && { echo "$(date '+%F %T') n1b already done"; exit 0; }
free=$(df --output=avail -BG /root | tail -1 | tr -dc 0-9)
[ "$free" -ge 6 ] || { echo "$(date '+%F %T') only ${free} GB free: stop"; exit 1; }
C=$(cat /root/n1b_commit)
if [ ! -d $W ]; then
    cd /root/magnum.np && git fetch -q --depth 50 origin feature/newell-pbc && git worktree add -q -f --detach $W "$C" || {
        echo "$(date '+%F %T') worktree at $C failed: stop"; exit 1; }
fi
cd $W || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git -C /root/magnum.np fetch -q --depth 50 origin feature/newell-pbc && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
if [ ! -L $T/runs ]; then
    [ -e $T/runs ] && mv $T/runs $T/runs_worktree
    ln -s $B/runs $T/runs
fi
cd $T || exit 1
L=$B/queue_n1b.log
echo "$(date '+%F %T') n1b chain start at $(git log --oneline -1)" >> $L
python run_queue.py jobs/final_n1b.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 30 --gate_maxiter 0 \
    --status status_n1b.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rq=$?; echo "$(date '+%F %T') n1b queue exit code $rq" >> $L
python eval_v3.py > runs/v3_n1b_report.md 2>> $L
e=$?; echo "$(date '+%F %T') n1b report exit code $e" >> $L
[ $rq -eq 0 ] && [ $e -eq 0 ] && touch runs/N1B_DONE
