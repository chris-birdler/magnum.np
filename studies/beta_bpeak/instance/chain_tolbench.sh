#!/bin/bash
# Cost of the cycle tolerance and of the final stray-field setting on the mesh of the main study (decision Chris
# 2026-10-17: measure, then decide): jobs/tolbench.txt (corners K1_0, K2_0; variants a5 / a6 / tight) on the 4 GPUs
# after the V3g queue (chain_n1g.sh) has ended, then runs/tolbench_report.md (eval_tolbench.py). Own worktree
# /root/magnum.np-tb at /root/tb_commit, runs/ is a link to the shared run folder. Logs every exit code; TOLBENCH_DONE
# only if all are 0.
#   start: echo <commit> > /root/tb_commit; cp chain_tolbench.sh /root/; add "nohup bash /root/chain_tolbench.sh >>
#   /root/tb_chain.out 2>&1 &" to onstart; start it the same way.
B=/root/magnum.np/studies/beta_bpeak
W=/root/magnum.np-tb
T=$W/studies/beta_bpeak
exec 3> /root/tb_chain.lock
flock -n 3 || { echo "$(date '+%F %T') tb chain already running"; exit 0; }
[ -f /root/tb_commit ] || { echo "$(date '+%F %T') no /root/tb_commit: stop"; exit 1; }
[ -f $B/runs/TOLBENCH_DONE ] && { echo "$(date '+%F %T') tolbench already done"; exit 0; }
C=$(cat /root/tb_commit)
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
L=$B/queue_tolbench.log
echo "$(date '+%F %T') tb chain start at $(git log --oneline -1); waits for the end of the V3g queue" >> $L
until grep -q "[n]1g queue exit code" $B/queue_n1g.log 2> /dev/null; do sleep 300; done
python run_queue.py jobs/tolbench.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 6 --gate_maxiter 0 \
    --status status_tolbench.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rq=$?; echo "$(date '+%F %T') tolbench queue exit code $rq" >> $L
python eval_tolbench.py > runs/tolbench_report.md 2>> $L
e=$?; echo "$(date '+%F %T') tolbench report exit code $e" >> $L
[ $rq -eq 0 ] && [ $e -eq 0 ] && touch runs/TOLBENCH_DONE
