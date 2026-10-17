#!/bin/bash
# Corner test (decision Chris 2026-10-17: first the corners, then decide on a new main study): the 12 corner runs of the
# main study with the final stray-field setting (jobs/corners_nx.txt; only the solver changes) on the 4 GPUs, then
# runs/corners_nx_report.md (eval_corners_nx.py). Own worktree /root/magnum.np-cx at /root/cx_commit, runs/ is a link
# to the shared run folder. Logs every exit code; CORNERS_NX_DONE only if all are 0.
#   start: echo <commit> > /root/cx_commit; cp chain_corners_nx.sh /root/; add "nohup bash /root/chain_corners_nx.sh >>
#   /root/cx_chain.out 2>&1 &" to onstart; start it the same way.
B=/root/magnum.np/studies/beta_bpeak
W=/root/magnum.np-cx
T=$W/studies/beta_bpeak
exec 3> /root/cx_chain.lock
flock -n 3 || { echo "$(date '+%F %T') cx chain already running"; exit 0; }
[ -f /root/cx_commit ] || { echo "$(date '+%F %T') no /root/cx_commit: stop"; exit 1; }
[ -f $B/runs/CORNERS_NX_DONE ] && { echo "$(date '+%F %T') corners_nx already done"; exit 0; }
C=$(cat /root/cx_commit)
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
L=$B/queue_corners_nx.log
echo "$(date '+%F %T') cx chain start at $(git log --oneline -1)" >> $L
python run_queue.py jobs/corners_nx.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 12 --gate_maxiter 0 \
    --status status_corners_nx.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rq=$?; echo "$(date '+%F %T') corners_nx queue exit code $rq" >> $L
python eval_corners_nx.py > runs/corners_nx_report.md 2>> $L
e=$?; echo "$(date '+%F %T') corners_nx report exit code $e" >> $L
[ $rq -eq 0 ] && [ $e -eq 0 ] && touch runs/CORNERS_NX_DONE
