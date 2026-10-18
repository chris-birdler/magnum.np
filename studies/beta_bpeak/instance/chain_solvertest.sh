#!/bin/bash
# Solver test at 150 mT (decision Chris 2026-10-18): jobs/solvertest.txt (FD against Newell from the state of K1_s_nx
# after the 100 mT stage), then runs/solvertest_report.md. Starts when the corner queue of chain_corners_nx.sh has
# ended (the precision queue then runs beside it, so that no GPU is idle). Own worktree /root/magnum.np-st at
# /root/st_commit; runs/ is a link to the shared run folder.
#   start: echo <commit> > /root/st_commit; cp chain_solvertest.sh /root/; add "nohup bash /root/chain_solvertest.sh >>
#   /root/st_chain.out 2>&1 &" to /root/onstart.sh; start it the same way.
#   bash chain_solvertest.sh now 1,3   starts at once on the GPUs 1 and 3 (idle GPUs, Chris 2026-10-18)
B=/root/magnum.np/studies/beta_bpeak
W=/root/magnum.np-st
T=$W/studies/beta_bpeak
exec 3> /root/st_chain.lock
flock -n 3 || { echo "$(date '+%F %T') st chain already running"; exit 0; }
[ -f /root/st_commit ] || { echo "$(date '+%F %T') no /root/st_commit: stop"; exit 1; }
[ -f $B/runs/SOLVERTEST_DONE ] && { echo "$(date '+%F %T') solvertest already done"; exit 0; }
C=$(cat /root/st_commit)
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
L=$B/queue_solvertest.log
echo "$(date '+%F %T') st chain (re)start at $(git log --oneline -1 | cut -c1-60)" >> $L
# only the lines after the last (re)start of the corner chain count (the stopped start of 2026-10-17 is also in the log)
GPUS=${2:-0,1,2,3}
[ "$1" = "now" ] || until tac $B/queue_corners_nx.log 2> /dev/null | sed '/cx chain (re)start/q' | grep -q "corners_nx queue exit code"; do
    sleep 120
done
echo "$(date '+%F %T') start on the GPUs $GPUS" >> $L
python run_queue.py jobs/solvertest.txt --gpus $GPUS --min_free_gb 3 --retries 1 --job_timeout_h 24 \
    --status status_solvertest.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rq=$?; echo "$(date '+%F %T') solvertest queue exit code $rq" >> $L
python eval_solvertest.py > runs/solvertest_report.md 2>> $L
e=$?; echo "$(date '+%F %T') solvertest report exit code $e" >> $L
[ $rq -eq 0 ] && [ $e -eq 0 ] && touch runs/SOLVERTEST_DONE
