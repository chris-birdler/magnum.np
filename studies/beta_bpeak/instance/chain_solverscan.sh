#!/bin/bash
# Solver scan (option (a), decision Chris 2026-10-18): jobs/solverscan.txt (FD at all 7 stages of K1_0_nx … K4_0_nx
# from the same states as the corner runs), then runs/solverscan_report.md. Starts when the corner queue of chain_corners_nx.sh has
# ended (the precision queue then runs beside it, so that no GPU is idle). Own worktree /root/magnum.np-ss at
# /root/ss_commit; runs/ is a link to the shared run folder.
#   start: echo <commit> > /root/ss_commit; cp chain_solverscan.sh /root/; add "nohup bash /root/chain_solverscan.sh >>
#   /root/ss_chain.out 2>&1 &" to /root/onstart.sh; start it the same way.
#   bash chain_solverscan.sh now 1,3   starts at once on the GPUs 1 and 3 (idle GPUs, Chris 2026-10-18)
B=/root/magnum.np/studies/beta_bpeak
W=/root/magnum.np-ss
T=$W/studies/beta_bpeak
exec 3> /root/ss_chain.lock
flock -n 3 || { echo "$(date '+%F %T') ss chain already running"; exit 0; }
[ -f /root/ss_commit ] || { echo "$(date '+%F %T') no /root/ss_commit: stop"; exit 1; }
[ -f $B/runs/SOLVERSCAN_DONE ] && { echo "$(date '+%F %T') solverscan already done"; exit 0; }
C=$(cat /root/ss_commit)
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
L=$B/queue_solverscan.log
echo "$(date '+%F %T') ss chain (re)start at $(git log --oneline -1 | cut -c1-60)" >> $L
# only the lines after the last (re)start of the corner chain count (the stopped start of 2026-10-17 is also in the log)
GPUS=${2:-0,1,2,3}
[ "$1" = "now" ] || until tac $B/queue_corners_nx.log 2> /dev/null | sed '/cx chain (re)start/q' | grep -q "corners_nx queue exit code"; do
    sleep 120
done
echo "$(date '+%F %T') start on the GPUs $GPUS" >> $L
python run_queue.py jobs/solverscan.txt --gpus $GPUS --min_free_gb 3 --retries 1 --job_timeout_h 24 \
    --status status_solverscan.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rq=$?; echo "$(date '+%F %T') solverscan queue exit code $rq" >> $L
python eval_solverscan.py > runs/solverscan_report.md 2>> $L
e=$?; echo "$(date '+%F %T') solverscan report exit code $e" >> $L
[ $rq -eq 0 ] && [ $e -eq 0 ] && touch runs/SOLVERSCAN_DONE
