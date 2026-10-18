#!/bin/bash
# Start "each particle uniform" (Chris 2026-10-18, Bloch points, way 2): jobs/initmode.txt (K2_0_nx, K2_1_nx, K4_1_nx, K3_1_nx
# virgin relaxation from this start, 9 mT stage), then runs/initmode_report.md. Starts when the corner queue of chain_corners_nx.sh has
# ended (the precision queue then runs beside it, so that no GPU is idle). Own worktree /root/magnum.np-ip at
# /root/ip_commit; runs/ is a link to the shared run folder.
#   start: echo <commit> > /root/ip_commit; cp chain_initmode.sh /root/; add "nohup bash /root/chain_initmode.sh >>
#   /root/ip_chain.out 2>&1 &" to /root/onstart.sh; start it the same way.
#   bash chain_initmode.sh now 1,3   starts at once on the GPUs 1 and 3 (idle GPUs, Chris 2026-10-18)
B=/root/magnum.np/studies/beta_bpeak
W=/root/magnum.np-ip
T=$W/studies/beta_bpeak
exec 3> /root/ip_chain.lock
flock -n 3 || { echo "$(date '+%F %T') ip chain already running"; exit 0; }
[ -f /root/ip_commit ] || { echo "$(date '+%F %T') no /root/ip_commit: stop"; exit 1; }
[ -f $B/runs/INITMODE_DONE ] && { echo "$(date '+%F %T') initmode already done"; exit 0; }
C=$(cat /root/ip_commit)
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
L=$B/queue_initmode.log
echo "$(date '+%F %T') ip chain (re)start at $(git log --oneline -1 | cut -c1-60)" >> $L
# only the lines after the last (re)start of the corner chain count (the stopped start of 2026-10-17 is also in the log)
GPUS=${2:-0,1,2,3}
[ "$1" = "now" ] || until tac $B/queue_corners_nx.log 2> /dev/null | sed '/cx chain (re)start/q' | grep -q "corners_nx queue exit code"; do
    sleep 120
done
echo "$(date '+%F %T') start on the GPUs $GPUS" >> $L
python run_queue.py jobs/initmode.txt --gpus $GPUS --min_free_gb 3 --retries 1 --job_timeout_h 24 \
    --status status_initmode.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rq=$?; echo "$(date '+%F %T') initmode queue exit code $rq" >> $L
python eval_initmode.py > runs/initmode_report.md 2>> $L
e=$?; echo "$(date '+%F %T') initmode report exit code $e" >> $L
[ $rq -eq 0 ] && [ $e -eq 0 ] && touch runs/INITMODE_DONE
