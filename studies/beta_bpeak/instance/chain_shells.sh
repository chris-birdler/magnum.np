#!/bin/bash
# Shell convergence of the reference solver (decision Chris 2026-10-16, option A of PROTOKOLL 7.108): mode "shells" of
# check_mesh_energy.py on the CPU (CUDA_DEVICE=-1; the 4 GPUs run check 1 and the integrator test), one process per
# mesh (dx 3, 1.5, 0.75; 8 threads each), then the report runs/shells_<prefix>.md. Own worktree /root/magnum.np-shells at
# /root/shells_commit, runs/ is a link to the shared run folder. Logs every exit code; SHELLS_DONE only if all are 0.
#   start: echo <commit> > /root/shells_commit; cp chain_shells.sh /root/; setsid nohup bash /root/chain_shells.sh N0b
#   >> /root/shells_chain.out 2>&1 &
P=${1:-N0b}
B=/root/magnum.np/studies/beta_bpeak
W=/root/magnum.np-shells
T=$W/studies/beta_bpeak
exec 3> /root/shells_chain.lock
flock -n 3 || { echo "$(date '+%F %T') shells chain already running"; exit 0; }
[ -f /root/shells_commit ] || { echo "$(date '+%F %T') no /root/shells_commit: stop"; exit 1; }
[ -f $B/runs/SHELLS_DONE ] && { echo "$(date '+%F %T') shells already done"; exit 0; }
C=$(cat /root/shells_commit)
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
L=$B/queue_shells.log
echo "$(date '+%F %T') shells chain start at $(git log --oneline -1)" >> $L
pids=""
for dx in 3 1.5 0.75; do
    CUDA_DEVICE=-1 OMP_NUM_THREADS=8 nice -n 5 python check_mesh_energy.py shells $P $dx > runs/shells_${P}_dx$dx.log 2>&1 &
    pids="$pids $!"
done
rc=0
for p in $pids; do wait $p; r=$?; echo "$(date '+%F %T') shells process $p exit code $r" >> $L; [ $r -eq 0 ] || rc=$r; done
python check_mesh_energy.py shells_report $P > runs/shells_$P.md 2>> $L
e=$?; echo "$(date '+%F %T') shells report exit code $e" >> $L
[ $rc -eq 0 ] && [ $e -eq 0 ] && touch runs/SHELLS_DONE
