#!/bin/bash
# PLAN 4.12 V1 (decision Chris 2026-10-17): validation of the stray-field candidates against the alias sum. Own worktree
# /root/magnum.np-v1 at /root/v1_commit, runs/ is a link to the shared run folder. The alias reference runs on the CPU,
# the Newell candidates on GPUs 0 … 2. Then the report runs/v1_N0b.md. Logs every exit code; V1_DONE only if all are 0.
#   start: echo <commit> > /root/v1_commit; cp chain_v1.sh /root/; setsid nohup bash /root/chain_v1.sh >> /root/v1_chain.out 2>&1 &
B=/root/magnum.np/studies/beta_bpeak
W=/root/magnum.np-v1
T=$W/studies/beta_bpeak
exec 3> /root/v1_chain.lock
flock -n 3 || { echo "$(date '+%F %T') v1 chain already running"; exit 0; }
[ -f /root/v1_commit ] || { echo "$(date '+%F %T') no /root/v1_commit: stop"; exit 1; }
[ -f $B/runs/V1_DONE ] && { echo "$(date '+%F %T') v1 already done"; exit 0; }
C=$(cat /root/v1_commit)
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
L=$B/queue_v1.log
echo "$(date '+%F %T') v1 chain start at $(git log --oneline -1)" >> $L
CUDA_DEVICE=-1 OMP_NUM_THREADS=8 python check_demag_alias.py v1_ref N0b > runs/v1_ref_N0b.log 2>&1 &
p0=$!
CUDA_DEVICE=0 python check_demag_alias.py v1_newell N0b 1 2 4 > runs/v1_newell_a.log 2>&1 &
p1=$!
CUDA_DEVICE=1 python check_demag_alias.py v1_newell N0b 8 > runs/v1_newell_b.log 2>&1 &
p2=$!
CUDA_DEVICE=2 python check_demag_alias.py v1_newell N0b 2x 4x > runs/v1_newell_c.log 2>&1 &
p3=$!
rc=0
for p in $p0 $p1 $p2 $p3; do wait $p; r=$?; echo "$(date '+%F %T') v1 process $p exit code $r" >> $L; [ $r -eq 0 ] || rc=$r; done
python check_demag_alias.py v1_report N0b > runs/v1_N0b.md 2>> $L
e=$?; echo "$(date '+%F %T') v1 report exit code $e" >> $L
[ $rc -eq 0 ] && [ $e -eq 0 ] && touch runs/V1_DONE
