#!/bin/bash
# Option 1' of PROTOKOLL 7.108 (decision Chris 2026-10-16): (1) static check of the extrapolated kernel on the CPU
# (check_mesh_energy.py shells_ext N0b <dx> for dx 3, 1.5, 0.75, then the report runs/shells_N0b.md with the column
# E(4x)); (2) the 8 dx 1.5 runs with the extrapolated kernel in the drive cycles (jobs/newell_kext15.txt) on the 4 GPUs
# after check 1 and the integrator test have ended; report runs/kext15_report.md (eval_kext.py). Own worktree
# /root/magnum.np-kext at /root/kext_commit, runs/ is a link to the shared run folder. Logs every exit code; KEXT_DONE
# only if all are 0.
#   start: echo <commit> > /root/kext_commit; cp chain_kext.sh /root/; add "nohup bash /root/chain_kext.sh >>
#   /root/kext_chain.out 2>&1 &" to onstart; start it the same way.
B=/root/magnum.np/studies/beta_bpeak
W=/root/magnum.np-kext
T=$W/studies/beta_bpeak
exec 3> /root/kext_chain.lock
flock -n 3 || { echo "$(date '+%F %T') kext chain already running"; exit 0; }
[ -f /root/kext_commit ] || { echo "$(date '+%F %T') no /root/kext_commit: stop"; exit 1; }
[ -f $B/runs/KEXT_DONE ] && { echo "$(date '+%F %T') kext already done"; exit 0; }
C=$(cat /root/kext_commit)
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
L=$B/queue_kext.log
echo "$(date '+%F %T') kext chain start at $(git log --oneline -1)" >> $L
rs=0
if ! grep -q "Reading: extrapolated kernel valid\|Reading: not valid" runs/shells_N0b.md 2> /dev/null; then
    pids=""
    for dx in 3 1.5 0.75; do
        CUDA_DEVICE=-1 OMP_NUM_THREADS=8 nice -n 5 python check_mesh_energy.py shells_ext N0b $dx > runs/shells_ext_N0b_dx$dx.log 2>&1 &
        pids="$pids $!"
    done
    for p in $pids; do wait $p; r=$?; echo "$(date '+%F %T') shells_ext process $p exit code $r" >> $L; [ $r -eq 0 ] || rs=$r; done
    python check_mesh_energy.py shells_report N0b > runs/shells_N0b.md 2>> $L
    r=$?; echo "$(date '+%F %T') shells report exit code $r" >> $L; [ $r -eq 0 ] || rs=$r
fi
until grep -q "[e]ps075 queue exit code" $B/queue_eps075.log 2> /dev/null && \
      grep -q "[t]ol075 queue exit codes" $B/queue_tol075.log 2> /dev/null; do sleep 120; done
python run_queue.py jobs/newell_kext15.txt --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 6 --gate_maxiter 0 \
    --status status_kext.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rq=$?; echo "$(date '+%F %T') kext queue exit code $rq" >> $L
python eval_kext.py > runs/kext15_report.md 2>> $L
e=$?; echo "$(date '+%F %T') kext report exit code $e" >> $L
[ $rs -eq 0 ] && [ $rq -eq 0 ] && [ $e -eq 0 ] && touch runs/KEXT_DONE
