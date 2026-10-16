#!/bin/bash
# Integrator test on dx 0.75 (decision Chris 2026-10-16, option C2, PROTOKOLL 7.108): the F runs of seeds 941 and 942
# from their relaxed state with a 10 x tighter RKF45 tolerance of the drive cycles (atol 1e-7, rtol 1e-6;
# jobs/newell_tol075.txt). s941 on GPU 3 at once, s942 on GPU 2 after check 2 (check_mesh_energy.py solver N0b) has
# ended. GPUs 0 and 1 run check 1 (chain_eps075.sh) from /root/magnum.np at its own commit; this chain uses its own
# worktree /root/magnum.np-tol at /root/tol075_commit, with runs/ as a link to the shared run folder.
# Abort rule (fixed before the runs): a run whose stdout.log has no new line for 3 x 4 700 s (one cycle 3 x longer than
# in stage 2) is stopped by the kill of its PID; run_queue does not repeat it (--retries 0).
# Report runs/tol075_report.md (eval_eps075.py tol). Logs every exit code; TOL075_DONE only if all end with code 0.
#   install: echo <commit> > /root/tol075_commit; cp chain_tol075.sh /root/; add "nohup bash /root/chain_tol075.sh >>
#   /root/tol075_chain.out 2>&1 &" to onstart; start it the same way.
B=/root/magnum.np/studies/beta_bpeak
W=/root/magnum.np-tol
T=$W/studies/beta_bpeak
RUNS="N0b_F_s941_atol1e-07_rtol1e-06 N0b_F_s942_atol1e-07_rtol1e-06"
LIMIT=14100
exec 3> /root/tol075_chain.lock
flock -n 3 || { echo "$(date '+%F %T') tol075 chain already running"; exit 0; }
[ -f /root/tol075_commit ] || { echo "$(date '+%F %T') no /root/tol075_commit: stop"; exit 1; }
[ -f $B/runs/TOL075_DONE ] && { echo "$(date '+%F %T') tol075 already done"; exit 0; }
free=$(df --output=avail -BG /root | tail -1 | tr -dc 0-9)
[ "$free" -ge 6 ] || { echo "$(date '+%F %T') only ${free} GB free: stop"; exit 1; }
C=$(cat /root/tol075_commit)
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
L=$B/queue_tol075.log
echo "$(date '+%F %T') tol075 chain start at $(git log --oneline -1)" >> $L
Q="--min_free_gb 3 --retries 0 --job_timeout_h 24 --gate_maxiter 0 --stop_file runs/STOP_NEW_JOBS"
python run_queue.py jobs/newell_tol075.txt --gpus 3 --shard 0/2 --status status_tol075_a.json $Q >> $L 2>&1 &
qa=$!
qb=""
started_b=0
while true; do
    if [ $started_b -eq 0 ] && ! pgrep -f "[c]heck_mesh_energy.py solver" > /dev/null; then
        python run_queue.py jobs/newell_tol075.txt --gpus 2 --shard 1/2 --status status_tol075_b.json $Q >> $L 2>&1 &
        qb=$!
        started_b=1
    fi
    for r in $RUNS; do
        f=$B/runs/$r/stdout.log
        { [ -f $f ] && [ ! -f $B/runs/$r/DONE ] && [ ! -f $B/runs/$r/ABORTED_SLOW ]; } || continue
        age=$(( $(date +%s) - $(stat -c %Y $f) ))
        if [ $age -gt $LIMIT ]; then
            pid=$(pgrep -f "[r]un_loops\.py .*runs/$r\$")
            echo "$(date '+%F %T') ABORT $r: no new line in stdout.log for $age s (limit $LIMIT s), kill PID $pid" >> $L
            touch $B/runs/$r/ABORTED_SLOW
            [ -n "$pid" ] && kill $pid
        fi
    done
    a_on=0; b_on=0
    kill -0 $qa 2> /dev/null && a_on=1
    [ -n "$qb" ] && kill -0 $qb 2> /dev/null && b_on=1
    [ $a_on -eq 0 ] && [ $started_b -eq 1 ] && [ $b_on -eq 0 ] && break
    sleep 300
done
wait $qa; ra=$?
wait $qb; rb=$?
echo "$(date '+%F %T') tol075 queue exit codes $ra $rb" >> $L
python eval_eps075.py tol > runs/tol075_report.md 2>> $L
e=$?; echo "$(date '+%F %T') tol075 report exit code $e" >> $L
[ $ra -eq 0 ] && [ $rb -eq 0 ] && [ $e -eq 0 ] && touch runs/TOL075_DONE
