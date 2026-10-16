#!/bin/bash
# Abort rule of chain_tol075.sh (PROTOKOLL 7.108): a run whose stdout.log has no new line for 3 x 4 700 s is stopped.
# magnum.np sets the process title without the arguments, so this script finds the PID of a run by the link of its
# stdout (/proc/<pid>/fd/1 = runs/<run>/stdout.log). Ends when chain_tol075.sh has logged the exit codes of its queues.
#   install: cp watchdog_tol075.sh /root/; add "nohup bash /root/watchdog_tol075.sh >> /root/tol075_watchdog.out 2>&1 &"
#   to onstart; start it the same way.
B=/root/magnum.np/studies/beta_bpeak
L=$B/queue_tol075.log
LIMIT=14100
RUNS="N0b_F_s941_atol1e-07_rtol1e-06 N0b_F_s942_atol1e-07_rtol1e-06"
exec 3> /root/tol075_watchdog.lock
flock -n 3 || { echo "$(date '+%F %T') tol075 watchdog already running"; exit 0; }
echo "$(date '+%F %T') tol075 watchdog start (limit $LIMIT s)" >> $L
until grep -q "[t]ol075 queue exit codes" $L 2> /dev/null; do
    for r in $RUNS; do
        f=$B/runs/$r/stdout.log
        { [ -f $f ] && [ ! -f $B/runs/$r/DONE ]; } || continue
        pid=""
        for p in $(pgrep -f "[m]agnumnp .*run_loops.py"); do
            [ "$(readlink /proc/$p/fd/1)" = "$f" ] && pid=$p
        done
        [ -n "$pid" ] || continue
        age=$(( $(date +%s) - $(stat -c %Y $f) ))
        if [ $age -gt $LIMIT ]; then
            echo "$(date '+%F %T') WATCHDOG ABORT $r: no new line in stdout.log for $age s (limit $LIMIT s), kill PID $pid" >> $L
            kill $pid
        fi
    done
    sleep 300
done
echo "$(date '+%F %T') tol075 watchdog end" >> $L
