#!/bin/bash
# The box dispatcher (Chris 2026-10-18: one queue per box, one job per GPU). It runs run_queue.py --spool until
# spool/EXIT exists and takes the job files that submit.sh puts into spool/ (priority in the file name). Post-processing
# of the snapshots and accumulators under a supervisor loop for the production names (<design name>_nx). Own worktree
# /root/magnum.np-disp at /root/disp_commit; runs/ and spool/ are links to the shared folders.
#   install: echo <commit> > /root/disp_commit; cp dispatcher.sh submit.sh /root/; add "nohup bash /root/dispatcher.sh >>
#   /root/disp.out 2>&1 &" to /root/onstart.sh; start it the same way.
#   new code: echo <commit> > /root/disp_commit; touch spool/EXIT (the running jobs end, the queued stay in spool/);
#   when the dispatcher has ended: rm spool/EXIT and start it again.
B=/root/magnum.np/studies/beta_bpeak
W=/root/magnum.np-disp
T=$W/studies/beta_bpeak
exec 3> /root/disp.lock
flock -n 3 || { echo "$(date '+%F %T') dispatcher already running"; exit 0; }
[ -f /root/disp_commit ] || { echo "$(date '+%F %T') no /root/disp_commit: stop"; exit 1; }
C=$(cat /root/disp_commit)
if [ ! -d $W ]; then
    cd /root/magnum.np && git fetch -q --depth 50 origin main && git worktree add -q -f --detach $W "$C" || {
        echo "$(date '+%F %T') worktree at $C failed: stop"; exit 1; }
fi
cd $W || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git -C /root/magnum.np fetch -q --depth 50 origin main && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
mkdir -p $B/spool
for x in runs spool; do
    if [ ! -L $T/$x ]; then
        [ -e $T/$x ] && mv $T/$x $T/${x}_worktree
        ln -s $B/$x $T/$x
    fi
done
cd $T || exit 1
L=$B/dispatcher.log
echo "$(date '+%F %T') dispatcher (re)start at $(git log --oneline -1 | cut -c1-60)" >> $L
( while [ ! -f spool/EXIT ]; do
      python instance_postproc.py --runs runs --stop spool/EXIT --workers 4 --pattern '^(S\d{3}|C\d|K\d_\d)_(nx|hs|ad|as)$' \
          >> $B/postproc_dispatcher.log 2>&1 && break
      echo "$(date '+%F %T') instance_postproc (dispatcher) restarted" >> runs/ALERT
      sleep 60
  done ) &
python run_queue.py --spool spool --gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 24 --exclusive_gpu \
    --status status_dispatcher.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
echo "$(date '+%F %T') dispatcher ended, exit code $?" >> $L
wait
