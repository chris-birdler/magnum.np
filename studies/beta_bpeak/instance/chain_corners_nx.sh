#!/bin/bash
# Corner test and precision test of PLAN 4.13 (approved by Chris 2026-10-17). Production setup of a new main study:
# the 12 corner runs (jobs/corners_nx.txt, d/l_ex 212 and 300 alternate), then the precision test (jobs/precision_nx.txt,
# from the init.pt of K1_0_nx … K4_0_nx), then runs/corners_nx_report.md and runs/precision_nx_report.md.
# As chain_sobol.sh: post-processing of the snapshots and accumulators under a supervisor loop (instance_postproc.py,
# pattern of the corner runs), gate retry with 200 000 iterations (run_queue default), one retry, job timeout 24 h
# (PLAN 4.13 S11), --min_free_gb 3 with the disk reservation of run_queue.expected_bytes.
# Own worktree /root/magnum.np-cx at /root/cx_commit; runs/ is a link to the shared run folder. Every exit code is
# logged; CORNERS_NX_DONE only if all are 0.
#   install: echo <commit> > /root/cx_commit; cp chain_corners_nx.sh /root/; add "nohup bash /root/chain_corners_nx.sh
#   >> /root/cx_chain.out 2>&1 &" to /root/onstart.sh; start it the same way.
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
echo "$(date '+%F %T') cx chain (re)start at $(git log --oneline -1)" >> $L
# post-processing under a supervisor loop: restarted if it dies, until CORNERS_NX_QUEUE_DONE and no run left
( while true; do
      python instance_postproc.py --runs runs --stop runs/CORNERS_NX_QUEUE_DONE --workers 4 \
          --pattern '^K\d_\d_nx$' >> $B/postproc_corners_nx.log 2>&1 && break
      echo "$(date '+%F %T') instance_postproc exited with an error: restart in 60 s" >> $B/postproc_corners_nx.log
      echo "$(date '+%F %T') instance_postproc (corners_nx) restarted" >> runs/ALERT
      sleep 60
  done ) &
PP=$!
Q="--gpus 0,1,2,3 --min_free_gb 3 --retries 1 --job_timeout_h 24 --stop_file runs/STOP_NEW_JOBS"
python run_queue.py jobs/corners_nx.txt $Q --status status_corners_nx.json >> $L 2>&1
rq=$?; echo "$(date '+%F %T') corners_nx queue exit code $rq" >> $L
python run_queue.py jobs/precision_nx.txt $Q --status status_precision_nx.json >> $L 2>&1
rp=$?; echo "$(date '+%F %T') precision_nx queue exit code $rp" >> $L
touch runs/CORNERS_NX_QUEUE_DONE
wait $PP
echo "$(date '+%F %T') corners_nx post-processing finished" >> $L
python eval_corners_nx.py > runs/corners_nx_report.md 2>> $L
e1=$?; echo "$(date '+%F %T') corners_nx report exit code $e1" >> $L
python eval_precision.py > runs/precision_nx_report.md 2>> $L
e2=$?; echo "$(date '+%F %T') precision_nx report exit code $e2" >> $L
[ $rq -eq 0 ] && [ $rp -eq 0 ] && [ $e1 -eq 0 ] && [ $e2 -eq 0 ] && touch runs/CORNERS_NX_DONE
