#!/bin/bash
# Two checks after PLAN 4.11 (decision Chris 2026-10-16, PROTOKOLL 7.108): (1) the energy balance on dx 0.75: the F runs
# of seeds 941 and 942 again from their relaxed state with 1024 samples per cycle (jobs/newell_eps075.txt, GPUs 0 and
# 1), report runs/eps075_report.md (eval_eps075.py); (2) at the same time on GPU 2: the stray-field energies of the
# N0b states with both solvers (check_mesh_energy.py solver N0b → runs/demag_solver_N0b.md). Pins its own commit on
# feature/newell-pbc (/root/eps075_commit). Logs every exit code; EPS075_DONE only if all end with code 0.
#   install: echo <commit> > /root/eps075_commit; cp chain_eps075.sh /root/; add "nohup bash /root/chain_eps075.sh >>
#   /root/eps075_chain.out 2>&1 &" to onstart; start it the same way.
B=/root/magnum.np/studies/beta_bpeak
exec 3> /root/eps075_chain.lock
flock -n 3 || { echo "$(date '+%F %T') eps075 chain already running"; exit 0; }
[ -f /root/eps075_commit ] || { echo "$(date '+%F %T') no /root/eps075_commit: stop"; exit 1; }
[ -f $B/runs/EPS075_DONE ] && { echo "$(date '+%F %T') eps075 already done"; exit 0; }
sleep 60
while pgrep -f "[c]hain_newell" > /dev/null; do sleep 120; done
free=$(df --output=avail -BG /root | tail -1 | tr -dc 0-9)
[ "$free" -ge 6 ] || { echo "$(date '+%F %T') only ${free} GB free: stop"; exit 1; }
C=$(cat /root/eps075_commit)
cd /root/magnum.np || exit 1
if [ "$(git rev-parse HEAD)" != "$C" ]; then
    git fetch -q --depth 50 origin feature/newell-pbc && git checkout -q -f "$C" || {
        echo "$(date '+%F %T') checkout of $C failed: stop"; exit 1; }
fi
cd $B || exit 1
L=queue_eps075.log
echo "$(date '+%F %T') eps075 chain start at $(git log --oneline -1)" >> $L
if [ ! -f runs/demag_solver_N0b.md ] || ! grep -q "Convergence" runs/demag_solver_N0b.md; then
    CUDA_DEVICE=2 python check_mesh_energy.py solver N0b > runs/demag_solver_N0b.md 2> runs/demag_solver_N0b.log &
    pe=$!
fi
python run_queue.py jobs/newell_eps075.txt --gpus 0,1 --min_free_gb 3 --retries 1 --job_timeout_h 24 --gate_maxiter 0 \
    --status status_eps075.json --stop_file runs/STOP_NEW_JOBS >> $L 2>&1
rc=$?; echo "$(date '+%F %T') eps075 queue exit code $rc" >> $L
e2=0
if [ -n "$pe" ]; then wait $pe; e2=$?; fi
echo "$(date '+%F %T') solver N0b exit code $e2" >> $L
python eval_eps075.py > runs/eps075_report.md 2>> $L
e=$?; echo "$(date '+%F %T') eps075 report exit code $e" >> $L
[ $rc -eq 0 ] && [ $e -eq 0 ] && [ $e2 -eq 0 ] && touch runs/EPS075_DONE
