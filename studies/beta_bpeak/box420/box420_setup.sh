#!/bin/bash
# PLAN 4.20 box set-up (2026-11-04; from box418_setup.sh of PLAN 4.18 step 3): code of diag/sliver-coupling at the pinned commit, install, GPU check, then the box
# speed check of the plan (prof_step.py --skip 30 --plain 96 --window 0 --trace, configuration off156 of staircase
# PLAN 16 from K4_2_sat/init.pt, 1 run on each GPU). Usage: box420_setup.sh <commit>. Log: /root/progress.txt.
C=$1
P=/root/progress.txt
log() { echo "$(date -u '+%F %T') $*" >> $P; }
[ -n "$C" ] || { echo "usage: box420_setup.sh <commit>"; exit 1; }
log "setup start, commit $C"
if [ ! -d /root/magnum.np/.git ]; then
    git clone -q https://github.com/chris-birdler/magnum.np.git /root/magnum.np
    rc=$?; [ $rc -eq 0 ] || { log "clone failed rc $rc"; exit 1; }
fi
cd /root/magnum.np && git fetch -q origin && git checkout -q -f "$C"
rc=$?; [ $rc -eq 0 ] || { log "checkout of $C failed rc $rc"; exit 1; }
pip install -q --no-deps -e . > /root/pip.log 2>&1
command -v rsync > /dev/null || { apt-get update -q && apt-get install -y -q rsync; } >> /root/pip.log 2>&1   # fetch420.sh
pip install -q numpy scipy pyvista setproctitle tqdm torchdiffeq xitorch matplotlib pytest statsmodels >> /root/pip.log 2>&1
rc=$?; [ $rc -eq 0 ] || { log "pip failed rc $rc"; exit 1; }
python - >> $P 2>&1 <<'EOF'
import torch
print("torch", torch.__version__, "cuda", torch.version.cuda, "GPUs", torch.cuda.device_count(),
      [torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())])
EOF
nvidia-smi --query-gpu=index,name,memory.total,clocks.max.sm,driver_version --format=csv >> $P
log "setup done at $(git log --oneline -1 | cut -c1-60)"
# dispatcher worktree at the same commit (dispatcher.sh fetches only main); runs/ and spool/ (dispatcher makes no runs/)
[ -d /root/magnum.np-disp ] || git -C /root/magnum.np worktree add -q -f --detach /root/magnum.np-disp "$C"
mkdir -p /root/magnum.np/studies/beta_bpeak/runs /root/magnum.np/studies/beta_bpeak/spool
echo "$C" > /root/disp_commit
cp /root/magnum.np/studies/beta_bpeak/instance/dispatcher.sh /root/magnum.np/studies/beta_bpeak/instance/submit.sh /root/
log "dispatcher worktree at $(git -C /root/magnum.np-disp log --oneline -1 | cut -c1-60)"

# speed check (PLAN 4.16): the off156 line of box 11 (d24.sh), AS5 code = magnumnp of this commit
A="--d_lex 300.0 --dx_lex 3.0 --phi 0.55 --Leff_lex 18.0 --r_p 1.5 --alpha 0.1 --n_amp 1 --samples 256 --precision single --protocol virgin_asc --init_block 6.0 --relax_mode converge --relax_maxiter 80000 --relax_chunk 500 --relax_tol_E 1e-06 --relax_confirm 2 --strayfield newell_ext --strayfield_shells 4 --atol 1e-06 --relax_tol_E_abs 3e-10 --atol_cycles 1e-07 --rtol_cycles 1e-06 --cycles_per_amp 7 --min_calm 2 --max_cycles_per_amp 15 --b_tol 0.01 --init_mode saturated --predemag_h 0 --h_sat 2 --bloch_free_tries 3 --seed 2043 --sat_static --max_cycles 1"
W=/root/speed416; mkdir -p $W
[ -f /root/states/K4_2_sat/init.pt ] || { log "no /root/states/K4_2_sat/init.pt: no speed check"; exit 1; }
cd /root/magnum.np/studies/beta_bpeak
# one run on each GPU in turn (2026-11-01: box 58529564 failed on GPU 0 only); clocks, temperature and throttle
# reasons each second in smi_g<i>.csv
for g in $(nvidia-smi --query-gpu=index --format=csv,noheader); do
    log "start speed off156 GPU $g"
    nvidia-smi -i $g --query-gpu=timestamp,clocks.sm,temperature.gpu,power.draw,utilization.gpu,clocks_event_reasons.active \
        --format=csv -lms 1000 > $W/smi_g$g.csv &
    S=$!
    CUDA_DEVICE=$g TORCHINDUCTOR_CACHE_DIR=/root/ic/off156 PYTHONPATH=/root/magnum.np python ../staircase_correction/prof_step.py \
        --skip 30 --plain 96 --window 0 --trace --kcache /root/kc/off156 --json $W/off156_g$g.json -- $A \
        --init_from /root/states/K4_2_sat/init.pt --b_list 0.006 --out $W/run_off156_g$g > $W/off156_g$g.log 2>&1
    rc=$?
    kill $S
    log "end speed off156 GPU $g exit $rc"
    python - $g >> $P 2>&1 <<'EOF'
import json, sys
g = sys.argv[1]
try:
    j = json.load(open("/root/speed416/off156_g%s.json" % g))["plain"]
    r = j["ms_per_stage"] / 3.699
    # the check finds slow boxes: only a ratio above 1.10 fails; below 0.90 is a finding, not a fail (2026-11-01)
    print("SPEED gpu %s off156 %.3f ms per stage, ratio to box 11 %.3f -> %s%s" % (g, j["ms_per_stage"], r,
          "PASS" if r <= 1.10 else "FAIL", " (faster than 0.90: finding)" if r < 0.90 else ""))
except Exception as e:
    print("SPEED gpu %s no result (%r) -> FAIL" % (g, e))
EOF
done
log "DONE setup420"
