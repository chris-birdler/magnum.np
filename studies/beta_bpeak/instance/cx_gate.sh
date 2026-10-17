#!/bin/bash
# Smoke check of PLAN 4.13 (part 5, step 2) as an early warning. Chris 2026-10-18: idle GPUs waste money, so the corner
# chain starts at once and does not wait for this check. The script waits for the smoke run runs/SMOKE_nx_K1, checks
# its pass criteria before the post-processing, runs the post-processing once on it, and writes the result to its log;
# a failure goes to runs/ALERT (the corner chain goes on).
#   start: setsid nohup bash /root/cx_gate.sh > /root/cx_gate.log 2>&1 < /dev/null &
W=/root/magnum.np-cx/studies/beta_bpeak
R=/root/magnum.np/studies/beta_bpeak/runs
S=$R/SMOKE_nx_K1
log() { echo "$(date '+%F %T') $*"; }
log "gate start: wait for $S"
while [ ! -f $S/DONE ] && pgrep -f "[m]agnumnp run_loops" > /dev/null; do sleep 60; done
sleep 30
cd $W || exit 1
python - "$S" <<'EOF'
import json, pathlib, re, sys
import numpy as np
d = pathlib.Path(sys.argv[1])
res = [("DONE", (d / "DONE").exists())]
if res[0][1]:
    s = json.loads((d / "summary.json").read_text())
    st = [x for x in s["stages"] if x.get("measure")][0]
    n = len(st["cycles"])
    cyc = sorted({int(re.search(r"_c(\d+)_ph", f.name).group(1)) for f in d.glob("m_amp00_c*_ph*.vti")})
    a = np.load(d / "acc_00_amp00.npz") if (d / "acc_00_amp00.npz").exists() else None
    res += [("relaxation converged", bool(s["init"].get("relax_converged"))),
            ("stage longer than the minimum (cleanup ran)", n >= 5),
            ("snapshots only of the last 2 cycles", cyc == [n - 2, n - 1]
             and len(list(d.glob("m_amp00_c*_ph*.vti"))) == 32),
            ("acc_00 with mk2 (7 harmonics)", a is not None and "mk2" in a.files and a["mk2"].shape[1] == 7),
            ("ckpt_00 written", (d / "ckpt_00_amp00.pt").exists()),
            ("energy terms per cycle and at the stage end", all("E_start" in c for c in st["cycles"]) and "E_end" in st)]
    print("cycles %d, snapshot cycles %s, relaxation iterations %s" % (n, cyc, s["init"].get("relax_iterations")))
for k, v in res:
    print("%-48s %s" % (k, "PASS" if v else "FAIL"))
sys.exit(0 if all(v for _, v in res) else 1)
EOF
c1=$?
log "file checks exit code $c1"
if [ $c1 -eq 0 ]; then
    touch $R/SMOKE_PP_STOP
    python instance_postproc.py --runs runs --stop runs/SMOKE_PP_STOP --once --workers 1 --pattern '^SMOKE_nx_K1$'
    log "post-processing exit code $?"
    [ -f $S/POSTPROC_DONE ] && [ ! -f $S/POSTPROC_FAILED ] && c2=0 || c2=1
    log "post-processing check $c2"
else
    c2=1
fi
if [ $c1 -eq 0 ] && [ $c2 -eq 0 ]; then
    log "smoke PASS"
else
    echo "$(date '+%F %T') smoke test of PLAN 4.13 FAILED (see /root/cx_gate.log); the corner chain goes on" >> $R/ALERT
    log "smoke FAIL (corner chain goes on)"
fi
