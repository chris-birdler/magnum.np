#!/usr/bin/env python
"""
Check C5 of PLAN 4.15 on one finished "_sat" run folder, before its post-processing (instance_postproc deletes the
snapshots and the accumulator files): DONE; summary.json with the static record and the passes; 7 complete stages and
7 stage states; the energy terms (E_start) in each cycle record; snapshots only of the last 2 cycles (16 phases) at
amp00, amp03 and amp06; accumulators with the harmonics (mk2, 7 columns) and the sum check p_box / p_dis = 1.
Run on the instance (2026-10-20 on K1_0_sat: all PASS, PROTOKOLL 7.134):

    python3 check_sat_run.py runs/K1_0_sat
"""
import json, pathlib, re, sys
import numpy as np

d = pathlib.Path(sys.argv[1])
ok = True
def chk(name, cond, info=""):
    global ok
    ok &= bool(cond)
    print("%-4s %s %s" % ("PASS" if cond else "FAIL", name, info))

chk("DONE", (d / "DONE").exists())
s = json.loads((d / "summary.json").read_text())
i = s["init"]
pd = i.get("predemag", {})
st = pd.get("static")
chk("static record", isinstance(st, dict) and st.get("relax_converged") is True,
    "iter %s, criterion %s, alignment %s, n_bloch %s, wall %.0f s" % (
        st.get("relax_iterations"), st.get("criterion"), st.get("alignment"), st.get("n_bloch"), st.get("wall_s", 0)))
pp = i.get("predemag_passes", [])
chk("passes", len(pp) >= 1 and all("static" in p for p in pp),
    "%d pass(es), BPs after each pass %s, bloch_free %s" % (len(pp), [p.get("bloch_points") for p in pp], i.get("bloch_free")))
stages = s["stages"]
chk("7 stages complete", len(stages) == 7 and all(x.get("complete") for x in stages),
    "h_amp %s" % [round(x["h_amp"], 4) for x in stages])
ck = sorted(d.glob("ckpt_*_amp*.pt"))
chk("7 stage states", len(ck) == 7, " ".join(f.name for f in ck))
# energy terms in each cycle record
miss = [(x["name"], c.get("cycle")) for x in stages for c in x["cycles"] if not c.get("E_start")]
e0 = stages[0]["cycles"][0]["E_start"]
chk("E_start in each cycle", not miss, "missing %s; example keys %s" % (miss[:5], list(e0.keys()) if isinstance(e0, dict) else e0))
# snapshots: only the last 2 cycles of each snapshot stage, 16 phases each
snaps = {}
for f in d.glob("m_*_c*_ph*.vti"):
    m = re.match(r"m_(amp\d+)_c(\d+)_ph(\d+)\.vti", f.name)
    snaps.setdefault(m.group(1), {}).setdefault(int(m.group(2)), set()).add(m.group(3))
for x in stages:
    tag = x["name"].split("_")[-1]
    if tag in snaps:
        cyc = sorted(c["cycle"] for c in x["cycles"])
        got = snaps[tag]
        chk("snapshots %s" % tag, sorted(got) == cyc[-2:] and all(len(v) == 16 for v in got.values()),
            "cycles in record %s, snapshot cycles %s, phases %s" % (cyc, sorted(got), [len(v) for v in got.values()]))
chk("snapshot stages", sorted(snaps) == ["amp00", "amp03", "amp06"], "%s" % sorted(snaps))
# accumulators with the harmonics
for f in sorted(d.glob("acc_*_amp*.npz")):
    z = np.load(f)
    r = float(z["p_box"]) / float(z["p_dis_samples"]) if "p_dis_samples" in z else float("nan")
    chk("acc %s" % f.name, "mk2" in z and z["mk2"].ndim == 2 and z["mk2"].shape[1] == 7 and abs(r - 1) < 1e-4,
        "mk2 %s, cycles kept %s, p_box/p_dis %.6f" % (z["mk2"].shape if "mk2" in z else None, list(z["cycles_kept"]), r))
print("stage summary: name, b_peak last cycle, n cycles, calm_ok, BP end, steady")
for x in stages:
    c = x["cycles"][-1]
    print("  %s  b %.4f  n %d  calm %s  BP %s  steady %s" % (x["name"], c["b_peak"], len(x["cycles"]), x.get("calm_ok"),
                                                         x.get("bloch_points_end"), x.get("steady")))
print("C5 data part:", "PASS" if ok else "FAIL")
