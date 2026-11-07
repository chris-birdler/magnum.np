#!/usr/bin/env python
"""
Numbers of gate 1 of PLAN 4.21 again with V1 (D0 of Chris 2026-11-06), written to results/plan421_gate1_v1.json:
  runs     per K5 run of runs/plan421_gate1_v1 (L_eff/l_ex 6 and 12): the tie counts, the passes of the preparation and
           the Bloch points after each, the preparation half cycles that ended, the hours of the set-up (init.wall_s:
           from the start of the solver to the end of the first diagnostics; the saturated start has no relaxation
           there), of each pass of the preparation and of the stages, the Bloch samples per stage, and the Bessel
           numbers of acc.json (bessel_max, the share of the cells with sum_k p_k above 1.001 p) with the identity
           sum_k R_k + R_rest = 1 (true by the definition of R_rest)
  control  the same Bessel numbers for the voxel run K5v of the first gate 1 (runs/plan421_gate1: the same K5 line
           without the fraction surface, the correction and a tie), from acc_reduce on a copy of its accumulators (the
           original folder is not changed)
  phase1   bessel_max of all stages of the corner runs K<c>_<k>_sat of the study worktree (4.15 corner line, voxel, no
           tie; acc.json): the number of stages above 1.01, in all and at amp00, and the largest value
    python eval_gate421_v1.py runs/plan421_gate1_v1 runs/plan421_gate1 \
        ../../../magnum.np-beta_bpeak/studies/beta_bpeak/runs
"""
import json
import pathlib
import shutil
import sys
import tempfile

import acc_reduce
from eval_gate421 import half_cycles

HERE = pathlib.Path(__file__).resolve().parent


def bessel(acc):
    return {st["stage"]: {"bessel_max": st["harm"]["bessel_max"],
                          "share_above_1.001": st["harm"]["bessel_frac_over_1e-3"],
                          "identity_dev": sum(st["harm"]["R_k"]) + st["harm"]["R_rest"] - 1.0} for st in acc["stages"]}


def bloch_samples(f):
    rows = [l.split(",") for l in f.read_text().splitlines() if l and not l.startswith("#")]
    return {"samples": len(rows), "with_bloch_point": sum(int(r[2]) + int(r[3]) > 0 for r in rows)}


def run_numbers(d):
    cfg = json.loads((d / "config.json").read_text())
    s = json.loads((d / "summary.json").read_text())
    pp = s["init"]["predemag_passes"]
    return {"Leff_lex": cfg["Leff_lex"], "tie": cfg["diag"]["tie"], "passes": len(pp),
            "bloch_points_after_pass": [len(p["bloch_points"]) for p in pp],
            "half_cycles": half_cycles(d.parent / (d.name + ".log")),
            "setup_h": s["init"]["wall_s"] / 3600.0, "pass_h": [p["wall_s"] / 3600.0 for p in pp],
            "stage_h": [sum(c["wall_s"] for c in st["cycles"]) / 3600.0 for st in s["stages"]],
            "bloch": {f.name: bloch_samples(f) for f in sorted(d.glob("bloch_[0-9]*_amp*.csv"))},
            "bessel": bessel(json.loads((d / "acc.json").read_text()))}


def control(d):
    with tempfile.TemporaryDirectory() as t:
        c = pathlib.Path(t) / d.name
        c.mkdir()
        for f in list(d.glob("acc_*_*.npz")) + list(d.glob("samples_*.csv")) + [d / "config.json", d / "geometry.npz",
                                                                                d / "summary.json"]:
            shutil.copy2(f, c / f.name)
        acc_reduce.run(c)
        return bessel(json.loads((c / "acc.json").read_text()))


def phase1(runs):
    v = [(f.parent.name, st["stage"], st["harm"]["bessel_max"]) for f in sorted(runs.glob("K*_sat/acc.json"))
         for st in json.loads(f.read_text())["stages"] if st.get("harm")]
    a0 = [x for x in v if x[1] == "amp00"]
    return {"n_runs": len({x[0] for x in v}), "n_stages": len(v), "stages_above_1.01": sum(x[2] > 1.01 for x in v),
            "amp00_stages": len(a0), "amp00_above_1.01": sum(x[2] > 1.01 for x in a0), "max": max(x[2] for x in v)}


if __name__ == "__main__":
    v1, g1, p1 = (pathlib.Path(a).resolve() for a in sys.argv[1:4])
    res = {"runs": {n: run_numbers(v1 / n) for n in ("K5f", "K5f12")}, "control_K5v": control(g1 / "K5v"),
           "phase1": phase1(p1)}
    (HERE / "results" / "plan421_gate1_v1.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))
