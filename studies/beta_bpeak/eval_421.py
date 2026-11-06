#!/usr/bin/env python
"""
PLAN 4.21 (gap test of V1 with the production line; D0 of Chris 2026-11-06): the readings R1 ... R5, the information
and the verdict, fixed in PLAN 4.21 before the data. Input: the downloaded run folders K<c>_<k>_v16 (runs/plan421_box)
and the phase 1 partners K<c>_<k>_sat (runs of the study worktree).
  R1  Bloch samples per run and stage (bloch_<i>_<stage>.csv, all cycles from cycle 0), split by the place of the Bloch
      points (bloch_pos_<i>_<stage>.csv of E3): a sample with a Bloch point in an interior cube (the 8 cells have the
      centre inside the sphere) and a sample with one in a cut-layer cube (a counted cube with a cell whose centre is
      outside); a sample can be in both classes. For the cut layer also the longest run of consecutive Bloch samples.
      R1 holds at K1, K3 or K4 if no run of that corner has a Bloch sample and at least 2 of its 3 runs ran all 7
      stages. K2: information (D2 (a)).
  R2  passes of the preparation and the Bloch points after each pass (summary.json); R2 holds at K1, K3 or K4 if each
      run of that corner has a clean prepared state within 3 passes; it fails for a run with 3 converged passes and a
      Bloch point after the third; a run that stopped before (no summary.json, a relaxation without a stop) is not
      evaluable. K2: information.
  R3  stages with a value: a calm tail of at least 2 cycles within 15 (summary.json calm_ok; the production rule);
      a finding if more than 2 of the 84 stages have none.
  R4  the first-run check (R4_RESULT of chain421.sh, check_sat_run.txt).
  R5  resources: passes, the hours of the preparation and of the cycles, the peak GPU and host memory (summary.json
      resources, PLAN 4.21 K4). The times include E3 (an upper value of the production time).
  Information: ln P and β per run (eval_corners_sat.outputs, plus ln P at 150 mT) against the phase 1 partner, paired
  by seed per corner (n 3, mean ± SE, 95 % t-interval; factors on P). This is the whole surface model against the
  voxel surface, not V1 alone.
  Verdict: "the production line with V1 passes" if R1 and R2 hold at K1, K3 and K4, at most 2 of the 84 stages have no
  value and no R4 FAIL; "fails" if a run of K1, K3 or K4 has a Bloch sample or no clean prepared state after 3 passes;
  else "not determined".

    python eval_421.py runs/plan421_box --phase1 ../../../magnum.np-beta_bpeak/studies/beta_bpeak/runs \
        --out results/plan421
"""
import argparse
import json
import math
import pathlib
import warnings

import numpy as np
from scipy import stats

import analyze
import analyze_design as AD
import domains
import eval_corners_nx as EN
import eval_corners_sat as ES
import units as U

CORNERS = ("K1", "K2", "K3", "K4")
STRICT = ("K1", "K3", "K4")
N_STAGES = 7
MAX_NO_VALUE = 2


def read_rows(path):
    """Numeric rows of a csv file with # comments ([] if missing or empty)."""
    if not path.exists():
        return []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)              # a file with the header only
        a = np.loadtxt(path, delimiter=",", comments="#", ndmin=2)
    return a.tolist() if a.size else []


def cube_classes(geo):
    """(interior, cut layer) masks of the cube corners [N,N,N]: the counted cubes (8 cells of one particle,
    run_loops.bloch_mask_t) with all 8 cells inside the sphere, and the counted cubes with a cell outside."""
    ids, inside = geo["ids"], geo["inside"]
    same, allin = ids >= 0, inside.copy()
    for s in ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)):
        same &= np.roll(ids, (-s[0], -s[1], -s[2]), axis=(0, 1, 2)) == ids
        allin &= np.roll(inside, (-s[0], -s[1], -s[2]), axis=(0, 1, 2))
    return same & allin, same & ~allin


def longest_run(flags):
    best = cur = 0
    for f in flags:
        cur = cur + 1 if f else 0
        best = max(best, cur)
    return best


def r1_run(d, geo):
    """Per stage: samples, Bloch samples, interior and cut-layer Bloch samples, the longest cut-layer run."""
    inner, cutl = cube_classes(geo)
    out = []
    for f in sorted(d.glob("bloch_[0-9][0-9]_*.csv")):
        if f.name.startswith("bloch_pos"):
            continue
        si, stage = f.stem.split("_", 2)[1:]
        rows = read_rows(f)
        keys = [(int(r[0]), int(r[1])) for r in rows]
        nb = {k: r[2] + r[3] for k, r in zip(keys, rows)}
        cls = {}
        for r in read_rows(d / ("bloch_pos_%s_%s.csv" % (si, stage))):
            k = (int(r[0]), int(r[1]))
            i, j, kk = int(r[2]), int(r[3]), int(r[4])
            c = cls.setdefault(k, [False, False])
            c[0] |= bool(inner[i, j, kk])
            c[1] |= bool(cutl[i, j, kk])
        cut_flags = [cls.get(k, [False, False])[1] for k in keys]
        out.append({"stage": stage, "samples": len(keys), "bloch_samples": int(sum(v > 0 for v in nb.values())),
                    "interior_samples": int(sum(v[0] for v in cls.values())),
                    "cut_layer_samples": int(sum(v[1] for v in cls.values())),
                    "longest_cut_layer_run": longest_run(cut_flags),
                    "bloch_pos_missing": bool(sum(v > 0 for v in nb.values()) and not cls)})
    return out


def run_reading(d):
    """The readings of one run folder. Without summary.json (a run that stopped in the preparation or in stage 0) the
    run has 0 stages and its preparation has not ended; R1 comes from the Bloch files."""
    sf = d / "summary.json"
    s = json.loads(sf.read_text()) if sf.exists() else {}
    i = s.get("init", {})
    passes = i.get("predemag_passes") or []
    nb = [len(p.get("bloch_points") or []) for p in passes]            # the Bloch points after each pass
    conv = [bool(p.get("relax_converged", True)) for p in passes]
    prep_fail = len(passes) >= 3 and all(conv) and nb[-1] > 0            # R2 fails: 3 converged passes, not clean
    prep_clean = bool(passes) and all(conv) and nb[-1] == 0 and bool(i.get("bloch_free", True))
    prep_h = (i.get("wall_s", 0.0) + sum(p.get("wall_s", 0.0) for p in passes)) / 3600.0
    stages = s.get("stages", [])
    cyc_h = sum(c.get("wall_s", 0.0) for st in stages for c in st.get("cycles", [])) / 3600.0
    r4 = (d / "R4_RESULT").read_text().strip() if (d / "R4_RESULT").exists() else None
    res = [st.get("resources") for st in stages if st.get("resources")] + [i.get("resources")]
    res = [x for x in res if x]
    gpu = [x["gpu_peak_reserved_GiB"] for x in res if x.get("gpu_peak_reserved_GiB") is not None]
    host = [x["host_peak_rss_GiB"] for x in res if x.get("host_peak_rss_GiB") is not None]
    return {"run": d.name, "corner": d.name[:2], "k": int(d.name.split("_")[1]), "done": (d / "DONE").exists(),
            "n_stages": len(stages), "complete_7": len(stages) == N_STAGES and all(st.get("complete") for st in stages),
            "calm_ok": [bool(st.get("calm_ok")) for st in stages], "passes": len(passes), "bp_after_pass": nb,
            "prep_clean": prep_clean, "prep_fail": prep_fail, "prep_evaluable": prep_clean or prep_fail,
            "prep_h": prep_h, "cycles_h": cyc_h, "R4": r4,
            "gpu_peak_GiB": max(gpu) if gpu else None, "host_peak_GiB": max(host) if host else None,
            "R1": r1_run(d, domains.run_geometry(d))}


def pair_stat(diff):
    """(mean, SE, CI low, CI high, n) of paired differences with the 95 % t-interval; None for fewer than 2."""
    v = np.array([x for x in diff if np.isfinite(x)])
    if len(v) < 2:
        return None
    se = v.std(ddof=1) / math.sqrt(len(v))
    h = stats.t.ppf(0.975, len(v) - 1) * se
    return float(v.mean()), float(se), float(v.mean() - h), float(v.mean() + h), len(v)


def information(runs_dir, phase1, names, lnF):
    """ln P and β of the test runs and of their phase 1 partners; paired differences per corner."""
    rows = {}
    for n in names:
        for key, d in (("v16", runs_dir / n), ("p1", phase1 / n.replace("_v16", "_sat"))):
            if not (d / "DONE").exists():
                continue
            lnw = EN.lnw(d)
            o = ES.outputs(lnw, AD.run_record(d), lnF)
            o["lnP@150mT"] = float(lnw[ES.AMPS.index(150.0)] + lnF)
            rows.setdefault(n, {})[key] = o
    keys = ES.OUTPUTS + ["lnP@150mT"]
    out = {}
    for c in CORNERS:
        out[c] = {}
        for k in keys:
            diff = [rows[n]["v16"][k] - rows[n]["p1"][k] for n in names
                    if n.startswith(c) and "v16" in rows.get(n, {}) and "p1" in rows.get(n, {})]
            ps = pair_stat(diff)
            out[c][k] = ps
    return rows, out


def verdict(R):
    """The verdict of PLAN 4.21 and its reasons."""
    why, fail = [], False
    for c in STRICT:
        rs = [r for r in R if r["corner"] == c]
        bloch = [r["run"] for r in rs if any(st["bloch_samples"] for st in r["R1"])]
        noprep = [r["run"] for r in rs if r["prep_fail"]]
        if bloch:
            fail = True
            why.append("%s: Bloch samples in %s" % (c, ", ".join(bloch)))
        if noprep:
            fail = True
            why.append("%s: no clean prepared state after 3 passes in %s" % (c, ", ".join(noprep)))
        n7 = sum(r["complete_7"] for r in rs)
        if n7 < 2:
            why.append("%s: only %d runs ran all 7 stages" % (c, n7))
        unev = [r["run"] for r in rs if not r["prep_evaluable"]]
        if unev:
            why.append("%s: R2 not evaluable in %s" % (c, ", ".join(unev)))
    no_value = sum(not x for r in R for x in r["calm_ok"]) + sum(N_STAGES - r["n_stages"] for r in R)
    if no_value > MAX_NO_VALUE:
        why.append("%d stages without a value (more than %d)" % (no_value, MAX_NO_VALUE))
    r4f = [r["run"] for r in R if r["R4"] != "PASS" and r["done"]]
    if r4f:
        why.append("R4 not PASS: %s" % ", ".join(r4f))
    if fail:
        return "the production line with V1 fails", why, no_value
    return ("the production line with V1 passes" if not why else "not determined"), why, no_value


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs")
    ap.add_argument("--phase1", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    runs_dir, phase1 = pathlib.Path(a.runs), pathlib.Path(a.phase1)
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    u = U.units()
    lnF = math.log(u["Kd"] * U.F_REL_DEFAULT * u["f_M"] / 1e6)           # ln w -> ln P (W/cm³), as eval_corners_sat
    names = sorted((d.name for d in runs_dir.iterdir() if d.is_dir() and d.name.endswith("_v16")),
                   key=lambda n: (n[:2], int(n.split("_")[1])))
    R = [run_reading(runs_dir / n) for n in names]
    v, why, no_value = verdict(R)
    rows, info = information(runs_dir, phase1, names, lnF)
    res = {"verdict": v, "reasons": why, "stages_without_value": no_value, "runs": R, "info_rows": rows, "info": info}
    pathlib.Path(a.out + ".json").write_text(json.dumps(res, indent=1))
    L = ["# PLAN 4.21: gap test of V1 with the production line (eval_421.py)", "",
         "**Verdict: %s.**%s" % (v, (" Reasons: " + "; ".join(why) + ".") if why else ""), "",
         "| run | stages | passes (BPs after each) | R4 | Bloch samples per stage (interior / cut layer, longest cut run) |"
         " prep / cycles (h) | GPU / host peak (GiB) |", "|---|---|---|---|---|---|---|"]
    for r in R:
        b = "; ".join("%s %d (%d / %d, %d)" % (st["stage"], st["bloch_samples"], st["interior_samples"],
                                              st["cut_layer_samples"], st["longest_cut_layer_run"])
                      for st in r["R1"] if st["bloch_samples"]) or "0"
        L.append("| %s | %d (%d with a value) | %d (%s) | %s | %s | %.2f / %.2f | %s / %s |" % (
            r["run"], r["n_stages"], sum(r["calm_ok"]), r["passes"], " / ".join(str(x) for x in r["bp_after_pass"]),
            r["R4"], b, r["prep_h"], r["cycles_h"],
            "%.1f" % r["gpu_peak_GiB"] if r["gpu_peak_GiB"] is not None else "–",
            "%.1f" % r["host_peak_GiB"] if r["host_peak_GiB"] is not None else "–"))
    L += ["", "## Information: against the phase 1 partner (paired by seed; factor on P = exp(mean Δ ln P))", "",
          "| corner | quantity | mean Δ ± SE | 95 % t-interval | factor on P | n |", "|---|---|---|---|---|---|"]
    for c in CORNERS:
        for k, s in info[c].items():
            if s is None:
                continue
            fac = ("%.3f [%.3f, %.3f]" % (math.exp(s[0]), math.exp(s[2]), math.exp(s[3]))) if k.startswith("lnP") else "–"
            L.append("| %s | %s | %+.3f ± %.3f | [%+.3f, %+.3f] | %s | %d |" % (c, k, s[0], s[1], s[2], s[3], fac, s[4]))
    pathlib.Path(a.out + ".md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
