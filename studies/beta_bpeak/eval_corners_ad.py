#!/usr/bin/env python
"""
Corner evaluation with three arms (Chris 2026-10-19; method agreed with Chris before the evaluation, PROTOKOLL 7.119).
  old  the 12 corner runs of the old main study (DemagFieldPBC, atol 1e-5; eval_corners_nx.old_dir: the _mc2 rerun
       where it exists)
  nx   the corner test runs <run>_nx (setup of PLAN 4.13; relaxed random start, with Bloch points; PROTOKOLL 7.112)
  ad   option C <run>_ad (setup of PLAN 4.13; AC from 300 mT, at most 3 attempts; no Bloch point in the prepared state;
       PROTOKOLL 7.118)
Paired per microstructure (n 12: 4 corners × 3).
β per run = slope of ln w over ln B at the nominal amplitudes (main-fit definition, steady-tail rule ≥ 2 calm cycles),
in three windows: full 9 … 150 mT (7 amplitudes, at least 6 with a value), low 9 … 50 mT (all 4), high 70 … 150 mT
(all 3; information only). A comparison fits both runs on the amplitudes with a value in both; the per-corner tables
fit all three arms on the amplitudes with a value in all three.
Readings (rule of eval_corners_nx, fixed before the evaluation) for Δβ(ad − old) and Δβ(ad − nx) in the windows full
and low:
  "robust"          the 90 % CI of the mean Δβ lies inside ± 0.10, and no corner has |mean Δβ| > 0.20 (n 3 per corner);
  "shifted"         the 90 % CI of the mean Δβ excludes 0 and |mean Δβ| > 0.10;
  "not determined"  all other cases, or fewer than 9 pairs.
Information: β per corner and arm (mean ± SE, n 3); levers β(Ki) − β(K1) per arm and window (SE from the two corner
SE); P per corner, arm and amplitude (W/cm³, geometric mean over the 3 microstructures) and the factors
exp(mean Δ ln P) (ad/old, ad/nx; SE of Δ ln P, n 12); stages with a Bloch point; fields (R1, hottest 10 %, outer shell,
jumps, β_f; eval_corners_fields.record); quality of the ad runs (eval_corners_nx.info_lines). The numbers for the
figures go to results/corners_ad_data.json (plot_corners_ad.py).

    python eval_corners_ad.py [--runs DIR] [--json FILE] > results/corners_ad_report.md
"""
import argparse
import json
import math
import pathlib

import numpy as np

import analyze
import analyze_design as AD
import conv_jobs as CJ
import eval_corners_fields as EF
import eval_corners_nx as EN
import units as U

HERE = pathlib.Path(__file__).resolve().parent
ARMS = ("old", "nx", "ad")
AMPS = np.asarray(AD.AMPS_MT, float)
WINDOWS = {"full": (tuple(range(0, 7)), 6), "low": (tuple(range(0, 4)), 4), "high": (tuple(range(4, 7)), 3)}
WIN_LABEL = {"full": "9 … 150 mT", "low": "9 … 50 mT", "high": "70 … 150 mT"}
KEY = {"β robust to the setup change": "robust", "β shifted": "shifted", "not determined": "not determined"}
RUNS = EF.RUNS


def arm_dir(root, n, arm):
    return EN.old_dir(root, n) if arm == "old" else root / (n + "_" + arm)


def window_beta(lnws, win):
    """β of each ln w array over the amplitudes of window win with a value in all arrays (at least the window minimum);
    NaN for all if fewer."""
    idx, nmin = WINDOWS[win]
    idx = np.array(idx)
    a = np.array([np.asarray(v, float)[idx] for v in lnws])
    m = np.all(np.isfinite(a), axis=0)
    if m.sum() < nmin:
        return [float("nan")] * len(lnws)
    x = np.log(AMPS[idx][m])
    return [float(np.polyfit(x, v[m], 1)[0]) for v in a]


def corner_of(n):
    return n.split("_")[0]


def bp_stages(root, n):
    """Stages (indices 0 … 6) with a Bloch point, per arm. ad: count at each sample (bloch_<i>_<stage>.csv); nx: census
    at the start or the end of the stage (census_bloch_nx_v2.json); old: census of the start state and the end state
    only (census_bloch_old_v2.json) -> ([], start, end)."""
    out = {}
    d = root / (n + "_ad")
    s = []
    for i in range(7):
        f = sorted(d.glob("bloch_%02d_*.csv" % i))
        if f:
            a = np.loadtxt(f[0], delimiter=",", comments="#", ndmin=2)
            if len(a) and (a[:, 2] + a[:, 3]).max() > 0:
                s.append(i)
    out["ad"] = s
    cn = root / "census_bloch_nx_v2.json"
    if cn.exists():
        c = json.loads(cn.read_text()).get(n + "_nx", {})
        keys = ["init.pt"] + ["ckpt_%02d_amp%02d.pt" % (i, i) for i in range(7)]
        nb = [len(c[k]["bloch_points"]) if k in c else None for k in keys]
        out["nx"] = [i for i in range(7) if (nb[i] or 0) > 0 or (nb[i + 1] or 0) > 0]
    co = root / "census_bloch_old_v2.json"
    if co.exists():
        c = json.loads(co.read_text()).get(n, {})
        out["old"] = (len(c["init.pt"]["bloch_points"]) if "init.pt" in c else None,
                      len(c["checkpoint.pt"]["bloch_points"]) if "checkpoint.pt" in c else None)
    return out


def fstat(s, f="%+.3f"):
    return ("%s ± %.3f" % (f % s[0], s[1])) if s else "–"


def main(root=None, json_out=None):
    root = pathlib.Path(root) if root else HERE / "runs"
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    u = U.units()
    lnF = math.log(u["Kd"] * U.F_REL_DEFAULT * u["f_M"] / 1e6)          # ln w -> ln P (W/cm³)
    W = {n: {a: EN.lnw(arm_dir(root, n, a)) for a in ARMS} for n in RUNS}
    data = {"runs": RUNS, "amps_mT": AMPS.tolist(), "lnP": {n: {a: (W[n][a] + lnF).tolist() for a in ARMS} for n in RUNS}}
    L = ["# Corner evaluation with three arms: old main study, new solver with Bloch points, new solver without Bloch "
         "points (eval_corners_ad.py)", "",
         "old: DemagFieldPBC, atol 1e-5 (the _mc2 rerun where it exists); nx: setup of PLAN 4.13, relaxed random start "
         "(7.112); ad: setup of PLAN 4.13, option C (7.118). β = slope of ln w over ln B (main-fit definition). Windows: "
         + "; ".join("%s %s" % (k, WIN_LABEL[k]) for k in WINDOWS) + ". Mean ± SE; n 3 per corner, n 12 paired.", ""]

    # --- per run β, all three arms on common amplitudes --------------------------------------------------------------
    B3 = {w: {n: dict(zip(ARMS, window_beta([W[n][a] for a in ARMS], w))) for n in RUNS} for w in WINDOWS}
    data["beta"] = B3
    L += ["## β per run (all three arms on the same amplitudes)", "",
          "| run (old partner) | " + " | ".join("%s old / nx / ad" % WIN_LABEL[w] for w in WINDOWS) + " |",
          "|---|" + "---|" * len(WINDOWS)]
    for n in RUNS:
        L.append("| %s (%s) | " % (n, EN.old_dir(root, n).name) + " | ".join(
            " / ".join("%.3f" % B3[w][n][a] if np.isfinite(B3[w][n][a]) else "–" for a in ARMS) for w in WINDOWS) + " |")

    # --- per corner β and levers -------------------------------------------------------------------------------------
    CS = {w: {c: {a: EN.stat([B3[w][n][a] for n in RUNS if corner_of(n) == c]) for a in ARMS} for c in CJ.CORNERS_NX}
          for w in WINDOWS}
    data["corner_beta"] = {w: {c: {a: list(CS[w][c][a][:2]) if CS[w][c][a] else None for a in ARMS}
                               for c in CJ.CORNERS_NX} for w in WINDOWS}
    for w in WINDOWS:
        L += ["", "## β per corner, %s (mean ± SE, n 3)" % WIN_LABEL[w], "", "| corner | old | nx | ad |",
              "|---|---|---|---|"]
        for c in CJ.CORNERS_NX:
            L.append("| %s | %s |" % (c, " | ".join(fstat(CS[w][c][a], "%.3f") for a in ARMS)))
        lev = []
        for c in CJ.CORNERS_NX[1:]:
            cells = []
            for a in ARMS:
                x, y = CS[w][c][a], CS[w]["K1"][a]
                cells.append("%+.3f ± %.3f" % (x[0] - y[0], math.hypot(x[1], y[1])) if x and y else "–")
            lev.append("| %s − K1 | %s |" % (c, " | ".join(cells)))
        L += ["", "Levers β(Ki) − β(K1), %s:" % WIN_LABEL[w], "", "| lever | old | nx | ad |", "|---|---|---|---|"] + lev

    # --- paired comparisons and readings ------------------------------------------------------------------------------
    data["paired"] = {}
    L += ["", "## Paired Δβ and readings (rule fixed before the evaluation)", "",
          "| comparison | window | mean Δβ (SE, n; 90 % CI) | per corner (n 3) | reading |", "|---|---|---|---|---|"]
    for a, b in (("ad", "old"), ("ad", "nx"), ("nx", "old")):
        for w in WINDOWS:
            db = {n: (lambda p: p[1] - p[0])(window_beta([W[n][b], W[n][a]], w)) for n in RUNS}
            by_c = {c: [db[n] for n in RUNS if corner_of(n) == c] for c in CJ.CORNERS_NX}
            key, s = EN.reading(list(db.values()), by_c)
            rd = KEY[key] if w != "high" and (a, b) != ("nx", "old") else "(information)"
            data["paired"]["%s-%s_%s" % (a, b, w)] = {"runs": db, "stat": list(s) if s else None}
            L.append("| %s − %s | %s | %s | %s | %s |" % (
                a, b, WIN_LABEL[w],
                "%+.3f ± %.3f (n %d; %+.3f … %+.3f)" % (s[0], s[1], s[4], s[2], s[3]) if s else "–",
                "; ".join("%s %s" % (c, fstat(EN.stat(v))) for c, v in by_c.items()), rd))

    # --- loss density P and factors -------------------------------------------------------------------------------------
    L += ["", "## Loss density P (W/cm³; geometric mean over the 3 microstructures of a corner)", "",
          "| corner | arm | " + " | ".join("%g mT" % b for b in AMPS) + " |", "|---|---|" + "---|" * len(AMPS)]
    for c in CJ.CORNERS_NX:
        for a in ARMS:
            v = np.nanmean([W[n][a] for n in RUNS if corner_of(n) == c], 0) + lnF
            L.append("| %s | %s | %s |" % (c, a, " | ".join("%.3g" % math.exp(x) if np.isfinite(x) else "–" for x in v)))
    data["factor"] = {}
    for a, b in (("ad", "old"), ("ad", "nx"), ("nx", "old")):
        cells = []
        data["factor"]["%s/%s" % (a, b)] = []
        for k, amp in enumerate(AMPS):
            s = EN.stat([W[n][a][k] - W[n][b][k] for n in RUNS])
            data["factor"]["%s/%s" % (a, b)].append(list(s) if s else None)
            cells.append("%g mT × %.3f (Δ ln P %+.3f ± %.3f, n %d)" % (amp, math.exp(s[0]), s[0], s[1], s[4]) if s
                         else "%g mT –" % amp)
        L += ["", "Factor on P, %s / %s (exp of the mean paired Δ ln P): " % (a, b) + "; ".join(cells) + "."]

    # --- Bloch points -----------------------------------------------------------------------------------------------------
    L += ["", "## Stages with a Bloch point", "",
          "ad: count at each sample; nx: census at the start or the end of a stage (7.112); old: census of the start "
          "state and the end state only (no stage information; the _mc2 reruns have no census).", "",
          "| run | ad: stages with a BP | nx: stages with a BP at start or end | old: BPs at start / end |",
          "|---|---|---|---|"]
    data["bp"] = {}
    for n in RUNS:
        bp = bp_stages(root, n)
        data["bp"][n] = bp
        nm = lambda s: ", ".join("%g" % AMPS[i] for i in s) + " mT" if s else "none"
        L.append("| %s | %s | %s | %s |" % (n, nm(bp.get("ad", [])), nm(bp.get("nx", [])) if "nx" in bp else "–",
                                            "%s / %s" % bp["old"] if "old" in bp else "–"))

    # --- fields ---------------------------------------------------------------------------------------------------------
    D = {n: {a: EF.record(arm_dir(root, n, a), lnF) for a in ARMS} for n in RUNS}
    for key, title, f in (("R1", "Fundamental share R1 of the dissipation", "%.2f"),
                          ("top", "Share of the loss in the hottest 10 % of the cells", "%.2f"),
                          ("shell", "Share of the loss in the outer shell r ≥ 0.8 (49 % of the volume)", "%.2f"),
                          ("jump", "Share of the dissipation in jumps (KH 8, 5 sigma)", "%.3f")):
        L += ["", "## " + title + " (mean over the 3 microstructures; 9 … 150 mT)", "", "| corner | old | nx | ad |",
              "|---|---|---|---|"]
        for c in CJ.CORNERS_NX:
            L.append("| %s | %s |" % (c, " | ".join(EF.fmt(np.nanmean([D[n][a][key] for n in RUNS if corner_of(n) == c], 0), f)
                                                for a in ARMS)))
        for a, b in (("ad", "old"), ("ad", "nx")):
            d = np.array([D[n][a][key] - D[n][b][key] for n in RUNS])
            m = np.nanmean(d, 0)
            se = np.nanstd(d, 0, ddof=1) / np.sqrt(np.sum(np.isfinite(d), 0))
            L += ["", "Paired %s − %s (mean ± SE, n 12): " % (a, b) + "; ".join(
                "%g mT %+.3f ± %.3f" % (amp, x, e) for amp, x, e in zip(AMPS, m, se)) + "."]
    L += ["", "## Exponents over 9 … 150 mT: β and the exponent of the fundamental channel β_f = β + d ln R1 / d ln B",
          "", "| corner | β old / nx / ad | β_f old / nx / ad |", "|---|---|---|"]
    for c in CJ.CORNERS_NX:
        rr = [D[n] for n in RUNS if corner_of(n) == c]
        L.append("| %s | %s | %s |" % (c, " / ".join("%.2f" % np.nanmean([r[a]["beta"] for r in rr]) for a in ARMS),
                                       " / ".join("%.2f" % np.nanmean([r[a]["beta_f"] for r in rr]) for a in ARMS)))

    # --- quality of the ad runs -----------------------------------------------------------------------------------------
    L += [x.replace("## Information (PLAN 4.13)", "## Quality of the ad runs (PLAN 4.13 S5, S9; Δβ here: ad − old, "
                    "9 … 150 mT)") for x in EN.info_lines(root, RUNS, {n: (W[n]["old"], W[n]["ad"]) for n in RUNS}, "_ad")]
    print("\n".join(L))
    if json_out:
        pathlib.Path(json_out).write_text(json.dumps(data, indent=1, default=lambda o: o.tolist()
                                                     if isinstance(o, np.ndarray) else float(o)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=None)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    main(a.runs, a.json)
