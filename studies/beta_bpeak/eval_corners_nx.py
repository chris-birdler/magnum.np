#!/usr/bin/env python
"""
Corner test (decision Chris 2026-10-17: first the corners, then decide on a new main study; rules fixed before the
runs, commit a3363b4). The 12 corner runs of the main study (K1 … K4, 3 seeds each, dx 3 l_ex, 9 … 150 mT) again with
the production setup of PLAN 4.13 (conv_jobs.py --newell corners; runs <run>_nx): final stray-field setting
(newell_ext, 4 shells), tight cycle tolerance, relaxation with atol 1e-6 and the absolute stop threshold, online stop,
second drive correction. Δβ measures the setup change as a whole; the test cannot split it (PLAN 4.13 D3).
Value per run: β = slope of ln w over ln B at the nominal amplitudes (main-fit definition; analyze_design.run_record
with the steady-tail rule ≥ 2 calm cycles). The old partner of a new run is its _mc2 rerun where it exists (K1_0,
K4_0, K4_1), as in the old main result (PLAN 4.5 C). Paired per run: both fits use the amplitudes with a value in the
old and the new run, at least 6 of 7; Δβ = β(new) − β(old).
Reading:
  "β robust to the setup change"  the 90 % CI of the mean Δβ over the 12 runs lies inside ± 0.10, and no corner has
                                  |mean Δβ| > 0.20 (n 3 per corner);
  "β shifted"                     the 90 % CI of the mean Δβ excludes 0 and |mean Δβ| > 0.10;
  "not determined"                all other cases, or fewer than 9 runs with β in both.
Information: per corner β old and new (mean ± SE, n 3) and the corner differences β(Ki) − β(K1) old and new;
per amplitude the mean Δ ln w (new − old) with SE (n 12); PLAN 4.13 S9 flags of the new stages (mean |ε| of the calm
tail > 0.02, SD of ln w_loop in the calm tail > 0.05) and the mean Δβ without the flagged stages; the energy-terms check
of S5 (d = ε + ΔE_int / w_loop per calm cycle, the numerical error of the energy balance); relaxation iterations.

    python eval_corners_nx.py > runs/corners_nx_report.md
"""
import json
import math
import pathlib

import numpy as np
from scipy import stats

import analyze
import analyze_design as AD
import conv_jobs as CJ

R = pathlib.Path(__file__).resolve().parent / "runs"
LIMIT, CORNER_LIMIT, MIN_RUNS = 0.10, 0.20, 9


def lnw(d):
    r = AD.run_record(d)
    return np.array([r["lnw@%gmT" % b] for b in AD.AMPS_MT]) if r else np.full(len(AD.AMPS_MT), np.nan)


MIN_AMPS = 6
FLAG_EPS, FLAG_W = 0.02, 0.05          # PLAN 4.13 S9


def old_dir(root, n):
    """The old partner of corner run n: its _mc2 rerun if finished (PLAN 4.5 C), else the original run."""
    return AD.rerun_dir(root / n) or root / n


def measured_stages(d):
    """The measured stages of a run folder (summary.json), in the order of the amplitudes; [] without DONE."""
    if not (d / "DONE").exists():
        return []
    return [s for s in json.loads((d / "summary.json").read_text())["stages"] if s.get("measure")]


def calm_tail(st):
    return analyze.steady_tail(st["cycles"][analyze.first_measured(st):], AD.BAL_EVENT)


def stage_flags(st):
    """PLAN 4.13 S9: (flag ε, flag w) of one stage on its calm tail; (False, False) for a tail < 2 cycles."""
    tail = calm_tail(st)
    if len(tail) < 2:
        return False, False
    eps = np.array([c["w_dis"] / c["w_loop"] - 1.0 for c in tail])
    lw = np.log([c["w_loop"] for c in tail])
    return bool(np.abs(eps).mean() > FLAG_EPS), bool(lw.std(ddof=1) > FLAG_W)


def balance(st):
    """PLAN 4.13 S5 energy-terms check of one stage: per cycle (cycle, ε, −ΔE_int / w_loop, d = ε + ΔE_int / w_loop).
    ΔE_int: change over the cycle of the sum of all terms without the drive (E_start of the next cycle, E_end for the
    last). [] if the run has no energy terms."""
    cyc, out = st["cycles"], []
    for i, c in enumerate(cyc):
        if "E_start" not in c:
            return []
        e1 = cyc[i + 1]["E_start"] if i + 1 < len(cyc) else st.get("E_end")
        if e1 is None or c["w_loop"] == 0:
            continue
        dE = sum(v for k, v in e1.items() if k != "drive") - sum(v for k, v in c["E_start"].items() if k != "drive")
        eps = c["w_dis"] / c["w_loop"] - 1.0
        out.append((c["cycle"], eps, -dE / c["w_loop"], eps + dE / c["w_loop"]))
    return out


def beta_pair(o, n):
    """(β old, β new) over the amplitudes with a value in both runs (at least MIN_AMPS), else NaN."""
    m = np.isfinite(o) & np.isfinite(n)
    if m.sum() < MIN_AMPS:
        return float("nan"), float("nan")
    x = np.log(np.asarray(AD.AMPS_MT))[m]
    return float(np.polyfit(x, o[m], 1)[0]), float(np.polyfit(x, n[m], 1)[0])


def stat(v):
    v = np.array([x for x in v if np.isfinite(x)])
    if len(v) < 2:
        return None
    se = v.std(ddof=1) / math.sqrt(len(v))
    h = stats.t.ppf(0.95, len(v) - 1) * se
    return v.mean(), se, v.mean() - h, v.mean() + h, len(v)


def reading(db, by_corner):
    s = stat(db)
    if s is None or s[4] < MIN_RUNS:
        return "not determined", s
    if LIMIT > s[3] and s[2] > -LIMIT and all(abs(np.nanmean(v)) <= CORNER_LIMIT for v in by_corner.values()):
        return "β robust to the setup change", s
    if (s[2] > 0 or s[3] < 0) and abs(s[0]) > LIMIT:
        return "β shifted", s
    return "not determined", s


def main(root=R):
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    runs = [n.replace(CJ.NX_SUFFIX, "") for n, _ in CJ.corner_newell_lines()]
    W = {n: (lnw(old_dir(root, n)), lnw(root / (n + CJ.NX_SUFFIX))) for n in runs}
    B = {n: beta_pair(o, nw) for n, (o, nw) in W.items()}
    db = [B[n][1] - B[n][0] for n in runs]
    by_c = {c: [B[n][1] - B[n][0] for n in runs if n.startswith(c + "_")] for c in CJ.CORNERS_NX}
    key, s = reading(db, by_c)
    L = ["# Corner test: the 12 corner runs with the production setup of PLAN 4.13 (eval_corners_nx.py)", "",
         "β = slope of ln w over ln B, 9 … 150 mT (main-fit definition). Old: DemagFieldPBC, atol 1e-5 (the _mc2 rerun "
         "where it exists); new: PLAN 4.13. Δβ = new − old, paired per run.", "",
         "| run (old partner) | β old | β new | Δβ |", "|---|---|---|---|"]
    for n in runs:
        L.append("| %s (%s) | %.3f | %.3f | %+.3f |" % (n, old_dir(root, n).name, B[n][0], B[n][1], B[n][1] - B[n][0]))
    L += ["", "## Per corner (mean ± SE, n 3)", "", "| corner | β old | β new | Δβ |", "|---|---|---|---|"]
    co = {}
    for c in CJ.CORNERS_NX:
        bo = stat([B[n][0] for n in runs if n.startswith(c + "_")])
        bn = stat([B[n][1] for n in runs if n.startswith(c + "_")])
        d = stat(by_c[c])
        co[c] = (bo, bn)
        f = lambda x: "%.3f ± %.3f" % (x[0], x[1]) if x else "–"
        L.append("| %s | %s | %s | %s |" % (c, f(bo), f(bn), "%+.3f ± %.3f" % (d[0], d[1]) if d else "–"))
    L += ["", "Corner differences β(Ki) − β(K1), old / new: " + "; ".join(
        "%s %+.3f / %+.3f" % (c, co[c][0][0] - co["K1"][0][0], co[c][1][0] - co["K1"][1][0])
        for c in CJ.CORNERS_NX[1:] if co[c][0] and co[c][1] and co["K1"][0] and co["K1"][1]) + "."]
    L += ["", "Mean Δ ln w per amplitude (new − old, SE, n): " + "; ".join(
        "%g mT %s" % (b, "%+.3f ± %.3f (%d)" % (x[0], x[1], x[4]) if (x := stat([W[n][1][k] - W[n][0][k] for n in runs]))
                       else "–") for k, b in enumerate(AD.AMPS_MT)) + "."]
    L += ["", "Mean Δβ over the runs: %s. **Reading: %s.**" % (
        "%+.3f ± %.3f (SE, n %d; 90 %% CI %+.3f … %+.3f)" % (s[0], s[1], s[4], s[2], s[3]) if s else "–", key)]
    L += info_lines(root, runs, W)
    print("\n".join(L))


def info_lines(root, runs, W):
    """PLAN 4.13 information parts: S9 flags, Δβ without the flagged stages, S5 energy-terms check, relaxation."""
    L = ["", "## Information (PLAN 4.13)", ""]
    fl, Wf, bal = {}, {}, {k: [] for k in range(len(AD.AMPS_MT))}
    for n in runs:
        d = root / (n + CJ.NX_SUFFIX)
        sts = measured_stages(d)
        f = [stage_flags(s) for s in sts]
        fl[n] = f
        w = W[n][1].copy()
        for k, (fe, fw) in enumerate(f):
            if (fe or fw) and k < len(w):
                w[k] = np.nan
        Wf[n] = w
        for k, s in enumerate(sts):
            tail = {c["cycle"] for c in calm_tail(s)}
            bal[k] += [r for r in balance(s) if r[0] in tail]
    ne = sum(fe for f in fl.values() for fe, _ in f)
    nw = sum(fw for f in fl.values() for _, fw in f)
    ns = sum(len(f) for f in fl.values())
    db = [beta_pair(W[n][0], Wf[n])[1] - beta_pair(W[n][0], Wf[n])[0] for n in runs]
    s = stat(db)
    L.append("- Flags (S9) of the new stages: flag ε %d, flag w %d of %d stages. Mean Δβ without the flagged stages: %s." % (
        ne, nw, ns, "%+.3f ± %.3f (SE, n %d)" % (s[0], s[1], s[4]) if s else "–"))
    L += ["", "Energy-terms check (S5) on the calm-tail cycles of the new runs: d = ε + ΔE_int / w_loop.", "",
          "| amplitude | cycles | median \\|ε\\| | median \\|d\\| | max \\|d\\| |", "|---|---|---|---|---|"]
    for k, b in enumerate(AD.AMPS_MT):
        r = np.array([x[1:] for x in bal[k]]) if bal[k] else None
        L.append("| %g mT | %d | %s | %s | %s |" % (b, len(bal[k]), *(("%.4f" % np.median(np.abs(r[:, 0])),
                 "%.4f" % np.median(np.abs(r[:, 2])), "%.4f" % np.abs(r[:, 2]).max()) if r is not None else ("–",) * 3)))
    it = []
    for n in runs:
        d = root / (n + CJ.NX_SUFFIX) / "summary.json"
        if d.exists():
            i = json.loads(d.read_text()).get("init", {})
            it.append("%s %s" % (n, i.get("relax_iterations", "–")))
    L += ["", "Relaxation iterations of the new runs: " + "; ".join(it) + "."]
    return L


if __name__ == "__main__":
    main()
