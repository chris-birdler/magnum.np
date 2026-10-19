#!/usr/bin/env python
"""
Damping α 0.01 against α 0.1 at the centre point, β fits (request of Chris 2026-10-13, PROTOKOLL 7.79). Two variants for
α 0.01:
  converged  only stages with a calm tail ≥ 2 cycles (steady rule 7.54): 9 and 20 mT from the 30-cycle runs
             cond_C<k>_a001_b<mT>, 35 … 150 mT from the design runs C<k>_a001;
  all        every seed and stage: mean w_dis of the last cycles (30-cycle runs: last 10, design runs: last 3), also when
             the stage does not settle (a time mean of the loss).
α 0.1: design values (steady rule ≥ 2 with reruns; all 8 seeds converged at all 7 amplitudes).
Fits of ln w over ln B (nominal B): with 9 mT (9 … 150), without 9 mT (20 … 150), and 3-point windows. Per seed where the
seed has all points; else a model with one offset per seed and one common slope (uses all values, gaps allowed).
Δβ = β(α 0.01) − β(α 0.1) from the paired values (seed offsets). Mean ± SE, SD, 90 % CI (t).

    python eval_alpha_fits.py > results/alpha_fits.md       # also results/alpha_fits.png
"""
import json
import math
import pathlib

import numpy as np
from scipy import stats

import analyze
import analyze_design as AD
import eval_alpha_cond as EAC

R = pathlib.Path(__file__).resolve().parent / "runs"
AMPS = list(AD.AMPS_MT)
SEEDS = range(8)
FITS = (("with 9 mT: 9 … 150", tuple(AMPS)), ("without 9 mT: 20 … 150", tuple(AMPS[1:])),
        ("9 / 20 / 35", (9, 20, 35)), ("20 / 35 / 50", (20, 35, 50)), ("35 / 50 / 70", (35, 50, 70)),
        ("50 / 70 / 100", (50, 70, 100)), ("70 / 100 / 150", (70, 100, 150)))


def nominal(bmT):
    return min(AMPS, key=lambda a: abs(a / bmT - 1))


def stage_by_amp(name):
    st = json.loads((R / name / "summary.json").read_text())["stages"]
    return {nominal(1500.0 * s["cycles"][-1]["b_peak"]): s for s in st if s.get("measure") and s["cycles"]}


def last_mean(stage, n):
    return float(np.mean([c["w_dis"] for c in stage["cycles"][-n:]]))


def data():
    L0, Lc, La = {}, {}, {}
    for s in SEEDS:
        ref = AD.run_record(AD.rerun_dir(R / ("C%d" % s)) or R / ("C%d" % s))
        des = AD.run_record(R / ("C%d_a001" % s))
        sb = stage_by_amp("C%d_a001" % s)
        L0[s] = {b: ref["lnw@%gmT" % b] for b in AMPS}
        Lc[s] = {b: des["lnw@%gmT" % b] for b in AMPS}
        La[s] = {b: math.log(last_mean(sb[b], 3)) for b in AMPS if b in sb}
        for b in (9, 20):
            name = "cond_C%d_a001_b%d" % (s, b)
            w = EAC.calm_value(name)[0]
            Lc[s][b] = math.log(w) if w is not None else float("nan")
            st = json.loads((R / name / "summary.json").read_text())["stages"][-1]
            La[s][b] = math.log(last_mean(st, 10))
    return L0, Lc, La


def fe(rows):
    """y = a_seed + slope · ln B; rows (seed, B, y). (slope, SE, n values, n seeds, dof) or None."""
    seeds = sorted({r[0] for r in rows})
    if len(rows) - len(seeds) - 1 < 1:
        return None
    X = np.zeros((len(rows), len(seeds) + 1))
    y = np.array([r[2] for r in rows])
    for i, (s, b, _) in enumerate(rows):
        X[i, seeds.index(s)] = 1.0
        X[i, -1] = math.log(b)
    coef = np.linalg.lstsq(X, y, rcond=None)[0]
    dof = len(rows) - X.shape[1]
    cov = float(((y - X @ coef) ** 2).sum()) / dof * np.linalg.pinv(X.T @ X)
    return coef[-1], math.sqrt(cov[-1, -1]), len(rows), len(seeds), dof


def per_seed(L, amps):
    v = [np.polyfit(np.log(amps), [L[s][b] for b in amps], 1)[0] for s in SEEDS
         if all(np.isfinite(L[s].get(b, np.nan)) for b in amps)]
    return np.array(v)


def fmt_ps(v):
    if len(v) < 3:
        return "n %d: no statement" % len(v)
    h = stats.t.ppf(0.95, len(v) - 1) * v.std(ddof=1) / math.sqrt(len(v))
    return "%.2f / %.2f / %.2f … %.2f (n %d)" % (v.mean(), v.std(ddof=1), v.mean() - h, v.mean() + h, len(v))


def fmt_fe(r):
    if r is None:
        return "too few values"
    h = stats.t.ppf(0.95, r[4]) * r[1]
    return "%+.2f ± %.2f (90 %% CI %+.2f … %+.2f; %d values, %d seeds)" % (r[0], r[1], r[0] - h, r[0] + h, r[2], r[3])


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    L0, Lc, La = data()
    out = ["# Damping α 0.01 against α 0.1: β fits (PROTOKOLL 7.79)", "",
           "Centre point, seeds C0 … C7, P ∝ w (loss per cycle). Variant **converged**: only stages with ≥ 2 calm "
           "cycles. Variant **all**: every stage, time mean of the dissipation of the last cycles (10 at 9 and 20 mT, "
           "3 at 35 … 150 mT). α 0.1: all 8 seeds converged at all B. Per-seed fits: mean / SD / 90 % CI; Δβ from the "
           "paired values with one offset per seed (± SE, 90 % CI).", ""]
    out += ["## Converged values at α 0.01", "", "| seed | " + " | ".join("%g" % b for b in AMPS) + " |",
            "|---|" + "---|" * len(AMPS)]
    for s in SEEDS:
        out.append("| C%d | %s |" % (s, " | ".join("✔" if np.isfinite(Lc[s][b]) else "–" for b in AMPS)))
    out += ["", "## Loss ratio P(α 0.01) / P(α 0.1) per B (geometric mean, 90 % CI, n)", "",
            "| B (mT) | converged | all |", "|---|---|---|"]
    for b in AMPS:
        cells = []
        for L in (Lc, La):
            r = np.array([L[s][b] - L0[s][b] for s in SEEDS if np.isfinite(L[s].get(b, np.nan))])
            if len(r) >= 3:
                h = stats.t.ppf(0.95, len(r) - 1) * r.std(ddof=1) / math.sqrt(len(r))
                cells.append("× %.2f (%.2f … %.2f, n %d)" % (math.exp(r.mean()), math.exp(r.mean() - h),
                                                             math.exp(r.mean() + h), len(r)))
            else:
                cells.append("n %d: no statement" % len(r))
        out.append("| %g | %s | %s |" % (b, *cells))
    for vname, L in (("converged", Lc), ("all", La)):
        out += ["", "## β, variant %s" % vname, "",
                "| fit (mT) | β α 0.1, per seed | β α 0.01, per seed | β α 0.01, seed offsets | Δβ (α 0.01 − α 0.1), seed offsets |",
                "|---|---|---|---|---|"]
        for lab, amps in FITS:
            r1 = [(s, b, L[s][b]) for s in SEEDS for b in amps if np.isfinite(L[s].get(b, np.nan))]
            rd = [(s, b, L[s][b] - L0[s][b]) for s in SEEDS for b in amps if np.isfinite(L[s].get(b, np.nan))]
            out.append("| %s | %s | %s | %s | %s |" % (lab, fmt_ps(per_seed(L0, amps)), fmt_ps(per_seed(L, amps)),
                                                       fmt_fe(fe(r1)), fmt_fe(fe(rd))))
    print("\n".join(out))
    plot(L0, Lc, La)


def plot(L0, Lc, La):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(15.5, 4.6))
    for ax, L, title in ((axs[0], Lc, "α 0.01 converged stages only"), (axs[1], La, "α 0.01 all stages (time mean)")):
        for s in SEEDS:
            ax.loglog(AMPS, [math.exp(L0[s][b]) for b in AMPS], "-", color="#9B111E", lw=0.8, alpha=0.5)
            y = [math.exp(L[s][b]) if np.isfinite(L[s].get(b, np.nan)) else np.nan for b in AMPS]
            ax.loglog(AMPS, y, "o-", color="#1F4E79", lw=0.8, ms=3, alpha=0.7)
        ax.plot([], [], "-", color="#9B111E", label="α 0.1, one line per seed")
        ax.plot([], [], "o-", color="#1F4E79", label="α 0.01, one line per seed")
        ax.set_xticks(AMPS)
        ax.set_xticklabels([str(b) for b in AMPS])
        ax.minorticks_off()
        ax.set_xlabel("peak flux density B (mT)")
        ax.set_ylabel("loss per cycle w (normalised)")
        ax.set_title(title, fontsize=11)
        ax.legend(fontsize=9)
    ax = axs[2]
    wins = [f for f in FITS if len(f[1]) == 3]
    xc = [math.exp(np.mean(np.log(a))) for _, a in wins]
    for L, col, mk, lab, off in ((L0, "#9B111E", "o", "α 0.1", 0.97), (Lc, "#1F4E79", "s", "α 0.01 converged", 1.0),
                                 (La, "#666666", "^", "α 0.01 all", 1.03)):
        m, e = [], []
        for _, amps in wins:
            r = fe([(s, b, L[s][b]) for s in SEEDS for b in amps if np.isfinite(L[s].get(b, np.nan))])
            m.append(r[0] if r else np.nan)
            e.append(2 * r[1] if r else np.nan)
        ax.errorbar([x * off for x in xc], m, yerr=e, fmt=mk + "-", color=col, capsize=3, label=lab)
    ax.axhline(2.0, color="black", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_xticks(AMPS)
    ax.set_xticklabels([str(b) for b in AMPS])
    ax.minorticks_off()
    ax.set_xlabel("peak flux density B (mT), centre of the 3-point window")
    ax.set_ylabel("β of the window (seed offsets, ± 2 SE)")
    ax.set_title("β in 3-point windows", fontsize=11)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(pathlib.Path(__file__).resolve().parent / "results" / "alpha_fits.png", dpi=130)


if __name__ == "__main__":
    main()
