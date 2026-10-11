#!/usr/bin/env python
"""
Size trend of the mesh error (PLAN 4.6, rules fixed before the result): base microstructure (L_eff 12 l_ex, r_p 0,
φ 0.65) at d/l_ex 150 / 212 / 300 (M150_*, M212_*, M300_*), dx 3 against dx 1.5 with the same seed (paired),
9 / 50 / 100 mT, steady-tail rule with ≥ 2 calm cycles, finished reruns X_mc2 replace X. K3 (d/l_ex 212, other
microstructure) is shown as a reference row; it is not part of the trend.

Per size: Δ ln P = ln[P(dx 3) / P(dx 1.5)] per seed and amplitude, Δβ = β(dx 3) − β(dx 1.5) per seed (power-law fit
over 9 / 50 / 100 mT); mean ± SE over the seeds and the 90 % CI (t, n − 1).
Trend label per output (150 → 212 → 300): "falls" if |mean| decreases at both steps by more than the SE of the step
(√(SE₁² + SE₂²)), "rises" if it increases at both steps by more than that, else "no clear trend".
Estimates for larger particles (not a measurement): least-squares line of the per-seed values over x = (l_ex/d)^k,
k = 1 (surface share) and k = 2 (vortex-core share), evaluated at d/l_ex 600 and d → ∞, with the SE of the fit.

    python eval_size.py > results/size_trend.md
"""
import json
import math
import pathlib

import numpy as np
from scipy import stats

import analyze
import analyze_design as AD
import eval_k3fill as KF

R = pathlib.Path(__file__).resolve().parent / "runs"
TREND = ((150, "M150", range(941, 945)), (212, "M212", range(931, 935)), (300, "M300", range(921, 925)))
REFERENCE = ("K3 (212; φ 0.69, L_eff 18, r_p 1.5, A 20 pJ/m)", "MK3", range(2031, 2037))
MESHES = ("dx15", "dx3")
AMPS = (9, 50, 100)
EXTRAP = (600.0, math.inf)


def run_dir(name):
    d = R / name
    return AD.rerun_dir(d) or d


def load(prefix, seeds):
    """{(mesh, seed, mT): (w, B in mT) or None}, d/l_ex used and gap in cells at dx 3 (None if no run exists)."""
    L, d_lex, gap = {}, None, None
    for mesh in MESHES:
        for s in seeds:
            for b in AMPS:
                d = run_dir("%s_%s_s%d_b%d" % (prefix, mesh, s, b))
                L[(mesh, s, b)] = KF.loss(d)
                if (d / "config.json").exists():
                    c = json.loads((d / "config.json").read_text())
                    d_lex = c["d_lex_used"]
                    if mesh == "dx3":
                        gap = c["gap_cells"]
    return L, d_lex, gap


def beta(L, mesh, s):
    pts = [L[(mesh, s, b)] for b in AMPS]
    if not all(pts):
        return None
    return np.polyfit(np.log([p[1] for p in pts]), np.log([p[0] for p in pts]), 1)[0]


def paired(L, seeds):
    """{output: {seed: value}} for the outputs 'lnP 9', 'lnP 50', 'lnP 100' and 'beta'."""
    out = {}
    for b in AMPS:
        out["lnP %d" % b] = {s: math.log(L[("dx3", s, b)][0] / L[("dx15", s, b)][0]) for s in L_seeds(L)
                             if L[("dx3", s, b)] and L[("dx15", s, b)]}
    bd = {}
    for s in L_seeds(L):
        b3, b15 = beta(L, "dx3", s), beta(L, "dx15", s)
        if b3 is not None and b15 is not None:
            bd[s] = b3 - b15
    out["beta"] = bd
    return out


def L_seeds(L):
    return sorted({k[1] for k in L})


def summary(v):
    """(mean, SE, n, half width of the 90 % CI) or None for fewer than 2 values."""
    v = np.asarray(list(v), float)
    if len(v) < 2:
        return None
    se = v.std(ddof=1) / math.sqrt(len(v))
    return v.mean(), se, len(v), stats.t.ppf(0.95, len(v) - 1) * se


def trend_label(sums):
    """'falls' / 'rises' / 'no clear trend' from the summaries in the order 150 → 212 → 300."""
    if any(x is None for x in sums):
        return "not determined (missing size)"
    steps = [abs(b[0]) - abs(a[0]) for a, b in zip(sums, sums[1:])]
    ses = [math.sqrt(a[1] ** 2 + b[1] ** 2) for a, b in zip(sums, sums[1:])]
    if all(st < -se for st, se in zip(steps, ses)):
        return "falls"
    if all(st > se for st, se in zip(steps, ses)):
        return "rises"
    return "no clear trend"


def extrapolate(points, k):
    """points: [(d/l_ex, value)]; line over x = (1/d)^k; returns {d: (value, SE)} for EXTRAP."""
    x = np.array([(1.0 / d) ** k for d, _ in points])
    y = np.array([v for _, v in points])
    if len(y) < 4 or len(set(x)) < 2:
        return None
    X = np.column_stack([np.ones_like(x), x])
    coef = np.linalg.lstsq(X, y, rcond=None)[0]
    dof = len(y) - 2
    s2 = float(((y - X @ coef) ** 2).sum()) / dof
    cov = s2 * np.linalg.inv(X.T @ X)
    out = {}
    for d in EXTRAP:
        g = np.array([1.0, 0.0 if math.isinf(d) else (1.0 / d) ** k])
        out[d] = (float(g @ coef), math.sqrt(float(g @ cov @ g)))
    return out


def fmt(sm, factor=False):
    if sm is None:
        return "–"
    m, se, n, h = sm
    if factor:
        return "× %.2f (ln %+.2f ± %.2f, CI %+.2f … %+.2f, n %d)" % (math.exp(m), m, se, m - h, m + h, n)
    return "%+.2f ± %.2f (CI %+.2f … %+.2f, n %d)" % (m, se, m - h, m + h, n)


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    to_P = None
    rows = []
    for d_nom, prefix, seeds in TREND + ((212, REFERENCE[1], REFERENCE[2]),):
        L, d_lex, gap = load(prefix, seeds)
        if to_P is None:
            for s in seeds:
                c = run_dir("%s_dx15_s%d_b9" % (prefix, s)) / "config.json"
                if c.exists():
                    to_P = json.loads(c.read_text())["units"]["Kd"] * 30e6 / 1e6
                    break
        rows.append(dict(prefix=prefix, d=d_lex or float(d_nom), gap=gap, L=L, out=paired(L, seeds)))
    trend, ref = rows[:len(TREND)], rows[len(TREND)]
    label = {r["prefix"]: "%s (d/l_ex %.0f)" % (r["prefix"], r["d"]) for r in trend}
    label[ref["prefix"]] = REFERENCE[0]

    print("# Size trend of the mesh error (dx 3 against dx 1.5)\n")
    print("Base microstructure (L_eff 12 l_ex, r_p 0, φ 0.65, A 10 pJ/m: d/l_ex 150 / 212 / 300 = 0.49 / 0.71 / "
          "1.0 µm). Steady loops (≥ 2 calm cycles), reruns replace their original. Errors: SE over the seeds, "
          "90 % CI with t(n − 1). n = seeds with a value on both meshes.\n")
    print("## Loss density P (W/cm³, geometric mean over the seeds)\n")
    print("| run group | gap at dx 3 (cells) | mesh | 9 mT | 50 mT | 100 mT |")
    print("|---|---|---|---|---|---|")
    for r in trend + [ref]:
        for mesh, mlab in zip(MESHES, ("dx 1.5", "dx 3")):
            cells = []
            for b in AMPS:
                v = [math.log(r["L"][(mesh, s, b)][0] * to_P) for s in L_seeds(r["L"]) if r["L"][(mesh, s, b)]]
                cells.append("%.3g (n %d)" % (math.exp(np.mean(v)), len(v)) if v else "–")
            print("| %s | %s | %s | %s |" % (label[r["prefix"]], "%.1f" % r["gap"] if r["gap"] else "–", mlab,
                                             " | ".join(cells)))
    print("\n## Mesh error per size (paired per seed)\n")
    print("| run group | P(dx 3)/P(dx 1.5) 9 mT | 50 mT | 100 mT | Δβ (9 / 50 / 100 mT) |")
    print("|---|---|---|---|---|")
    for r in trend + [ref]:
        o = r["out"]
        print("| %s | %s | %s | %s | %s |" % (label[r["prefix"]], fmt(summary(o["lnP 9"].values()), True),
                                              fmt(summary(o["lnP 50"].values()), True),
                                              fmt(summary(o["lnP 100"].values()), True),
                                              fmt(summary(o["beta"].values()))))
    print("\n## Trend over the size (150 → 212 → 300) and estimates for larger particles\n")
    print("Trend label: |mean| changes at both steps by more than the SE of the step. Estimates: line over "
          "(l_ex/d)^k through the per-seed values of the 3 sizes (estimate, not a measurement; valid only if no "
          "new loss mechanism appears above 1 µm).\n")
    print("| output | trend | k = 1: d/l_ex 600 | k = 1: d → ∞ | k = 2: d/l_ex 600 | k = 2: d → ∞ |")
    print("|---|---|---|---|---|---|")
    for key, name in (("lnP 9", "ln P ratio 9 mT"), ("lnP 50", "ln P ratio 50 mT"), ("lnP 100", "ln P ratio 100 mT"),
                      ("beta", "Δβ")):
        sums = [summary(r["out"][key].values()) for r in trend]
        pts = [(r["d"], v) for r in trend for v in r["out"][key].values()]
        cells = []
        for k in (1, 2):
            e = extrapolate(pts, k) if all(sums) else None
            cells += ["%+.2f ± %.2f" % e[d] if e else "–" for d in EXTRAP]
        print("| %s | %s | %s |" % (name, trend_label(sums), " | ".join(cells)))


if __name__ == "__main__":
    main()
