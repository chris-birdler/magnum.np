#!/usr/bin/env python
"""
Mesh error at d/l_ex 212 and 300 for the base microstructure (PLAN 4.6, option 2; rules fixed before the result).
Base microstructure: L_eff 12 l_ex, r_p 0, φ 0.65, α 0.1, converged virgin relaxation; M212_* (seeds 931 … 934) and
M300_* (mesh300, seeds 921 … 924); dx 3 against dx 1.5 with the same seed (paired); 9 / 50 / 100 mT; steady-tail rule
with ≥ 2 calm cycles; finished reruns X_mc2 replace X. K3 (d/l_ex 212, other microstructure) is a reference row.
The old pairs with a not converged relaxation (d/l_ex 96: P3_dx*/R_dx15, d/l_ex 212: R212_*; 9 and 50 mT) are read
in the audit mode and are shown as hints only (PROTOKOLL 7.64).

Mesh error per seed: Δ ln P = ln[P(dx 3) / P(dx 1.5)] at each amplitude; Δβ = β(dx 3) − β(dx 1.5), with β the slope
of the power-law fit of ln P over ln B (9 / 50 / 100 mT, and 9 → 50 mT for the comparison with the old pairs).
Per group: mean ± SE over the seeds and the 90 % CI (t, n − 1); fewer than 3 seeds: no statement.
Mesh share of the size lever (d/l_ex 212 → 300): Δ(300) − Δ(212) with SE √(SE₁² + SE₂²) and the 90 % CI
(Welch-Satterthwaite degrees of freedom). The size lever on each mesh is shown next to it.
Sign change: Δ(212) and Δ(300) have opposite signs and both 90 % CIs exclude 0.

    python eval_size.py > results/mesh212.md
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
MIN_N = 3
AMPS = (9, 50, 100)
GROUPS = (("M212", "d/l_ex 212 (0.71 µm)", range(931, 935)), ("M300", "d/l_ex 300 (1.0 µm)", range(921, 925)),
          ("MK3", "K3: d/l_ex 212, φ 0.69, L_eff 18, r_p 1.5", range(2031, 2037)))
# old pairs (relaxation not converged; one run per seed with the stages 9 → 50 mT): (label, fine, coarse, seeds)
OLD = (("d/l_ex 96, dx a = 3", "R_dx15_s%d", "P3_dx3_s%d", (3, 4, 5)),
       ("d/l_ex 96, dx a = 2", "R_dx15_s%d", "P3_dx2_s%d", (3, 4, 5)),
       ("d/l_ex 96, dx a = 1", "R_dx15_s%d", "P3_dx1_s%d", (3, 4, 5)),
       ("d/l_ex 212, dx a = 3", "R212_dx15_s%d", "R212_dx3_s%d", (911, 912, 913, 914)))
OUTPUTS = ("lnP 9", "lnP 50", "lnP 100", "beta", "beta 9-50")


def run_dir(name):
    d = R / name
    return AD.rerun_dir(d) or d


def nominal(B):
    return 9 if B < 15 else (50 if B < 75 else 100)


def stages(d):
    """{nominal mT: (w, B in mT)} of the steady stages of one run folder."""
    if not (d / "DONE").exists():
        return {}
    _, rows = analyze.load_run(d)
    return {nominal(1500.0 * r["b_peak"]): (r["w_loop"], 1500.0 * r["b_peak"]) for r in rows
            if np.isfinite(r["w_loop"]) and r["w_loop"] > 0}


def load_group(prefix, seeds):
    """{(mesh, seed): {mT: (w, B)}} for the jobs prefix_mesh_s<seed>_b<mT>."""
    out = {}
    for mesh in ("dx15", "dx3"):
        for s in seeds:
            st = {}
            for b in AMPS:
                st.update({k: v for k, v in stages(run_dir("%s_%s_s%d_b%d" % (prefix, mesh, s, b))).items() if k == b})
            out[(mesh, s)] = st
    return out


def beta(st, amps):
    if not all(b in st for b in amps):
        return None
    return np.polyfit(np.log([st[b][1] for b in amps]), np.log([st[b][0] for b in amps]), 1)[0]


def mesh_error(fine, coarse):
    """{output: value} of one seed from the stage dicts of dx 1.5 (fine) and the coarse mesh."""
    out = {}
    for b in AMPS:
        if b in fine and b in coarse:
            out["lnP %d" % b] = math.log(coarse[b][0] / fine[b][0])
    for key, amps in (("beta", AMPS), ("beta 9-50", (9, 50))):
        bc, bf = beta(coarse, amps), beta(fine, amps)
        if bc is not None and bf is not None:
            out[key] = bc - bf
    return out


def summary(v):
    """(mean, SE, n, half width of the 90 % CI); None for fewer than MIN_N values."""
    v = np.asarray(list(v), float)
    if len(v) < MIN_N:
        return None
    se = v.std(ddof=1) / math.sqrt(len(v))
    return v.mean(), se, len(v), stats.t.ppf(0.95, len(v) - 1) * se


def difference(a, b):
    """a − b of two summaries: (mean, SE, half width of the 90 % CI with Welch-Satterthwaite dof) or None."""
    if a is None or b is None:
        return None
    se = math.sqrt(a[1] ** 2 + b[1] ** 2)
    dof = se ** 4 / (a[1] ** 4 / (a[2] - 1) + b[1] ** 4 / (b[2] - 1)) if se > 0 else 1.0
    return a[0] - b[0], se, stats.t.ppf(0.95, dof) * se


def sign_change(a, b):
    return (a is not None and b is not None and a[0] * b[0] < 0 and abs(a[0]) > a[3] and abs(b[0]) > b[3])


def fmt(sm, factor=False, n=None):
    if sm is None:
        return "n %d: no statement" % n if n is not None else "–"
    m, se, k, h = sm
    if factor:
        return "× %.2f (ln %+.2f ± %.2f, CI %+.2f … %+.2f, n %d)" % (math.exp(m), m, se, m - h, m + h, k)
    return "%+.2f ± %.2f (CI %+.2f … %+.2f, n %d)" % (m, se, m - h, m + h, k)


def per_seed(group):
    seeds = sorted({s for _, s in group})
    vals = {k: [] for k in OUTPUTS}
    for s in seeds:
        for k, v in mesh_error(group[("dx15", s)], group[("dx3", s)]).items():
            vals[k].append(v)
    return vals


def table(rows, coarse="dx 3"):
    print("| group | P(%s)/P(dx 1.5) 9 mT | 50 mT | 100 mT | Δβ (9 / 50 / 100 mT) | Δβ (9 → 50 mT) |" % coarse)
    print("|---|---|---|---|---|---|")
    for lab, vals in rows:
        cells = [fmt(summary(vals[k]), k.startswith("lnP"), len(vals[k])) if vals[k] else "–" for k in OUTPUTS]
        print("| %s | %s |" % (lab, " | ".join(cells)))


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    data = {p: load_group(p, seeds) for p, _, seeds in GROUPS}
    kd = None
    for p, _, seeds in GROUPS:
        for s in seeds:
            c = run_dir("%s_dx15_s%d_b9" % (p, s)) / "config.json"
            if kd is None and c.exists():
                kd = json.loads(c.read_text())["units"]["Kd"]
    to_P = kd * 30e6 / 1e6
    print("# Mesh error at d/l_ex 212 and 300 (base microstructure), dx 3 against dx 1.5\n")
    print("Base microstructure: L_eff 12 l_ex, r_p 0, φ 0.65, α 0.1, A 10 pJ/m (l_ex 3.34 nm; the µm values hold for "
          "this A). Steady loops (≥ 2 calm cycles), reruns replace their original. dx in units of l_ex. Errors: SE "
          "over the seeds, 90 %% CI with t(n − 1); fewer than %d seeds: no statement.\n" % MIN_N)
    print("## Loss density P (W/cm³, geometric mean over the seeds)\n")
    print("| group | mesh | 9 mT | 50 mT | 100 mT |")
    print("|---|---|---|---|---|")
    for p, lab, seeds in GROUPS:
        for mesh, mlab in (("dx15", "dx 1.5"), ("dx3", "dx 3")):
            cells = []
            for b in AMPS:
                v = [math.log(data[p][(mesh, s)][b][0] * to_P) for s in seeds if b in data[p][(mesh, s)]]
                cells.append("%.3g (n %d)" % (math.exp(np.mean(v)), len(v)) if v else "–")
            print("| %s | %s | %s |" % (lab, mlab, " | ".join(cells)))
    vals = {p: per_seed(data[p]) for p, _, _ in GROUPS}
    print("\n## Mesh error (paired per seed)\n")
    table([(lab, vals[p]) for p, lab, _ in GROUPS])
    print("\n## Mesh share of the size lever (d/l_ex 212 → 300, base microstructure)\n")
    print("Size lever on one mesh: mean ln P (or β) at 300 minus at 212 (different seeds, Welch). Mesh share = lever "
          "at dx 3 − lever at dx 1.5 = Δ(300) − Δ(212). Sign change: opposite signs, both 90 % CIs exclude 0.\n")
    print("| output | lever at dx 1.5 | lever at dx 3 | mesh share Δ(300) − Δ(212) | sign change |")
    print("|---|---|---|---|---|")
    for k in OUTPUTS:
        lev = {}
        for mesh in ("dx15", "dx3"):
            v = {}
            for p in ("M212", "M300"):
                g = data[p]
                if k.startswith("lnP"):
                    b = int(k.split()[1])
                    v[p] = [math.log(g[(mesh, s)][b][0]) for (m, s) in g if m == mesh and b in g[(mesh, s)]]
                else:
                    amps = AMPS if k == "beta" else (9, 50)
                    v[p] = [x for x in (beta(g[(mesh, s)], amps) for (m, s) in g if m == mesh) if x is not None]
            lev[mesh] = difference(summary(v["M300"]), summary(v["M212"]))
        a, b = summary(vals["M300"][k]), summary(vals["M212"][k])
        sh = difference(a, b)
        cell = lambda d: "%+.2f ± %.2f (CI %+.2f … %+.2f)" % (d[0], d[1], d[0] - d[2], d[0] + d[2]) if d else "–"
        print("| %s | %s | %s | %s | %s |" % (k, cell(lev["dx15"]), cell(lev["dx3"]), cell(sh),
                                              "yes" if sign_change(b, a) else "no"))
    print("\n## Old pairs (relaxation not converged; hints only, PROTOKOLL 7.64)\n")
    analyze.ALLOW_UNCONVERGED = True              # audit mode for the old runs only
    rows = []
    for lab, f_fmt, c_fmt, seeds in OLD:
        v = {k: [] for k in OUTPUTS}
        for s in seeds:
            for k, x in mesh_error(stages(R / (f_fmt % s)), stages(R / (c_fmt % s))).items():
                v[k].append(x)
        rows.append((lab, v))
    analyze.ALLOW_UNCONVERGED = False
    table(rows, "dx a")


if __name__ == "__main__":
    main()
