"""Task 2: current values of every number in the deck slide_decks/beta_powder_sobol (audit 2026-10-12).
Recomputed with the audit code (main rule: steady tail >= 2 calm cycles, reruns), compared with the numbers in the
built deck (beta_powder_sobol.pptx, 11 Oct 2026) and with results/deck_data.json (11 Oct 16:42).
Output: deck_check.md
"""
import json
import math

import numpy as np
from scipy import stats

import audit_lib as al
import design_data as dd
import sobol_audit as sa


def welch(a, b):
    a, b = np.asarray(a), np.asarray(b)
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    return a.mean() - b.mean(), math.sqrt(va + vb), float(stats.ttest_ind(a, b, equal_var=False).pvalue)


def mesh_P(fmt, seeds, b, dx):
    out = []
    for s in seeds:
        n = al.use_name(fmt % (dx, s, b))
        _, _, st = al.stages(n, rule="steady", min_calm=2)
        if not st[0]["missing"]:
            out.append(st[0]["P"])
    return np.array(out)


def fit_seed(fmt, seeds, amps, dx):
    out = {}
    for s in seeds:
        v = [al.stages(al.use_name(fmt % (dx, s, b)), rule="steady", min_calm=2)[2][0] for b in amps]
        if all(not x["missing"] for x in v):
            out[s] = al.powerfit(np.log([x["B_mT"] for x in v]), np.log([x["P"] for x in v]))[0]
    return out


def main():
    D = json.loads((al.STUDY / "results" / "deck_data.json").read_text())
    recs = dd.load("steady2_mc2")
    per = [dd.per_run_outputs(r["lnw"]) for r in recs]
    beta = np.array([p["beta_all"] for p in per])
    L = ["# Audit: current values of the deck numbers", "",
         "Deck built 11 Oct 2026 16:45 from results/deck_data.json (11 Oct 16:42, analysis file %s, reruns used in the run "
         "list: %d). Current = audit recomputation with the main rule of PROTOKOLL 7.65." % (D["analysis_file"],
                                                                                         sum(r["folder"] != r["name"] for r in D["runs"])), ""]
    old_beta = np.array([r["beta_all"] for r in D["runs"]])
    L.append("- β all: deck %.2f … %.2f, mean %.2f (sd %.2f), > 2: %d of %d → current %.2f … %.2f, mean %.3f (sd %.3f), > 2: %d of %d" % (
        np.nanmin(old_beta), np.nanmax(old_beta), np.nanmean(old_beta), np.nanstd(old_beta, ddof=1), (old_beta > 2).sum(), len(old_beta),
        beta.min(), beta.max(), beta.mean(), beta.std(ddof=1), (beta > 2).sum(), len(beta)))
    # repeat scatter
    grp = {}
    for r, b in zip(recs, beta):
        tag = r["name"][:2] if r["name"].startswith("K") else ("C" if r["name"].startswith("C") else None)
        if tag:
            grp.setdefault(tag, []).append(b)
    num = sum((len(v) - 1) * np.var(v, ddof=1) for v in grp.values())
    dof = sum(len(v) - 1 for v in grp.values())
    sdp = math.sqrt(num / dof)
    sdt = beta.std(ddof=1)
    L.append("- repeat scatter of β all (pooled, %d repeats in %d groups, dof %d): deck %.2f → current %.3f; per group: %s" % (
        sum(len(v) for v in grp.values()), len(grp), dof, D["scatter"]["pooled"]["sd"], sdp,
        ", ".join("%s %.3f (n %d)" % (k, np.std(v, ddof=1), len(v)) for k, v in sorted(grp.items()))))
    L.append("- spread from the parameters √(sd² − sd_rep²): deck %.2f → current %.3f" % (
        math.sqrt(max(np.nanstd(old_beta, ddof=1) ** 2 - D["scatter"]["pooled"]["sd"] ** 2, 0)), math.sqrt(sdt ** 2 - sdp ** 2)))
    # loop scatter
    sdl = []
    for r in recs:
        _, _, summ = al.read(r["folder"])
        for s in summ["stages"]:
            t = al.tail(al.measured(s))
            if len(t) >= 3:
                sdl.append(np.std(np.log([c["w_loop"] for c in t]), ddof=1))
    L.append("- loop-to-loop scatter of P (median sd of ln P over calm tails with ≥ 3 loops): deck %.1f %% → current %.2f %% (n %d stages)" % (
        100 * D["loop_scatter"]["median"], 100 * np.median(sdl), len(sdl)))
    # main fit numbers (from the audit json)
    S = json.loads((al.HERE / "sobol_audit.json").read_text())["steady2_mc2|scatter|full"]
    L += ["", "Main fit (audit recomputation, identical to results/sobol_analysis_steady2_mc2.json to 4 digits):", ""]
    for o in ("beta_all", "lnw@50mT"):
        for f in dd.FACT:
            v = S["%s × %s" % (o, f)]
            sg = -1 if f == "lnd" else 1
            lo, hi = sorted((sg * v["ci95"][0], sg * v["ci95"][1]))
            if o == "beta_all":
                L.append("- β × %s (deck sign: A up = lnd down): %+.3f (%+.3f … %+.3f), p %.4f, class %s" % (f, sg * v["eff"], lo, hi, v["p"], v["class_p"]))
            else:
                pc = lambda x: 100 * (math.exp(x) - 1)
                L.append("- P@50 × %s: %+.1f %% (%+.0f … %+.0f %%), class %s" % (f, pc(sg * v["eff"]), pc(lo), pc(hi), v["class_p"]))
    # crossing of the lnd effect per amplitude
    X = dd.coded(recs)
    T = dd.terms(X)
    W = np.exp(np.array([r["lnw"] for r in recs]))
    vw, _, _ = sa.scatter_weights(recs, X)
    coefs = np.array([sa.irls_gamma(T, W[:, k], vw) for k in range(7)])
    hi, lo = X.copy(), X.copy()
    hi[:, 3], lo[:, 3] = 1, -1
    pa = -(dd.terms(hi) @ coefs.T - dd.terms(lo) @ coefs.T).mean(0)            # A up
    cr = None
    for i in range(6):
        if pa[i] * pa[i + 1] < 0:
            x0, x1 = dd.LNA[i], dd.LNA[i + 1]
            cr = math.exp(x0 + (x1 - x0) * pa[i] / (pa[i] - pa[i + 1]))
    L.append("- A effect on ln P per amplitude (A up): %s; sign change at ≈ %.0f mT (deck: ≈ 69 mT)" % (
        ", ".join("%g mT %+.2f" % (a, v) for a, v in zip(al.AMPS, pa)), cr))
    # mesh pair numbers
    k3 = ("MK3_dx%s_s%d_b%d", range(2031, 2037))
    ba = ("M300_dx%s_s%d_b%d", range(921, 925))
    L += ["", "Mesh pair K3 (A = 20) against base (A = 10) (K3 also differs in L_eff, r_p, φ):", ""]
    for dx in ("15", "3"):
        cells = []
        lr = []
        for b in (9, 50, 100):
            a_, c_ = np.log(mesh_P(k3[0], k3[1], b, dx)), np.log(mesh_P(ba[0], ba[1], b, dx))
            m, se, p = welch(a_, c_)
            lr.append(m)
            cells.append("%d mT × %.2f (ln %+.2f ± %.2f, p %.3f, n %d / %d)" % (b, math.exp(m), m, se, p, len(a_), len(c_)))
        c = None
        x = np.log([9, 50, 100])
        for i in range(2):
            if lr[i] * lr[i + 1] < 0:
                c = math.exp(x[i] + (x[i + 1] - x[i]) * lr[i] / (lr[i] - lr[i + 1]))
        L.append("- dx %s: %s; crossing ≈ %s mT" % ("1.5" if dx == "15" else dx, "; ".join(cells), "%.0f" % c if c else "none"))
    for dx in ("3", "15"):
        fk, fb = fit_seed(k3[0], k3[1], (9, 50, 100), dx), fit_seed(ba[0], ba[1], (9, 50, 100), dx)
        m, se, p = welch(list(fk.values()), list(fb.values()))
        L.append("- Δβ (9/50/100 fit) K3 − base, dx %s: %+.3f ± %.3f (n %d / %d, Welch p %.3f)" % (
            "1.5" if dx == "15" else dx, m, se, len(fk), len(fb), p))
    f3 = fit_seed(k3[0], k3[1], (9, 20, 35, 50, 70, 100), "3")
    f15 = fit_seed(k3[0], k3[1], (9, 20, 35, 50, 70, 100), "15")
    d = np.array([f3[s] - f15[s] for s in f3 if s in f15])
    L.append("- K3 β 9 … 100 mT: dx 3 %.3f, dx 1.5 %.3f (means of per-seed fits), paired %+.3f ± %.3f (n %d, p %.3f)" % (
        np.mean(list(f3.values())), np.mean(list(f15.values())), d.mean(), d.std(ddof=1) / math.sqrt(len(d)), len(d),
        stats.ttest_1samp(d, 0).pvalue))
    (al.HERE / "deck_check.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
