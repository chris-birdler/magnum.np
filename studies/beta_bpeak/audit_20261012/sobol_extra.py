"""Task 3 extras (audit 2026-10-12):
(1) conditional range effects (interactions) from the main fit with bootstrap CIs;
(2) Monte-Carlo error of the bootstrap: r_p p value and phi 90 % edge for 5 bootstrap seeds;
(3) check: the 'without flagged runs' column of the steady analysis equals an unweighted refit;
(4) predictions of the main fit at the mesh-check points K3 and base (factor values of the mesh jobs), for task 5.
Output: sobol_extra.md / .json
"""
import json
import math

import numpy as np

import audit_lib as al
import design_data as dd
import sobol_audit as sa

OUTS = ("beta_all", "beta@9-35mT", "beta@35-70mT", "beta@70-150mT", "lnw@10mT", "lnw@50mT", "lnw@100mT")


def coefs_all(T, W, vw):
    return np.array([sa.irls_gamma(T, W[:, k], vw) for k in range(7)])


def cond_effect(coefs, X, fj, fk, lev):
    j, k = dd.FACT.index(fj), dd.FACT.index(fk)
    hi, lo = X.copy(), X.copy()
    hi[:, j], lo[:, j] = 1.0, -1.0
    hi[:, k] = lo[:, k] = lev
    d = (dd.terms(hi) @ coefs.T - dd.terms(lo) @ coefs.T).mean(0)
    return {o: float(d @ dd.OUT[o]) for o in OUTS}


def point_pred(coefs, x_real):
    """ln E[w] at the 7 amplitudes at one point (real factor values)."""
    xc = np.array([[2.0 * (x_real[f] - dd.RANGES[f][0]) / (dd.RANGES[f][1] - dd.RANGES[f][0]) - 1.0 for f in dd.FACT]])
    return (dd.terms(xc) @ coefs.T)[0]


def main(n_boot=1000):
    recs = dd.load("steady2_mc2")
    X = dd.coded(recs)
    T = dd.terms(X)
    W = np.exp(np.array([r["lnw"] for r in recs]))
    vw, _, _ = sa.scatter_weights(recs, X)
    c0 = coefs_all(T, W, vw)
    rng = np.random.default_rng(7)
    boots = []
    for _ in range(n_boot):
        i = rng.integers(0, len(recs), len(recs))
        boots.append(coefs_all(T[i], W[i], vw[i]))
    L = ["# Audit extras: interactions, bootstrap noise, flag column, mesh-point predictions", ""]
    # (1) conditional effects
    L += ["## (1) Conditional range effects (main fit, 95 %% bootstrap CI, %d resamples)" % n_boot, "",
          "Range effect of the first factor with the second factor held at its low / centre / high coded level "
          "(the other two factors averaged over the design points).", "",
          "| effect of | held | level | β all | β 9-35 | β 35-70 | β 70-150 | ln P 10 mT | ln P 50 mT | ln P 100 mT |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    js = {"cond": {}}
    for fj, fk in (("lnd", "lnQ"), ("lnd", "r_p"), ("lnd", "phi"), ("lnQ", "lnd"), ("r_p", "lnd"), ("phi", "lnd"),
                   ("r_p", "lnQ"), ("lnQ", "r_p")):
        for lev in (-1.0, 0.0, 1.0):
            e0 = cond_effect(c0, X, fj, fk, lev)
            eb = [cond_effect(c, X, fj, fk, lev) for c in boots]
            cells = []
            for o in OUTS:
                lo, hi = np.percentile([b[o] for b in eb], [2.5, 97.5])
                cells.append("%+.2f (%+.2f … %+.2f)" % (e0[o], lo, hi))
                js["cond"]["%s|%s|%+g|%s" % (fj, fk, lev, o)] = [e0[o], lo, hi]
            L.append("| %s | %s | %+g | %s |" % (fj, fk, lev, " | ".join(cells)))
        # difference high - low of the held factor (interaction contrast)
        dl = {o: cond_effect(c0, X, fj, fk, 1.0)[o] - cond_effect(c0, X, fj, fk, -1.0)[o] for o in OUTS}
        db = [{o: cond_effect(c, X, fj, fk, 1.0)[o] - cond_effect(c, X, fj, fk, -1.0)[o] for o in OUTS} for c in boots]
        cells = []
        for o in OUTS:
            lo, hi = np.percentile([b[o] for b in db], [2.5, 97.5])
            cells.append("**%+.2f (%+.2f … %+.2f)**" % (dl[o], lo, hi))
        L.append("| %s | %s | high − low | %s |" % (fj, fk, " | ".join(cells)))
    # (2) bootstrap Monte-Carlo error
    L += ["", "## (2) Bootstrap Monte-Carlo noise of the decision quantities (5 bootstrap seeds × 1000)", "",
          "| seed | β×r_p p (normal) | β×r_p p (percentile) | β×φ 90 % CI | β×lnQ p | Holm class r_p | class φ |",
          "|---|---|---|---|---|---|---|"]
    js["mc"] = []
    for s in (20261007, 1, 2, 3, 4):
        res, _ = sa.analyse(recs, 1000, "scatter", "full", seed=s)
        rp, ph, q = res[("beta_all", "r_p")], res[("beta_all", "phi")], res[("beta_all", "lnQ")]
        L.append("| %d | %.4f | %.4f | %+.3f … %+.3f | %.5f | %s | %s |" % (s, rp["p"], rp["p_pct"], *ph["ci90"], q["p"],
                                                                         rp["class_p"], ph["class_p"]))
        js["mc"].append([s, rp["p"], rp["p_pct"], ph["ci90"], q["p"], rp["class_p"], ph["class_p"]])
    # (3) flag column check
    S = json.loads((al.STUDY / "results" / "sobol_analysis_steady2_mc2.json").read_text())
    unw, _ = sa.analyse(recs, 50, "none", "full")
    L += ["", "## (3) 'without flagged runs' column of results/sobol_analysis_steady2_mc2 against an unweighted refit", "",
          "n flagged in the steady analysis: %d. If 0, the column is the main fit on all runs without the scatter weights "
          "(point estimates below)." % S["n_flagged"], "", "| output × factor | main (weighted) | 'without flagged' | my unweighted refit |",
          "|---|---|---|---|"]
    for o in ("beta@9-35mT", "beta@35-70mT", "beta@70-150mT", "lnw@10mT", "lnw@100mT"):
        for f in dd.FACT:
            k = "%s × %s" % (o, f)
            L.append("| %s | %+.3f | %+.3f | %+.3f |" % (k, S["main"][k]["eff"], S["noflag"][k]["eff"], unw[(o, f)]["eff"]))
    # (4) predictions at the mesh points
    pts = {"base (d/l_ex 300, L 12, r_p 0, φ 0.65)": None, "K3 (d/l_ex 212, L 18, r_p 1.5, φ 0.69)": None}
    for name, run in (("base (d/l_ex 300, L 12, r_p 0, φ 0.65)", "M300_dx3_s921_b9"),
                      ("K3 (d/l_ex 212, L 18, r_p 1.5, φ 0.69)", "MK3_dx3_s2031_b9")):
        _, cfg, _ = al.read(run)
        pts[name] = {"lnQ": math.log(1 / cfg["Leff_lex"] ** 2), "r_p": cfg["r_p"], "phi": cfg["phi_vox"],
                     "lnd": math.log(cfg["d_lex_used"])}
    L += ["", "## (4) Main-fit prediction at the two mesh-check points (dx 3 design model; E[P] of the design population)", "",
          "| point | P 9 / 50 / 100 mT (W/cm³) | β(9/50/100 fit of ln E[P]) | β all |", "|---|---|---|---|"]
    pf = recs[0]["pf"]
    js["pred"] = {}
    preds = {}
    for name, x in pts.items():
        lm = point_pred(c0, x)
        lb = [point_pred(c, x) for c in boots]
        sel = dd.ix(9, 50, 100)
        b3 = np.polyfit(dd.LNA[sel], lm[sel], 1)[0]
        b3b = [np.polyfit(dd.LNA[sel], v[sel], 1)[0] for v in lb]
        ba = float(lm @ dd.OUT["beta_all"])
        bab = [float(v @ dd.OUT["beta_all"]) for v in lb]
        preds[name] = (lm, lb)
        L.append("| %s | %s | %.2f (%.2f … %.2f) | %.2f (%.2f … %.2f) |" % (
            name, " / ".join("%.1f (%.1f … %.1f)" % (math.exp(lm[k]) * pf, *np.exp(np.percentile([v[k] for v in lb], [2.5, 97.5])) * pf)
                             for k in sel), b3, *np.percentile(b3b, [2.5, 97.5]), ba, *np.percentile(bab, [2.5, 97.5])))
        js["pred"][name] = {"lnP": (lm + math.log(pf)).tolist(), "beta_9_50_100": b3, "beta_all": ba,
                            "beta_9_50_100_ci": np.percentile(b3b, [2.5, 97.5]).tolist()}
    (nk, nb) = list(preds)[1], list(preds)[0]
    sel = dd.ix(9, 50, 100)
    d0 = preds[nk][0] - preds[nb][0]
    db = [a - b for a, b in zip(preds[nk][1], preds[nb][1])]
    L.append("| K3 − base | Δ ln P %s | Δβ %+.2f (%+.2f … %+.2f) | Δβ all %+.2f |" % (
        " / ".join("%+.2f (%+.2f … %+.2f)" % (d0[k], *np.percentile([v[k] for v in db], [2.5, 97.5])) for k in sel),
        np.polyfit(dd.LNA[sel], d0[sel], 1)[0], *np.percentile([np.polyfit(dd.LNA[sel], v[sel], 1)[0] for v in db], [2.5, 97.5]),
        float(d0 @ dd.OUT["beta_all"])))
    js["pred"]["K3-base"] = {"dlnP": d0.tolist(), "dbeta_9_50_100": float(np.polyfit(dd.LNA[sel], d0[sel], 1)[0]),
                             "ci": np.percentile([np.polyfit(dd.LNA[sel], v[sel], 1)[0] for v in db], [2.5, 97.5]).tolist(),
                             "dlnP_ci": [np.percentile([v[k] for v in db], [2.5, 97.5]).tolist() for k in range(7)]}
    # decomposition of the K3 - base contrast into the four factors (path base -> K3 one factor at a time, fit)
    xb, xk = pts[nb], pts[nk]
    L += ["", "Decomposition of the predicted K3 − base contrast (change one factor after the other, order lnd, lnQ, r_p, φ; "
          "path-dependent because of interactions):", ""]
    cur = dict(xb)
    prev = point_pred(c0, cur)
    for f in ("lnd", "lnQ", "r_p", "phi"):
        cur[f] = xk[f]
        now = point_pred(c0, cur)
        d = now - prev
        L.append("- %s: Δβ(9/50/100) %+.2f, Δ ln P 9/50/100 mT %s" % (f, np.polyfit(dd.LNA[sel], d[sel], 1)[0],
                                                                   " / ".join("%+.2f" % d[k] for k in sel)))
        prev = now
    (al.HERE / "sobol_extra.md").write_text("\n".join(L) + "\n")
    (al.HERE / "sobol_extra.json").write_text(json.dumps(js, indent=1, default=float))
    print("\n".join(L))


if __name__ == "__main__":
    main()
