"""Task 3: independent re-analysis of the Sobol design (audit 2026-10-12).

Own Gamma-GLM (log link) by IRLS in numpy, own design matrix and range-effect contrast, own scatter model, bootstrap
over the runs. Variants: cycle rules (design_data.RULES), weights (scatter model / none), model (15 terms / 5 terms),
estimand (E[w] by GLM / E[ln w] by OLS of the per-run outputs, HC3).

Output: audit_20261012/sobol_audit.json and sobol_audit.md
    python sobol_audit.py [--boot 1000]
"""
import argparse
import json
import math

import numpy as np
from scipy import stats

import audit_lib as al
import design_data as dd

PRIMARY = ("beta_all", "lnw@50mT")
MARGIN = 0.15


def irls_gamma(T, y, vw=None, it=100, tol=1e-10):
    """Gamma GLM, log link: working weights = var weights, z = eta + (y - mu)/mu."""
    ok = np.isfinite(y)
    T, y = T[ok], y[ok]
    vw = np.ones(len(y)) if vw is None else vw[ok]
    b = np.linalg.lstsq(T * np.sqrt(vw)[:, None], np.log(y) * np.sqrt(vw), rcond=None)[0]
    for _ in range(it):
        eta = T @ b
        mu = np.exp(eta)
        z = eta + (y - mu) / mu
        bn = np.linalg.lstsq(T * np.sqrt(vw)[:, None], z * np.sqrt(vw), rcond=None)[0]
        if np.max(np.abs(bn - b)) < tol:
            b = bn
            break
        b = bn
    return b


def contrasts(X, kind):
    out = {}
    for j, f in enumerate(dd.FACT):
        hi, lo = X.copy(), X.copy()
        hi[:, j], lo[:, j] = 1.0, -1.0
        out[f] = (dd.terms(hi, kind), dd.terms(lo, kind))
    return out


def glm_effects(T, W, C, vw=None):
    """W [n, 7] -> coefficients per amplitude -> range effects for every output and factor."""
    coefs = np.array([irls_gamma(T, W[:, k], vw) for k in range(7)])     # [7, p]
    eff = {}
    for f, (Th, Tl) in C.items():
        dlnmu = (Th @ coefs.T - Tl @ coefs.T).mean(0)                    # [7]: range effect on ln E[w] per amplitude
        for o, wv in dd.OUT.items():
            eff[(o, f)] = float(dlnmu @ wv)
    return eff, coefs


def ols_hc3(T, y):
    ok = np.isfinite(y)
    T, y = T[ok], y[ok]
    XtX_inv = np.linalg.inv(T.T @ T)
    b = XtX_inv @ T.T @ y
    e = y - T @ b
    h = np.einsum("ij,jk,ik->i", T, XtX_inv, T)
    S = (T * (e ** 2 / (1 - h) ** 2)[:, None]).T @ T                      # HC3 meat
    cov = XtX_inv @ S @ XtX_inv
    return b, cov, e, ok


def scatter_weights(recs, X, kind="full"):
    """Scatter model of PLAN 4.0 for ln w @ 50 mT: ln(resid^2) of the OLS fit on the 4 linear factors; Holm over 4 slopes.
    Returns weights (mean 1) if a slope is significant, else None; and the p values."""
    T = dd.terms(X, kind)
    y = np.array([dd.per_run_outputs(r["lnw"])["lnw@50mT"] for r in recs])
    b, cov, e, ok = ols_hc3(T, y)
    z = np.log(e ** 2 + 1e-12)
    D = np.column_stack([np.ones(ok.sum()), X[ok]])
    bz = np.linalg.lstsq(D, z, rcond=None)[0]
    res = z - D @ bz
    s2 = res @ res / (len(z) - D.shape[1])
    se = np.sqrt(np.diag(np.linalg.inv(D.T @ D)) * s2)
    p = 2 * stats.t.sf(np.abs(bz / se), len(z) - D.shape[1])[1:]
    order = np.argsort(p)
    sig = np.zeros(4, bool)
    for rank, i in enumerate(order):
        if p[i] <= 0.05 / (4 - rank):
            sig[i] = True
        else:
            break
    w = 1.0 / np.exp(np.column_stack([np.ones(len(recs)), X]) @ bz)
    return (w / w.mean() if sig.any() else None), dict(zip(dd.FACT, p.tolist())), dict(zip(dd.FACT, bz[1:].tolist()))


def holm(pv, level=0.05):
    order = np.argsort(pv)
    rej = [False] * len(pv)
    for rank, i in enumerate(order):
        if pv[i] <= level / (len(pv) - rank):
            rej[i] = True
        else:
            break
    return rej


def analyse(recs, n_boot=1000, weights="scatter", kind="full", seed=20261012):
    X = dd.coded(recs)
    T = dd.terms(X, kind)
    C = contrasts(X, kind)
    W = np.exp(np.array([r["lnw"] for r in recs]))
    vw, sp, sb = (scatter_weights(recs, X) if weights == "scatter" else (None, None, None))
    est, coefs = glm_effects(T, W, C, vw)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        i = rng.integers(0, len(recs), len(recs))
        try:
            boots.append(glm_effects(T[i], W[i], C, None if vw is None else vw[i])[0])
        except np.linalg.LinAlgError:
            continue
    res = {}
    for k, v in est.items():
        b = np.array([bb[k] for bb in boots])
        se = float(b.std(ddof=1))
        p_norm = float(2 * stats.norm.sf(abs(v) / se))
        p_pct = float(min(1.0, 2 * min((b <= 0).mean(), (b >= 0).mean())))
        res[k] = {"eff": v, "se": se, "ci95": np.percentile(b, [2.5, 97.5]).tolist(),
                  "ci90": np.percentile(b, [5, 95]).tolist(), "p": p_norm, "p_pct": max(p_pct, 1.0 / len(b)),
                  "boot": b}
    keys = [(o, f) for o in PRIMARY for f in dd.FACT]
    for pk in ("p", "p_pct"):
        rej = holm([res[k][pk] for k in keys])
        for k, r in zip(keys, rej):
            lo, hi = res[k]["ci90"]
            res[k]["class_" + pk] = "lever" if r else ("not a lever" if (-MARGIN < lo and hi < MARGIN) else "not determined")
    return res, {"weights": None if vw is None else "scatter", "scatter_p": sp, "scatter_slopes": sb, "n": len(recs),
                 "n_boot": len(boots)}


def fitA(recs, kind="full"):
    """OLS of the per-run outputs (E[ln w] estimand), HC3; range-effect contrasts."""
    X = dd.coded(recs)
    T = dd.terms(X, kind)
    C = contrasts(X, kind)
    out = {}
    per = [dd.per_run_outputs(r["lnw"]) for r in recs]
    for o in dd.OUT:
        y = np.array([p[o] for p in per])
        b, cov, e, ok = ols_hc3(T, y)
        for f, (Th, Tl) in C.items():
            c = (Th - Tl).mean(0)
            v, se = float(c @ b), float(math.sqrt(c @ cov @ c))
            out[(o, f)] = {"eff": v, "se": se, "ci95": [v - 1.96 * se, v + 1.96 * se], "p": float(2 * stats.norm.sf(abs(v) / se))}
    return out


def a_material(res, boots=True):
    """Full exchange effect at fixed d = 1 µm with Herzer K_eff ~ A^-3 (PLAN 4.3): for A x2:
    d beta = -0.5 ln2 * E_lnd / ln(300/212) - 3 ln2 * E_lnQ / ln(6.25); also for A 10 -> 20 at fixed K_eff (= -E_lnd
    scaled to ln 2 / ln(300/212)... the design range of d/l_ex is ln(300/212) = 0.347 = ln(sqrt(2)) in A)."""
    out = {}
    rd, rq = math.log(300 / 212.0), math.log(6.25)
    for o in ("beta_all", "lnw@50mT", "lnw@10mT", "lnw@100mT", "beta@9-35mT", "beta@35-70mT", "beta@70-150mT"):
        el, eq = res[(o, "lnd")], res[(o, "lnQ")]
        f = lambda a, b: -0.5 * math.log(2) * a / rd - 3 * math.log(2) * b / rq
        v = f(el["eff"], eq["eff"])
        bb = np.array([f(a, b) for a, b in zip(el["boot"], eq["boot"])])
        out[o] = {"eff": v, "ci95": np.percentile(bb, [2.5, 97.5]).tolist(), "se": float(bb.std(ddof=1))}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=1000)
    a = ap.parse_args()
    data = {r: dd.load(r) for r in dd.RULES}
    variants = [("steady2_mc2", "scatter", "full"), ("steady2_mc2", "none", "full"), ("steady2_mc2", "scatter", "lin"),
                ("steady2_mc2", "none", "lin"), ("steady2", "scatter", "full"), ("steady1", "scatter", "full"),
                ("all", "scatter", "full"), ("last2_mc2", "scatter", "full"), ("wdis_mc2", "scatter", "full")]
    allres = {}
    for rule, wt, kind in variants:
        recs = data[rule]
        if rule == "all":
            recs = [r for r in recs if r["name"] != "S029"]          # pre-registered: S029 dropped (mean w < 0 at 9 mT)
        res, info = analyse(recs, a.boot if (rule, wt, kind) == ("steady2_mc2", "scatter", "full") else max(400, a.boot // 2),
                            wt, kind)
        allres[(rule, wt, kind)] = (res, info)
        print("done", rule, wt, kind, info["weights"], info["n"], info["n_boot"], flush=True)
    A = {k: fitA(data[k]) for k in ("steady2_mc2",)}
    A_lin = fitA(data["steady2_mc2"], "lin")
    main_res = allres[("steady2_mc2", "scatter", "full")][0]
    amat = a_material(main_res)
    # ---- report
    L = ["# Audit: independent re-analysis of the Sobol design", "",
         "Own loader (audit_lib / design_data), own Gamma GLM (IRLS), own bootstrap (%d resamples for the main variant, "
         "%d for the others). Range effect = mean over the 116 design points of the prediction at the factor maximum minus "
         "at its minimum (coded ±1), as PLAN 4.0. Signs: as the reduced factors (lnQ up = K_eff up; lnd up = A down at d = 1 µm)."
         % (a.boot, max(400, a.boot // 2)), ""]
    info = allres[("steady2_mc2", "scatter", "full")][1]
    L += ["Scatter model (ln w @ 50 mT): p %s; slopes %s; weights used: %s" % (
        {k: round(v, 4) for k, v in info["scatter_p"].items()}, {k: round(v, 2) for k, v in info["scatter_slopes"].items()},
        info["weights"]), ""]
    L += ["## Primary outputs, all variants", "",
          "| output × factor | " + " | ".join("%s/%s/%s" % v for v in variants) + " |",
          "|---|" + "---|" * len(variants)]
    for o in PRIMARY:
        for f in dd.FACT:
            cells = []
            for v in variants:
                r = allres[v][0][(o, f)]
                cells.append("%+.3f (%+.3f … %+.3f) p %.4f [%s]" % (r["eff"], r["ci95"][0], r["ci95"][1], r["p"],
                                                                   r["class_p"]))
            L.append("| %s × %s | " % (o, f) + " | ".join(cells) + " |")
    L += ["", "## Main variant: all outputs (95 % percentile CI; p normal / p percentile; class with Holm on the 8 "
          "primary tests)", "", "| output | factor | effect | 95 % CI | 90 % CI | p | p_pct | class (p) | class (p_pct) | "
          "fit A (E[ln w], 15 terms, HC3) | fit A main effects only |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for o in dd.OUT:
        for f in dd.FACT:
            r = main_res[(o, f)]
            fa, fl = A["steady2_mc2"][(o, f)], A_lin[(o, f)]
            L.append("| %s | %s | %+.3f | %+.3f … %+.3f | %+.3f … %+.3f | %.2g | %.2g | %s | %s | %+.3f ± %.3f | %+.3f ± %.3f |" % (
                o, f, r["eff"], *r["ci95"], *r["ci90"], r["p"], r["p_pct"], r.get("class_p", "") if o in PRIMARY else "-",
                r.get("class_p_pct", "") if o in PRIMARY else "-", fa["eff"], fa["se"], fl["eff"], fl["se"]))
    L += ["", "## Exchange stiffness at fixed d = 1 µm with Herzer K_eff ∝ A⁻³ (A × 2; PLAN 4.3 formula, linear in the "
          "range effects; bootstrap CI)", "", "| output | effect of A × 2 | 95 % CI |", "|---|---|---|"]
    for o, v in amat.items():
        L.append("| %s | %+.3f | %+.3f … %+.3f |" % (o, v["eff"], *v["ci95"]))
    (al.HERE / "sobol_audit.md").write_text("\n".join(L) + "\n")
    js = {"%s|%s|%s" % v: {"%s × %s" % k: {kk: vv for kk, vv in r.items() if kk != "boot"} for k, r in allres[v][0].items()}
          for v in variants}
    js["fitA"] = {"%s × %s" % k: v for k, v in A["steady2_mc2"].items()}
    js["fitA_lin"] = {"%s × %s" % k: v for k, v in A_lin.items()}
    js["A_material"] = amat
    js["info"] = {"%s|%s|%s" % v: allres[v][1] for v in variants}
    (al.HERE / "sobol_audit.json").write_text(json.dumps(js, indent=1, default=float))
    np.savez_compressed(al.HERE / "sobol_audit_boot.npz", **{"%s__%s" % k: v["boot"] for k, v in main_res.items()})
    print("\n".join(L[:60]))


if __name__ == "__main__":
    main()
