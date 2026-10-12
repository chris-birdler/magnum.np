"""Task 4: beta from real power-law fits per run (audit 2026-10-12).

Data: the 116 design runs, steady tail with >= 2 calm cycles, reruns (main rule of PROTOKOLL 7.65); realised B
(mean of the kept cycles); P = f K_d w in W/cm^3.
Per run: least squares of ln P on ln B
  full      all 7 amplitudes: beta, SE (from the residuals, dof 5), rms residual, max |residual|
  quad      ln P = a + b x + c x^2 (x = ln B - mean): curvature dbeta/dlnB = 2c with SE, t (dof 4)
  windows   3 neighbouring amplitudes (5 windows), 4-amplitude windows, 2-point segments
  hinge     continuous broken power law with the knot at one inner amplitude (best of 20 ... 100 mT): beta below / above
  chi2      of the full fit against the cycle error of ln P (sd of the kept cycles / sqrt(n))
Comparison with the beta outputs of the study code: beta_all (nominal-interpolated ln w, equal weights, as
analyze_design.per_run_outputs) and beta_pow (analyze.features: weights 1/(max(s_cyc, 0.10)^2 / n_kept)).
Output: beta_fits.md, beta_fits.csv, figures/beta_*.png
"""
import csv
import math

import numpy as np
from scipy import stats

import audit_lib as al
import design_data as dd

WIN3 = ((9, 20, 35), (20, 35, 50), (35, 50, 70), (50, 70, 100), (70, 100, 150))
WIN4 = ((9, 20, 35, 50), (20, 35, 50, 70), (35, 50, 70, 100), (50, 70, 100, 150))


def run_fits(name):
    cfg, summ, st = al.stages(al.use_name(name), rule="steady", min_calm=2)
    B = np.array([s["B_mT"] for s in st])
    P = np.array([s["P"] for s in st])
    n = np.array([s["n_kept"] for s in st])
    sd = np.array([s["cyc_sd_ln"] for s in st])
    x, y = np.log(B), np.log(P)
    nom = [al.nominal(b) for b in B]
    out = {"name": name, "B": B, "P": P}
    b, se, res, a = al.powerfit(x, y)
    out.update(beta_full=b, se_full=se, rms=float(np.sqrt(np.mean(res ** 2))), maxres=float(np.max(np.abs(res))), res=res)
    xc = x - x.mean()
    Xq = np.c_[np.ones_like(xc), xc, xc ** 2]
    cq, *_ = np.linalg.lstsq(Xq, y, rcond=None)
    rq = y - Xq @ cq
    s2 = rq @ rq / (len(y) - 3)
    cov = s2 * np.linalg.inv(Xq.T @ Xq)
    out.update(curv=2 * cq[2], curv_se=2 * math.sqrt(cov[2, 2]), curv_t=cq[2] / math.sqrt(cov[2, 2]))
    for w in WIN3 + WIN4:
        idx = [nom.index(m) for m in w]
        out["b%d-%d%s" % (w[0], w[-1], "" if len(w) == 3 else "(4)")] = al.powerfit(x[idx], y[idx])[0]
    for k in range(6):
        out["seg%d-%d" % (nom[k], nom[k + 1])] = (y[k + 1] - y[k]) / (x[k + 1] - x[k])
    # hinge model
    best = None
    for k in range(1, 6):
        H = np.c_[np.ones_like(x), x, np.maximum(0.0, x - x[k])]
        c, *_ = np.linalg.lstsq(H, y, rcond=None)
        r = y - H @ c
        rss = float(r @ r)
        if best is None or rss < best[0]:
            best = (rss, nom[k], c[1], c[1] + c[2])
    rss1 = float(res @ res)
    out.update(hinge_B=best[1], hinge_lo=best[2], hinge_hi=best[3], hinge_F=((rss1 - best[0]) / 1) / (best[0] / 3) if best[0] > 0 else np.inf)
    # chi2 against the cycle error
    err = np.where(np.isfinite(sd) & (sd > 0), sd, np.nanmedian(sd)) / np.sqrt(n)
    out["chi2_dof"] = float(np.sum((res / err) ** 2) / 5)
    # study outputs
    rec = next(r for r in RECS if r["name"] == name)
    out["beta_all_code"] = dd.per_run_outputs(rec["lnw"])["beta_all"]
    s_cyc = math.sqrt(np.nansum(sd ** 2 * (n - 1)) / np.sum(n - 1)) if np.sum(n - 1) > 0 else 0.1
    wts = 1.0 / (max(s_cyc, 0.10) ** 2 / n)
    Xw = np.c_[np.ones_like(x), x]
    cw = np.linalg.solve(Xw.T @ (Xw * wts[:, None]), Xw.T @ (wts * y))
    out["beta_pow_code"] = float(cw[1])
    for f in dd.FACT:
        out[f] = rec[f]
    out["A_pJm"], out["K_eff"], out["Leff"] = rec["A_pJm"], rec["K_eff"], rec["Leff"]
    return out


RECS = dd.load("steady2_mc2")


def main():
    rows = [run_fits(r["name"]) for r in RECS]
    keys = [k for k in rows[0] if k not in ("B", "P", "res", "name")]
    with open(al.HERE / "beta_fits.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name"] + keys)
        for r in rows:
            w.writerow([r["name"]] + [r[k] for k in keys])
    g = lambda k: np.array([r[k] for r in rows], float)
    L = ["# Audit: β from power-law fits per run (116 design runs, steady loops ≥ 2, reruns)", ""]
    bf = g("beta_full")
    L += ["## Full range 9 … 150 mT (7 amplitudes, realised B)", "",
          "- β_full: mean %.3f, sd %.3f, SE %.3f, median %.3f, min %.2f (%s), max %.2f (%s); above 2: %d of %d." % (
              bf.mean(), bf.std(ddof=1), bf.std(ddof=1) / math.sqrt(len(bf)), np.median(bf), bf.min(),
              rows[int(bf.argmin())]["name"], bf.max(), rows[int(bf.argmax())]["name"], (bf > 2).sum(), len(bf)),
          "- Difference to the study output β_all (nominal-interpolated ln w, equal weights): max |Δ| %.4f (realised vs "
          "nominal B). To β_pow of analyze.features (weights ∝ number of kept cycles): mean %+.4f, max |Δ| %.3f." % (
              np.max(np.abs(bf - g("beta_all_code"))), np.mean(g("beta_pow_code") - bf), np.max(np.abs(g("beta_pow_code") - bf))),
          "- Fit quality: rms residual of ln P median %.3f (quartiles %.3f … %.3f), max |residual| median %.3f; "
          "SE of β_full from the residuals median %.3f." % (np.median(g("rms")), *np.percentile(g("rms"), [25, 75]),
                                                           np.median(g("maxres")), np.median(g("se_full"))),
          "- χ²/dof of the single power law against the cycle error of ln P: median %.0f (quartiles %.0f … %.0f): the cycle "
          "error (≈ 1 %%) is far smaller than the misfit; a single power law is not exact for any run." % (
              np.median(g("chi2_dof")), *np.percentile(g("chi2_dof"), [25, 75])),
          ""]
    cv, ct = g("curv"), g("curv_t")
    tcrit = stats.t.ppf(0.975, 4)
    L += ["## Curvature dβ/d ln B (quadratic fit in ln-ln, per run)", "",
          "- mean %+.3f ± %.3f (SE over runs, n %d); median %+.3f; runs with |t| > %.2f (95 %%, dof 4): %d of %d "
          "(%d negative, %d positive)." % (cv.mean(), cv.std(ddof=1) / math.sqrt(len(cv)), len(cv), np.median(cv), tcrit,
                                           (np.abs(ct) > tcrit).sum(), len(ct), ((ct < -tcrit)).sum(), (ct > tcrit).sum()),
          ""]
    # local beta per window
    L += ["## Local β in windows (per-run power-law fits; mean ± SE over the 116 runs, sd, fraction > 2)", "",
          "| window | mean ± SE | sd | median | fraction > 2 |", "|---|---|---|---|---|"]
    for w in WIN3 + WIN4:
        k = "b%d-%d%s" % (w[0], w[-1], "" if len(w) == 3 else "(4)")
        v = g(k)
        L.append("| %s mT | %.3f ± %.3f | %.3f | %.3f | %.2f |" % (k[1:], v.mean(), v.std(ddof=1) / math.sqrt(len(v)),
                                                                 v.std(ddof=1), np.median(v), (v > 2).mean()))
    for k in range(6):
        key = "seg%d-%d" % (al.AMPS[k], al.AMPS[k + 1])
        v = g(key)
        L.append("| segment %d → %d mT | %.3f ± %.3f | %.3f | %.3f | %.2f |" % (al.AMPS[k], al.AMPS[k + 1], v.mean(),
                                                                             v.std(ddof=1) / math.sqrt(len(v)), v.std(ddof=1),
                                                                             np.median(v), (v > 2).mean()))
    # paired differences between windows
    L += ["", "Paired differences of local β within the runs (mean ± SE, t-test p, n 116):", ""]
    for a_, b_ in (("b20-50", "b9-35"), ("b35-70", "b9-35"), ("b70-150", "b35-70"), ("b70-150", "b9-35"), ("b50-100", "b20-50")):
        d = g(a_) - g(b_)
        L.append("- %s − %s: %+.3f ± %.3f (p %.2g)" % (a_[1:], b_[1:], d.mean(), d.std(ddof=1) / math.sqrt(len(d)),
                                                       stats.ttest_1samp(d, 0).pvalue))
    # hinge
    hb = g("hinge_B")
    L += ["", "## Broken power law (continuous, knot at the best inner amplitude)", "",
          "- knot location (count of runs): " + ", ".join("%d mT: %d" % (b, (hb == b).sum()) for b in (20, 35, 50, 70, 100)),
          "- β below / above the knot: %.2f ± %.2f / %.2f ± %.2f (mean ± SE); F (1, 3) > %.1f (95 %%): %d of %d runs." % (
              g("hinge_lo").mean(), g("hinge_lo").std(ddof=1) / math.sqrt(116), g("hinge_hi").mean(),
              g("hinge_hi").std(ddof=1) / math.sqrt(116), stats.f.ppf(0.95, 1, 3), (g("hinge_F") > stats.f.ppf(0.95, 1, 3)).sum(), 116),
          ""]
    # dependence on factors: OLS main effects + lnd:lnQ on the per-run outputs
    import statsmodels.api as sm
    X = dd.coded(RECS)
    L += ["## Dependence of the per-run β on the factors (OLS, HC3; coded factors ±1 = range ends; range effect = 2 × slope)", "",
          "| output | lnQ (K_eff 1.0→6.2 kJ/m³) | r_p (0→3) | φ (0.55→0.69) | lnd (d/l_ex 212→300) | lnd×lnQ (interaction, coded) | R² |",
          "|---|---|---|---|---|---|---|"]
    D = np.column_stack([np.ones(len(X)), X, X[:, 0] * X[:, 3]])
    for k in ["beta_full"] + ["b%d-%d" % (w[0], w[-1]) for w in WIN3] + ["curv"]:
        m = sm.OLS(g(k), D).fit(cov_type="HC3")
        L.append("| %s | %s | %.2f |" % (k, " | ".join("%+.3f ± %.3f" % (2 * m.params[i], 2 * m.bse[i]) for i in range(1, 6)), m.rsquared))
    L += ["", "(The interaction column is 2 × the coefficient of x_lnQ·x_lnd: the change of the lnd range effect per unit of "
          "coded lnQ.)", ""]
    # beta by L_eff level and A thirds
    L += ["## Mean local β against B by groups (mean over the runs of the group)", "",
          "| group | n | " + " | ".join("%d-%d" % (w[0], w[-1]) for w in WIN3) + " | full |", "|---|---|" + "---|" * 6]
    A = g("A_pJm")
    q1, q2 = np.quantile(A, [1 / 3, 2 / 3])
    groups = [("L_eff %d l_ex" % L_, g("Leff") == L_) for L_ in (12, 18, 24, 30)]
    groups += [("A %.1f … %.1f pJ/m" % (A.min(), q1), A <= q1), ("A %.1f … %.1f" % (q1, q2), (A > q1) & (A <= q2)),
               ("A %.1f … %.1f" % (q2, A.max()), A > q2)]
    rp, ph = g("r_p"), g("phi")
    groups += [("r_p < 1", rp < 1), ("r_p 1 … 2", (rp >= 1) & (rp < 2)), ("r_p ≥ 2", rp >= 2)]
    groups += [("L_eff 30 and A > %.1f" % q2, (g("Leff") == 30) & (A > q2)), ("L_eff 30 and A ≤ %.1f" % q1, (g("Leff") == 30) & (A <= q1)),
               ("L_eff 12 and A > %.1f" % q2, (g("Leff") == 12) & (A > q2)), ("L_eff 12 and A ≤ %.1f" % q1, (g("Leff") == 12) & (A <= q1))]
    for lab, sel in groups:
        vals = ["%.2f ± %.2f" % (g("b%d-%d" % (w[0], w[-1]))[sel].mean(), g("b%d-%d" % (w[0], w[-1]))[sel].std(ddof=1) / math.sqrt(sel.sum()))
                for w in WIN3]
        L.append("| %s | %d | %s | %.2f ± %.2f |" % (lab, sel.sum(), " | ".join(vals), bf[sel].mean(), bf[sel].std(ddof=1) / math.sqrt(sel.sum())))
    (al.HERE / "beta_fits.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    figs(rows, g, groups)


def figs(rows, g, groups):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
    ax = axs[0]
    ax.hist(g("beta_full"), bins=np.arange(1.4, 2.9, 0.1), color="#1F4E79", edgecolor="white")
    ax.axvline(2, color="grey", ls="--")
    ax.set_xlabel("β, power-law fit 9 … 150 mT (per run)")
    ax.set_ylabel("number of runs")
    ax.set_title("116 design runs, dx 3, steady loops")
    ax = axs[1]
    xc = [math.exp(np.mean(np.log(w))) for w in WIN3]
    for lab, sel in groups[:4]:
        m = [g("b%d-%d" % (w[0], w[-1]))[sel].mean() for w in WIN3]
        e = [2 * g("b%d-%d" % (w[0], w[-1]))[sel].std(ddof=1) / math.sqrt(sel.sum()) for w in WIN3]
        ax.errorbar(xc, m, yerr=e, marker="o", capsize=3, label="%s (n %d)" % (lab, sel.sum()))
    ax.set_xscale("log")
    ax.axhline(2, color="grey", ls=":")
    ax.set_xlabel("window centre B (mT)")
    ax.set_ylabel("local β (3-point fit), mean ± 2 SE")
    ax.legend(fontsize=8)
    ax.set_title("by Herzer cell L_eff")
    ax = axs[2]
    for lab, sel in groups[4:7]:
        m = [g("b%d-%d" % (w[0], w[-1]))[sel].mean() for w in WIN3]
        e = [2 * g("b%d-%d" % (w[0], w[-1]))[sel].std(ddof=1) / math.sqrt(sel.sum()) for w in WIN3]
        ax.errorbar(xc, m, yerr=e, marker="o", capsize=3, label="%s (n %d)" % (lab, sel.sum()))
    ax.set_xscale("log")
    ax.axhline(2, color="grey", ls=":")
    ax.set_xlabel("window centre B (mT)")
    ax.set_ylabel("local β (3-point fit), mean ± 2 SE")
    ax.legend(fontsize=8)
    ax.set_title("by exchange stiffness A at d = 1 µm")
    fig.tight_layout()
    fig.savefig(al.HERE / "figures" / "beta_fits.png", dpi=150)
    # interaction figure: beta_full against A for the 4 L_eff levels
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for L_ in (12, 18, 24, 30):
        sel = g("Leff") == L_
        ax.scatter(g("A_pJm")[sel], g("beta_full")[sel], s=18, label="L_eff %d l_ex (K_eff %.1f kJ/m³)" % (L_, 895.2 / L_ ** 2))
        c = np.polyfit(np.log(g("A_pJm")[sel]), g("beta_full")[sel], 1)
        xx = np.linspace(9.7, 20.5, 20)
        ax.plot(xx, np.polyval(c, np.log(xx)), lw=1.5)
    ax.set_xlabel("exchange stiffness A at d = 1 µm (pJ/m)")
    ax.set_ylabel("β, power-law fit 9 … 150 mT")
    ax.legend(fontsize=8)
    ax.set_title("A lever depends on K_eff (dx 3, 116 runs)")
    fig.tight_layout()
    fig.savefig(al.HERE / "figures" / "beta_A_by_Leff.png", dpi=150)


if __name__ == "__main__":
    main()
