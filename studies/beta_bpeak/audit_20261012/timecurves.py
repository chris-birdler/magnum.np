"""Task 6: time curves and field data (audit 2026-10-12).

(A) integrity: W_loop and W_dis of every cycle recomputed from samples_*.csv (trapezoid, 257 samples) against summary.json
(B) loop shape per stage (kept cycles of the main rule): fundamental of H and <M> -> core permeability mu_r = 1 + |M1|/|H1|,
    loss angle delta (sin delta), harmonics |M2|/|M1|, |M3|/|M1|; identity check w = pi mu0 |H1||M1| sin delta;
    per run the slopes over ln B of ln sin(delta) and ln mu_r: beta = 2 + s_delta - s_mu (+ small B-hat/B1 term)
(C) per-particle loops: w_i = mu0 * loop integral of H dMp_i (per particle volume); particles grouped by the angle of their
    stress axis to the drive (cos theta 0.795 = 37 deg, 0.188 = 79 deg)
(D) local loss measures: R1 from acc.json (all cycles >= 2, events included) against R1 from phasemap.json (16 phases,
    last 2 cycles); shares core / shell / texture top 10 % / p top 10 %; factor effects (OLS HC3)
(E) switched volume V_sw (domains.json) and its exponent against beta
(F) events: energy sign of the events in the measured cycles of the 116 original design runs
(G) damping test (C#_a001 vs C#): loss ratio per amplitude, beta difference, with all cycles and with the steady rule >= 1
Output: timecurves.md, timecurves.json, figures/tc_*.png
"""
import json
import math

import numpy as np
from scipy import stats

import audit_lib as al
import design_data as dd

MU0 = 4e-7 * math.pi


def read_csv(path):
    a = np.loadtxt(path, delimiter=",", comments="#")
    return a


def fourier(y, k):
    """complex amplitude of harmonic k over one period (first 256 of the 257 samples): y ~ Re(c e^{i k w t})."""
    n = len(y)
    th = 2 * np.pi * np.arange(n) / n
    return 2.0 / n * np.sum(y * np.exp(-1j * k * th))


def stage_curves(run, si, stname, kept_cycles):
    a = read_csv(al.RUNS / run / ("samples_%02d_%s.csv" % (si, stname)))
    out = {}
    for c in kept_cycles:
        r = a[a[:, 0] == c]
        out[c] = r
    return a, out


def analyse_run(name):
    use = al.use_name(name)
    cfg, summ, st = al.stages(use, rule="steady", min_calm=2)
    stg = [s for s in summ["stages"] if s.get("measure")]
    cos_p = np.array(cfg["anisotropy"]["cos_theta_p"]) if isinstance(cfg["anisotropy"], dict) else None
    res = {"name": name, "folder": use, "stages": []}
    integ = []
    for si, s in enumerate(summ["stages"]):
        if not s.get("measure"):
            continue
        meas = al.measured(s)
        kept = al.tail(meas)
        a = read_csv(al.RUNS / use / ("samples_%02d_%s.csv" % (si, s["name"])))
        for c in s["cycles"]:
            r = a[a[:, 0] == c["cycle"]]
            Wl = np.trapezoid(r[:, 2], r[:, 7])
            Wd = np.trapezoid(r[:, 8], r[:, 1])
            integ.append((abs(Wl - c["W_loop"]) / abs(c["W_loop"]), abs(Wd - c["W_dis"]) / abs(c["W_dis"])))
        H1s, M1s, M2s, M3s, wl, wp, mp1 = [], [], [], [], [], [], []
        for c in kept:
            r = a[a[:, 0] == c["cycle"]]
            H, M, Mp = r[:256, 2], r[:256, 3], r[:, 9:13]
            H1, M1 = fourier(H, 1), fourier(M, 1)
            H1s.append(H1)
            M1s.append(M1)
            M2s.append(abs(fourier(M - M.mean(), 2)))
            M3s.append(abs(fourier(M, 3)))
            wl.append(c["W_loop"])
            wp.append([MU0 * np.trapezoid(r[:, 2], Mp[:, i]) for i in range(4)])          # mu0 * integral of H dMp_i over the cycle
            mp1.append([abs(fourier(Mp[:256, i], 1)) for i in range(4)])
        H1, M1 = np.mean(H1s), np.mean(M1s)
        chi = M1 / H1
        rec = {"B": 1000 * cfg["Js"] * np.mean([c["b_peak"] for c in kept]), "mu_r": 1 + abs(chi),
               "sin_d": -math.sin(np.angle(chi)), "M2": np.mean(M2s) / abs(M1), "M3": np.mean(M3s) / abs(M1),
               "W_loop": np.mean(wl), "W_fund": math.pi * MU0 * abs(H1) * abs(M1) * (-math.sin(np.angle(chi))),
               "w_part": np.mean(wp, 0), "chi_part": np.mean(mp1, 0) / abs(H1), "n_kept": len(kept),
               "n_events_meas": int(sum(al.is_event(c) for c in meas))}
        res["stages"].append(rec)
    res["stages"].sort(key=lambda r: r["B"])
    res["integ"] = np.max(np.array(integ), 0)
    res["cos_p"] = cos_p
    res["cfg"] = {k: cfg[k] for k in ("r_p", "Leff_lex", "phi_vox", "d_lex_used")}
    return res


def ols(y, X):
    import statsmodels.api as sm
    return sm.OLS(y, X).fit(cov_type="HC3")


def main():
    recs = dd.load("steady2_mc2")
    names = [r["name"] for r in recs]
    X = dd.coded(recs)
    D = np.column_stack([np.ones(len(X)), X, X[:, 0] * X[:, 3]])
    R = [analyse_run(n) for n in names]
    L = ["# Audit: time curves and field data (116 design runs; kept cycles of the main rule unless stated)", ""]
    js = {}
    integ = np.array([r["integ"] for r in R])
    L += ["## (A) Integrity of the loop integrals", "",
          "W_loop and W_dis of every cycle (all stages, all cycles) recomputed from samples_*.csv: largest relative "
          "difference to summary.json %.1e (W_loop) and %.1e (W_dis) over the 116 run folders used." % tuple(integ.max(0)), ""]
    # (B) loop shape
    amps = al.AMPS
    mu = np.array([[s["mu_r"] for s in r["stages"]] for r in R])
    sd = np.array([[s["sin_d"] for s in r["stages"]] for r in R])
    M3 = np.array([[s["M3"] for s in r["stages"]] for r in R])
    M2 = np.array([[s["M2"] for s in r["stages"]] for r in R])
    Bm = np.array([[s["B"] for s in r["stages"]] for r in R])
    idc = np.array([[s["W_fund"] / s["W_loop"] for s in r["stages"]] for r in R])
    L += ["## (B) Loop shape: permeability, loss angle, harmonics (median over 116 runs; quartiles)", "",
          "| B (mT) | μ_r | sin δ | δ (deg) | \\|M2\\|/\\|M1\\| | \\|M3\\|/\\|M1\\| | π μ0 H1 M1 sin δ / W_loop |", "|---|---|---|---|---|---|---|"]
    for k, b in enumerate(amps):
        q = lambda v: "%.3g (%.3g … %.3g)" % (np.median(v), *np.percentile(v, [25, 75]))
        L.append("| %g | %s | %s | %.2f | %s | %s | %.4f … %.4f |" % (b, q(mu[:, k]), q(sd[:, k]), math.degrees(math.asin(np.median(sd[:, k]))),
                                                                  q(M2[:, k]), q(M3[:, k]), idc[:, k].min(), idc[:, k].max()))
    lnB = np.log(Bm)
    s_d = np.array([np.polyfit(lnB[i], np.log(sd[i]), 1)[0] for i in range(len(R))])
    s_m = np.array([np.polyfit(lnB[i], np.log(mu[i]), 1)[0] for i in range(len(R))])
    lnP = np.array([[math.log(s["W_loop"]) for s in r["stages"]] for r in R])
    bfull = np.array([np.polyfit(lnB[i], lnP[i], 1)[0] for i in range(len(R))])
    resid = bfull - (2 + s_d - s_m)
    L += ["", "Decomposition β = 2 + s_δ − s_μ (s = slope over ln B, 9 … 150 mT, per run): s_δ mean %+.3f (sd %.3f), "
          "s_μ mean %+.4f (sd %.4f); remainder β − (2 + s_δ − s_μ): mean %+.4f, max |.| %.4f (the peak B differs slightly "
          "from the fundamental). Correlation of β with s_δ: r = %.3f; with s_μ: r = %.3f." % (
              s_d.mean(), s_d.std(ddof=1), s_m.mean(), s_m.std(ddof=1), resid.mean(), np.abs(resid).max(),
              np.corrcoef(bfull, s_d)[0, 1], np.corrcoef(bfull, s_m)[0, 1]), ""]
    L += ["Factor effects (OLS HC3, range effect = 2 × coefficient of the coded factor; last column: lnd×lnQ interaction):", "",
          "| quantity | lnQ | r_p | φ | lnd | lnd×lnQ | R² |", "|---|---|---|---|---|---|---|"]
    qty = {"β (full)": bfull, "s_δ": s_d, "−s_μ": -s_m, "ln μ_r @9": np.log(mu[:, 0]), "ln μ_r @150": np.log(mu[:, 6]),
           "ln sin δ @9": np.log(sd[:, 0]), "ln sin δ @50": np.log(sd[:, 3]), "ln sin δ @150": np.log(sd[:, 6]),
           "|M3|/|M1| @150": M3[:, 6]}
    for k, v in qty.items():
        m = ols(v, D)
        L.append("| %s | %s | %.2f |" % (k, " | ".join("%+.3f ± %.3f" % (2 * m.params[i], 2 * m.bse[i]) for i in range(1, 6)), m.rsquared))
    js["B"] = {"mu_median": np.median(mu, 0).tolist(), "sin_median": np.median(sd, 0).tolist(), "s_d": s_d.tolist(), "s_m": s_m.tolist()}
    # (C) per particle
    L += ["", "## (C) Loss per particle by stress-axis angle (runs with r_p > 0.5; per particle volume)", "",
          "Ratio Q = (w of the 2 particles at 37°) / (w of the 2 particles at 79°) per stage; mean of ln Q ± SE over runs; "
          "slope of ln Q against r_p (OLS). Also the fundamental susceptibility ratio.", "",
          "| B (mT) | n runs | ln Q mean ± SE | Q | slope d lnQ / d r_p ± SE | χ1 ratio 37°/79° (median) | check Σ particles / box |",
          "|---|---|---|---|---|---|---|"]
    for k, b in enumerate(amps):
        lq, rp, cr, chk = [], [], [], []
        for r in R:
            if r["cfg"]["r_p"] < 0.5:
                continue
            cp = r["cos_p"]
            hi, lo = np.abs(cp) > 0.5, np.abs(cp) < 0.5
            wpp = r["stages"][k]["w_part"]
            if np.all(wpp > 0):
                lq.append(math.log(wpp[hi].sum() / wpp[lo].sum()))
                rp.append(r["cfg"]["r_p"])
                cr.append(r["stages"][k]["chi_part"][hi].mean() / r["stages"][k]["chi_part"][lo].mean())
            chk.append(r["cfg"]["phi_vox"] * wpp.mean() / r["stages"][k]["W_loop"])
        lq, rp = np.array(lq), np.array(rp)
        sl = stats.linregress(rp, lq)
        L.append("| %g | %d | %+.3f ± %.3f | %.2f | %+.3f ± %.3f | %.3f | %.3f … %.3f |" % (
            b, len(lq), lq.mean(), lq.std(ddof=1) / math.sqrt(len(lq)), math.exp(lq.mean()), sl.slope, sl.stderr,
            np.median(cr), min(chk), max(chk)))
    # (D) local measures
    L += ["", "## (D) Local loss measures (acc.json: all cycles ≥ 2 incl. events; phasemap.json: 16 phases of the last 2 cycles)", ""]
    rows = []
    for r, rec in zip(R, recs):
        acc = json.loads((al.RUNS / r["name"] / "acc.json").read_text())["stages"]
        pm = json.loads((al.RUNS / r["name"] / "phasemap.json").read_text())["stages"]
        _, _, sm = al.read(r["name"])
        msr = [s for s in sm["stages"] if s.get("measure")]
        ev = [sum(al.is_event(c) for c in al.measured(s)) for s in msr]
        accs = sorted(acc, key=lambda s: s["stage"])
        rows.append({"acc": accs, "pm": {p["stage"]: p for p in pm}, "ev": ev})
    L += ["| B (mT) | R1 acc median | R1 acc, stages without events | R1 phasemap median | R1 acc − phasemap (median, n) | "
          "core share | shell share | texture-top-10 % share | p-top-10 % share |", "|---|---|---|---|---|---|---|---|---|"]
    for k, (b, stn) in enumerate(zip(amps, ["amp%02d" % i for i in range(7)])):
        r1 = np.array([x["acc"][k]["R1"] for x in rows])
        noev = np.array([x["ev"][k] == 0 for x in rows])
        pmv = np.array([x["pm"][stn]["R1"] if stn in x["pm"] else np.nan for x in rows])
        dif = r1 - pmv
        wh = lambda key: np.median([x["acc"][k]["where"][key]["p_share"] for x in rows])
        L.append("| %g | %.3f | %.3f (n %d) | %s | %s | %.3f | %.3f | %.3f | %.3f |" % (
            b, np.median(r1), np.median(r1[noev]), noev.sum(), "%.3f" % np.nanmedian(pmv) if np.isfinite(pmv).any() else "–",
            "%+.3f (%d)" % (np.nanmedian(dif), np.isfinite(dif).sum()) if np.isfinite(dif).any() else "–",
            wh("r<0.5"), wh("r>=0.8"), wh("texture top 10%"), wh("p top 10%")))
    # effects on phasemap R1 (calm last cycles)
    L += ["", "Factor effects on R1 (phasemap, last 2 cycles; OLS HC3, range effects):", "",
          "| B | lnQ | r_p | φ | lnd | lnd×lnQ | R² |", "|---|---|---|---|---|---|---|"]
    for stn, b in (("amp00", 9), ("amp03", 50), ("amp06", 150)):
        y = np.array([x["pm"][stn]["R1"] for x in rows])
        m = ols(y, D)
        L.append("| %d mT | %s | %.2f |" % (b, " | ".join("%+.3f ± %.3f" % (2 * m.params[i], 2 * m.bse[i]) for i in range(1, 6)), m.rsquared))
    # (E) switched volume
    L += ["", "## (E) Switched volume per half cycle (domains.json, last cycle at 9 / 50 / 150 mT) against β", ""]
    vs, ex = [], []
    for r in R:
        dom = json.loads((al.RUNS / r["name"] / "domains.json").read_text())
        v = []
        for b in (9, 50, 150):
            st = min(dom["stages"], key=lambda s: abs(math.log(s["B_peak_mT"] / b)))
            v.append(np.mean(st["half_cycle"]["V_sw"]) / st["states"][0]["V"])
        vs.append(v)
        ex.append(np.polyfit(np.log([9, 50, 150]), np.log(v), 1)[0])
    vs, ex = np.array(vs), np.array(ex)
    m1 = ols(bfull, D)
    m2 = ols(bfull, np.column_stack([D, (ex - ex.mean()) / ex.std()]))
    L += ["- V_sw/V: median %.2f / %.2f / %.2f %% at 9 / 50 / 150 mT; exponent median %.2f (%.2f … %.2f)." % (
        *(100 * np.median(vs, 0)), np.median(ex), ex.min(), ex.max()),
          "- Correlation β_full vs V_sw exponent: r = %.3f. Adding the standardised exponent to the factor model: lnd effect "
          "%+.3f → %+.3f, interaction %+.3f → %+.3f, exponent coefficient %+.3f ± %.3f (per sd), R² %.2f → %.2f." % (
              np.corrcoef(bfull, ex)[0, 1], 2 * m1.params[4], 2 * m2.params[4], 2 * m1.params[5], 2 * m2.params[5],
              m2.params[6], m2.bse[6], m1.rsquared, m2.rsquared),
          "- Caveat: the snapshots of the reruns stay from the originals (13 runs; 3 of them at 9 mT in event cycles)."]
    # (F) events energy sign in the 116 originals
    rel, upt = 0, 0
    for n in names:
        _, _, sm = al.read(n)
        for s in sm["stages"]:
            if s.get("measure"):
                for c in s["cycles"][2:]:
                    if al.is_event(c):
                        if c["w_loop"] <= 0 or c["w_dis"] > c["w_loop"]:
                            rel += 1
                        else:
                            upt += 1
    L += ["", "## (F) Events in the measured cycles (cycles 2 … 6) of the 116 original design runs", "",
          "- release of stored energy (w_dis > w_loop or w_loop ≤ 0): %d; uptake: %d." % (rel, upt)]
    # (G) damping test
    L += ["", "## (G) Damping test α 0.01 vs α 0.1 (C0 … C7, same start state)", "",
          "| rule | " + " | ".join("%g mT" % b for b in amps) + " | Δβ (7-amplitude fit) |", "|---|" + "---|" * 8]
    for rule, mc in (("all", 1), ("steady", 1), ("steady", 2)):
        lr, db = [], []
        for i in range(8):
            _, _, a1 = al.stages("C%d_a001" % i, rule=rule, min_calm=mc)
            _, _, a0 = al.stages("C%d" % i, rule=rule, min_calm=mc)
            p1, p0 = np.array([s["P"] for s in a1]), np.array([s["P"] for s in a0])
            lr.append(np.log(p1 / p0))
            ok1, ok0 = np.isfinite(p1), np.isfinite(p0)
            b1 = al.powerfit(np.log([s["B_mT"] for s in a1])[ok1], np.log(p1[ok1]))[0] if ok1.sum() >= 3 else np.nan
            b0 = al.powerfit(np.log([s["B_mT"] for s in a0])[ok0], np.log(p0[ok0]))[0] if ok0.sum() >= 3 else np.nan
            db.append(b1 - b0)
        lr = np.array(lr)
        cells = []
        for k in range(7):
            v = lr[:, k][np.isfinite(lr[:, k])]
            cells.append("× %.2f (ln %+.2f ± %.2f, n %d)" % (math.exp(v.mean()), v.mean(), v.std(ddof=1) / math.sqrt(len(v)), len(v))
                         if len(v) > 1 else "n %d" % len(v))
        db = np.array(db)
        dbf = db[np.isfinite(db)]
        L.append("| %s ≥%d | %s | %+.2f ± %.2f (n %d; max %+.2f) |" % (rule, mc, " | ".join(cells), dbf.mean(),
                                                                      dbf.std(ddof=1) / math.sqrt(len(dbf)), len(dbf), dbf.max()))
    (al.HERE / "timecurves.md").write_text("\n".join(L) + "\n")
    (al.HERE / "timecurves.json").write_text(json.dumps(js, indent=1, default=float))
    print("\n".join(L))
    # figure: mu and sin delta vs B, by A thirds
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    A = np.array([r["A_pJm"] for r in recs])
    q1, q2 = np.quantile(A, [1 / 3, 2 / 3])
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
    for sel, lab in ((A <= q1, "A ≤ %.1f pJ/m" % q1), ((A > q1) & (A <= q2), "A %.1f … %.1f" % (q1, q2)), (A > q2, "A > %.1f" % q2)):
        axs[0].plot(amps, np.median(mu[sel], 0), "o-", label=lab)
        axs[1].plot(amps, np.median(sd[sel], 0), "o-", label=lab)
        axs[2].plot(amps, np.median(M3[sel], 0), "o-", label=lab)
    for ax, yl in zip(axs, ("core permeability μ_r (median)", "sin δ of the fundamental (median)", "|M3|/|M1| (median)")):
        ax.set_xscale("log")
        ax.set_xlabel("B (mT)")
        ax.set_ylabel(yl)
        ax.legend(fontsize=8)
    axs[1].set_yscale("log")
    fig.tight_layout()
    fig.savefig(al.HERE / "figures" / "tc_loopshape.png", dpi=150)


if __name__ == "__main__":
    main()
