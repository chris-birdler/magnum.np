"""Task 5: dx 3 against dx 1.5 from all mesh pairs (audit 2026-10-12).

Valid pairs (converged relaxation, steady loops >= 2, reruns replace originals):
  base   mesh300: d/l_ex 300, L_eff 12, r_p 0, phi 0.65; seeds 921 ... 924; 9 / 50 / 100 mT
  K3     d/l_ex 212, L_eff 18, r_p 1.5, phi 0.69; seeds 2031 ... 2036; 9 / 20 / 35 / 50 / 70 / 100 mT (incl. fill-in)
  cg     same-state test: dx 3 started from the coarse-grained dx 1.5 state (K3, 9 / 50 / 100 mT)
Hints only (relaxation not converged, PROTOKOLL 7.13): R212 (d/l_ex 212, base microstructure, 9 / 50 mT, seeds 911 ... 914),
  P1_d300 (d/l_ex 300, one job per amplitude, 9 / 50 mT, seeds 1 ... 3), P3_dx3 / R_dx15 (d/l_ex 96, seeds 3 ... 5).
Mesh error = value(dx 3) - value(dx 1.5), paired per seed. Errors: mean ± SE over seeds, t-CI with n - 1 dof.
Output: mesh_audit.md / .json, figures/mesh_*.png
"""
import json
import math

import numpy as np
from scipy import stats

import audit_lib as al

POINTS = {
    "base": ("M300_dx%s_s%d_b%d", list(range(921, 925)), (9, 50, 100)),
    "K3": ("MK3_dx%s_s%d_b%d", list(range(2031, 2037)), (9, 20, 35, 50, 70, 100)),
}


def job_value(name, rule="steady", min_calm=2):
    """(P W/cm^3, B mT, n_kept, missing) of a one-amplitude job, reruns replace originals; None if absent."""
    use = al.use_name(name)
    d = al.RUNS / use
    if not (d / "DONE").exists():
        return None
    cfg, summ, st = al.stages(use, rule=rule, min_calm=min_calm)
    s = st[0]
    return {"P": s["P"], "B": s["B_mT"], "n": s["n_kept"], "miss": s["missing"], "E0": summ["init"].get("E_over_Kd"),
            "n60_0": summ["init"].get("n_pairs_gt60"), "folder": use}


def tci(v, q=0.95):
    v = np.asarray([x for x in v if np.isfinite(x)])
    n = len(v)
    if n < 2:
        return (float(v.mean()) if n else float("nan"), float("nan"), n, (float("nan"), float("nan")))
    m, se = float(v.mean()), float(v.std(ddof=1) / math.sqrt(n))
    t = stats.t.ppf(0.5 + q / 2, n - 1)
    return m, se, n, (m - t * se, m + t * se)


def welch(a, b, q=0.95):
    a, b = np.asarray(a), np.asarray(b)
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    se = math.sqrt(va + vb)
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    m = a.mean() - b.mean()
    t = stats.t.ppf(0.5 + q / 2, df)
    p = 2 * stats.t.sf(abs(m / se), df)
    return m, se, df, (m - t * se, m + t * se), p


def fmt(m, se, n, ci):
    return "%+.3f ± %.3f (n %d, 95 %% CI %+.2f … %+.2f)" % (m, se, n, ci[0], ci[1])


def main():
    L = ["# Audit: mesh dx 3 against dx 1.5", ""]
    js = {}
    data = {}
    for pt, (fmt_, seeds, amps) in POINTS.items():
        V = {(dx, s, b): job_value(fmt_ % (dx, s, b)) for dx in ("3", "15") for s in seeds for b in amps}
        data[pt] = (V, seeds, amps)
    # --- level
    L += ["## Loss density level: Δ ln P = ln P(dx 3) − ln P(dx 1.5), paired per seed", "",
          "| point | B (mT) | P dx 1.5 geo. mean (W/cm³) | P dx 3 geo. mean | Δ ln P per seed | mean ± SE (n, 95 % CI) | factor |",
          "|---|---|---|---|---|---|---|"]
    lev = {}
    for pt, (V, seeds, amps) in data.items():
        for b in amps:
            d = [math.log(V[("3", s, b)]["P"] / V[("15", s, b)]["P"]) for s in seeds
                 if V[("3", s, b)] and V[("15", s, b)] and not V[("3", s, b)]["miss"] and not V[("15", s, b)]["miss"]]
            g15 = math.exp(np.mean([math.log(V[("15", s, b)]["P"]) for s in seeds if V[("15", s, b)] and not V[("15", s, b)]["miss"]]))
            g3 = math.exp(np.mean([math.log(V[("3", s, b)]["P"]) for s in seeds if V[("3", s, b)] and not V[("3", s, b)]["miss"]]))
            m, se, n, ci = tci(d)
            lev[(pt, b)] = d
            L.append("| %s | %d | %.1f | %.1f | %s | %s | × %.2f (%.2f … %.2f) |" % (pt, b, g15, g3, " ".join("%+.2f" % x for x in d),
                                                                        fmt(m, se, n, ci), math.exp(m), math.exp(ci[0]), math.exp(ci[1])))
    js["level"] = {"%s|%d" % k: v for k, v in lev.items()}
    L += ["", "Dependence on the microstructure (K3 − base, Welch):", ""]
    for b in (9, 50, 100):
        m, se, df, ci, p = welch(lev[("K3", b)], lev[("base", b)])
        L.append("- %d mT: %+.2f ± %.2f (df %.1f, 95 %% CI %+.2f … %+.2f, p %.3f)" % (b, m, se, df, ci[0], ci[1], p))
    # --- beta
    L += ["", "## β: Δβ = β(dx 3) − β(dx 1.5), per seed from power-law fits over the realised B", "",
          "| point | range | β dx 3 (mean ± SE) | β dx 1.5 | Δβ paired | 95 % CI |", "|---|---|---|---|---|---|"]
    betas = {}
    for pt, (V, seeds, amps) in data.items():
        ranges = [(9, 50, 100), (9, 50), (50, 100)]
        if pt == "K3":
            ranges = [(9, 20, 35, 50, 70, 100), (9, 50, 100), (9, 20, 35), (20, 35, 50), (35, 50, 70), (50, 70, 100), (9, 50), (50, 100),
                      (9, 20, 35, 50), (50, 70, 100)]
        for rg in ranges:
            per = {}
            for dx in ("3", "15"):
                per[dx] = {}
                for s in seeds:
                    vv = [V[(dx, s, b)] for b in rg]
                    if all(v and not v["miss"] for v in vv):
                        per[dx][s] = al.powerfit(np.log([v["B"] for v in vv]), np.log([v["P"] for v in vv]))[0]
            both = [s for s in seeds if s in per["3"] and s in per["15"]]
            d = [per["3"][s] - per["15"][s] for s in both]
            m, se, n, ci = tci(d)
            b3, b15 = tci(list(per["3"].values())), tci(list(per["15"].values()))
            betas[(pt, rg)] = per
            L.append("| %s | %s | %.2f ± %.2f (%d) | %.2f ± %.2f (%d) | %+.3f ± %.3f (n %d) | %+.2f … %+.2f |" % (
                pt, "/".join(map(str, rg)), b3[0], b3[1], b3[2], b15[0], b15[1], b15[2], m, se, n, ci[0], ci[1]))
    js["beta"] = {"%s|%s" % (k[0], "/".join(map(str, k[1]))): {dx: {str(s): v for s, v in per[dx].items()} for dx in per}
                  for k, per in betas.items()}
    # --- K3 - base contrast (A effect, confounded) on each mesh
    L += ["", "## K3 − base on each mesh (the 'A lever' pair; K3 also differs in L_eff 18 vs 12, r_p 1.5 vs 0, φ 0.69 vs 0.65)", "",
          "| quantity | dx 3 | dx 1.5 | mesh share (dx 3 − dx 1.5 of the contrast) |", "|---|---|---|---|"]
    cont = {}
    for rg in ((9, 50, 100), (9, 50), (50, 100)):
        r = {}
        for dx in ("3", "15"):
            r[dx] = welch(list(betas[("K3", rg)][dx].values()), list(betas[("base", rg)][dx].values()))
        ms = r["3"][0] - r["15"][0]
        mse_ = math.hypot(r["3"][1], r["15"][1])
        cont[rg] = r
        L.append("| β %s mT | %+.2f ± %.2f (95 %% CI %+.2f … %+.2f) | %+.2f ± %.2f (%+.2f … %+.2f) | %+.2f ± %.2f |" % (
            "/".join(map(str, rg)), r["3"][0], r["3"][1], *r["3"][3], r["15"][0], r["15"][1], *r["15"][3], ms, mse_))
    for b in (9, 50, 100):
        r = {}
        for dx in ("3", "15"):
            V3 = data["K3"][0]
            Vb = data["base"][0]
            a = [math.log(V3[(dx, s, b)]["P"]) for s in data["K3"][1] if V3[(dx, s, b)] and not V3[(dx, s, b)]["miss"]]
            c = [math.log(Vb[(dx, s, b)]["P"]) for s in data["base"][1] if Vb[(dx, s, b)] and not Vb[(dx, s, b)]["miss"]]
            r[dx] = welch(a, c)
        L.append("| ln P %d mT | %+.2f ± %.2f (× %.2f, p %.3f) | %+.2f ± %.2f (× %.2f, p %.3f) | %+.2f ± %.2f |" % (
            b, r["3"][0], r["3"][1], math.exp(r["3"][0]), r["3"][4], r["15"][0], r["15"][1], math.exp(r["15"][0]), r["15"][4],
            r["3"][0] - r["15"][0], math.hypot(r["3"][1], r["15"][1])))
    rr = cont[(9, 50, 100)]
    ratio = rr["15"][0] / rr["3"][0]
    # delta method SE of the ratio
    rse = abs(ratio) * math.hypot(rr["15"][1] / rr["15"][0], rr["3"][1] / rr["3"][0])
    L += ["", "Fine/coarse ratio of the β contrast (9/50/100 mT): %.2f ± %.2f (delta method). Applied to the design range "
          "effect of d/l_ex on β all (−0.572): %+.2f ± %.2f (estimate; assumes the mesh error scales the whole contrast)." % (
              ratio, rse, -0.572 * ratio, 0.572 * rse), ""]
    js["k3_base_ratio"] = [ratio, rse]
    # --- same-state test
    L += ["## Same-state test (MK3cg: dx 3 started from the coarse-grained dx 1.5 state)", "",
          "| B (mT) | P cg geo. mean | ln[P(cg)/P(dx 1.5)] | ln[P(cg)/P(dx 3 own)] | first-loop w_dis/w_loop cg | own dx 3 |", "|---|---|---|---|---|---|"]
    V = data["K3"][0]
    for b in (9, 50, 100):
        cg = {s: job_value("MK3cg_dx3_s%d_b%d" % (s, b)) for s in data["K3"][1]}
        d15 = [math.log(cg[s]["P"] / V[("15", s, b)]["P"]) for s in cg if cg[s] and not cg[s]["miss"] and V[("15", s, b)] and not V[("15", s, b)]["miss"]]
        d3 = [math.log(cg[s]["P"] / V[("3", s, b)]["P"]) for s in cg if cg[s] and not cg[s]["miss"] and V[("3", s, b)] and not V[("3", s, b)]["miss"]]
        bal_cg, bal_own = [], []
        for s in data["K3"][1]:
            _, _, sm = al.read("MK3cg_dx3_s%d_b%d" % (s, b))
            c0 = sm["stages"][0]["cycles"][0]
            bal_cg.append(c0["w_dis"] / c0["w_loop"])
            _, _, sm = al.read(al.use_name("MK3_dx3_s%d_b%d" % (s, b)))
            c0 = sm["stages"][0]["cycles"][0]
            bal_own.append(c0["w_dis"] / c0["w_loop"])
        a, b_ = tci(d15), tci(d3)
        L.append("| %d | %.1f | %s | %s | %.0f … %.0f | %.1f … %.1f |" % (
            b, math.exp(np.mean([math.log(cg[s]["P"]) for s in cg if cg[s] and not cg[s]["miss"]])), fmt(*a), fmt(*b_),
            min(bal_cg), max(bal_cg), min(bal_own), max(bal_own)))
    # relaxed state energy per mesh
    L += ["", "Relaxed start states (E/K_d, pairs > 60°), min … max over seeds:", ""]
    for pt, (V, seeds, amps) in data.items():
        for dx in ("3", "15"):
            e = [V[(dx, s, 9)]["E0"] for s in seeds if V[(dx, s, 9)]]
            n60 = [V[(dx, s, 9)]["n60_0"] for s in seeds if V[(dx, s, 9)]]
            L.append("- %s dx %s: E/K_d %.5f … %.5f, n60 %d … %d" % (pt, "1.5" if dx == "15" else dx, min(e), max(e), min(n60), max(n60)))
    # --- old pairs (hints)
    L += ["", "## Old pairs (relaxation not converged, 7.13): hints only; steady rule ≥ 2 and all-cycle rule", "",
          "| pairs | B | rule | Δ ln P per seed | mean ± SE (n, 95 % CI) |", "|---|---|---|---|---|"]
    old = [("R212 d/l_ex 212 (base micro.)", lambda dx, s: "R212_dx%s_s%d" % (dx, s), range(911, 915), ("15", "3")),
           ("P1 d/l_ex 300 (one job per amplitude)", None, range(1, 4), None),
           ("d/l_ex 96 (P3_dx3 vs R_dx15)", lambda dx, s: ("P3_dx3_s%d" if dx == "3" else "R_dx15_s%d") % s, range(3, 6), None)]
    for rule in ("steady", "all"):
        for lab, f, seeds, _ in old:
            for k, b in enumerate((9, 50)):
                d = []
                for s in seeds:
                    if f is None:
                        n3, n15 = "P1_d300_dx3_s%d_b%d" % (s, b), "P1_d300_s%d_b%d" % (s, b)
                        try:
                            _, _, s3 = al.stages(n3, rule=rule, min_calm=2)
                            _, _, s15 = al.stages(n15, rule=rule, min_calm=2)
                        except FileNotFoundError:
                            continue
                        v3, v15 = s3[0], s15[0]
                    else:
                        try:
                            _, _, s3 = al.stages(f("3", s), rule=rule, min_calm=2)
                            _, _, s15 = al.stages(f("15", s), rule=rule, min_calm=2)
                        except FileNotFoundError:
                            continue
                        v3 = min(s3, key=lambda x: abs(math.log(x["B_mT"] / b)))
                        v15 = min(s15, key=lambda x: abs(math.log(x["B_mT"] / b)))
                    if np.isfinite(v3["P"]) and np.isfinite(v15["P"]) and v3["P"] > 0 and v15["P"] > 0:
                        d.append(math.log(v3["P"] / v15["P"]))
                m, se, n, ci = tci(d)
                L.append("| %s | %d | %s | %s | %s |" % (lab, b, rule, " ".join("%+.2f" % x for x in d), fmt(m, se, n, ci)))
    (al.HERE / "mesh_audit.md").write_text("\n".join(L) + "\n")
    (al.HERE / "mesh_audit.json").write_text(json.dumps(js, indent=1, default=float))
    print("\n".join(L))
    fig(data, lev)


def fig(data, lev):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    ax = axs[0]
    for pt, col in (("base", "#1F4E79"), ("K3", "#9B111E")):
        V, seeds, amps = data[pt]
        for dx, ls in (("3", "--"), ("15", "-")):
            g = [math.exp(np.mean([math.log(V[(dx, s, b)]["P"]) for s in seeds if V[(dx, s, b)] and not V[(dx, s, b)]["miss"]])) for b in amps]
            ax.loglog(amps, g, ls, marker="o", color=col, label="%s dx %s" % (pt, "1.5" if dx == "15" else "3"))
    ax.set_xlabel("B (mT)")
    ax.set_ylabel("P at 30 MHz (W/cm³), geometric mean over seeds")
    ax.legend()
    ax = axs[1]
    for pt, col in (("base", "#1F4E79"), ("K3", "#9B111E")):
        amps = data[pt][2]
        m = [np.mean(lev[(pt, b)]) for b in amps]
        e = [stats.t.ppf(0.975, len(lev[(pt, b)]) - 1) * np.std(lev[(pt, b)], ddof=1) / math.sqrt(len(lev[(pt, b)])) for b in amps]
        ax.errorbar(amps, m, yerr=e, marker="o", capsize=4, color=col, label="%s (n %d)" % (pt, len(lev[(pt, amps[0])])))
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("B (mT)")
    ax.set_ylabel("ln P(dx 3) − ln P(dx 1.5), mean ± 95 % t-CI")
    ax.legend()
    fig.tight_layout()
    fig.savefig(al.HERE / "figures" / "mesh_level.png", dpi=150)


if __name__ == "__main__":
    main()
