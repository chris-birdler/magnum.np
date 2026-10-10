#!/usr/bin/env python
"""
Stage 1 test of the regime change near 30 … 50 mT (PROTOKOLL 7.50, 7.51): does the dissipation become jumpy, and does
the amplitude of this change follow A or K_eff? Exploratory. Source: the time series of the α 0.1 Sobol runs
(samples_*.csv: p_dis(t) and n60(t), 256 samples per cycle) and the mesh pairs MK3_* / M300_* (dx 3 vs dx 1.5 l_ex).

Kept data (steady-tail rule, PROTOKOLL 7.54): a run-amplitude counts only if its last KEEP_LAST = 3 measured cycles
come after the last event cycle (|w_dis / w_loop − 1| > 0.20: a field-triggered state change, 7.53); then these 3
cycles are used. The others are left out and counted.
p_dis comes from the LLG damping term and is ≥ 0. n60 = number of neighbour cell pairs in one particle with an angle
> 60° (vortex cores, Bloch points, steep walls).

Measures per run and amplitude (P = p_dis samples of the kept cycles, Pm = mean cycle, Ps = Fourier part k ≤ 6 of Pm):
  RF    non-smooth share   Σ|P − Ps| / ΣP     all dissipation that the smooth periodic response does not describe
  RN    non-repeating      Σ|P − Pm| / ΣP     differs from cycle to cycle (stochastic jumps)
  RH    sharp repeating    Σ|Pm − Ps| / ΣPm   sharp features at the same phase in each cycle (deterministic jumps)
  dS5   burst excess       share of ∫p dt in the 5 % strongest samples, raw minus the same for Ps (0 = no bursts)
  rho   event coupling     Spearman correlation, within the same phase, between the extra dissipation in a sample
                           interval (P − its mean over the cycles at that phase) and the extra change of n60 in that
                           interval; independent of how often n60 changes; > 0: an extra defect event carries extra loss
  R_np  non-periodic share of the 16-phase snapshots (phasemap.json, 9 / 50 / 150 mT, last 2 cycles)
  B_pk  per run: the amplitude of the largest RF (parabola in ln B through the largest value and its neighbours)

    python eval_jumps_time.py > results/jumps_time.md      # also results/jumps_time.png
"""
import json
import math
import pathlib
import warnings

import numpy as np
from scipy import stats

import analyze
import analyze_design as AD

HERE = pathlib.Path(__file__).resolve().parent
R = HERE / "runs"
KEEP_LAST = 3
K_SMOOTH = 6
N_SAMPLES = 256
N_BOOT = 2000
LEX_A10 = 3.342                       # nm, l_ex at A = 10 pJ/m and Js = 1.5 T
MESH = {"K3 (A 20 pJ/m)": ("MK3_dx%s_s%d_b%d", range(2031, 2037)),
        "base (A 10 pJ/m)": ("M300_dx%s_s%d_b%d", range(921, 925))}
MESH_AMPS = (9, 50, 100)
KEYS = ("RF", "RN", "RH", "dS5", "rho", "ev")


def series(path):
    """p_dis and n60 of the last KEEP_LAST complete cycles (first row of each cycle = end of the previous one)."""
    a = np.loadtxt(path, delimiter=",", comments="#")
    cyc = a[:, 0].astype(int)
    keep = [a[cyc == c][1:] for c in sorted(set(cyc))]
    keep = [x for x in keep if len(x) == N_SAMPLES][-KEEP_LAST:]
    if len(keep) < KEEP_LAST:
        return None
    return np.array([x[:, 8] for x in keep]), np.array([x[:, -1] for x in keep])


def top5(x):
    s = np.sort(x.ravel())[::-1]
    return s[: len(s) // 20].sum() / s.sum()


def measures(P, N60):
    n = len(P)
    Pm = P.mean(0)
    F = np.fft.rfft(Pm)
    F[K_SMOOTH + 1:] = 0
    Ps = np.fft.irfft(F, n=N_SAMPLES)
    tot = P.sum()
    out = {"RF": np.abs(P - Ps).sum() / tot, "RN": np.abs(P - Pm).sum() / tot, "RH": np.abs(Pm - Ps).sum() / Pm.sum(),
           "dS5": top5(P) - top5(np.tile(np.clip(Ps, 0, None), n))}
    pf, nf = P.ravel(), N60.ravel()
    pint = 0.5 * (pf + np.concatenate([pf[:1], pf[:-1]])).reshape(n, -1)
    dn = np.abs(np.diff(nf, prepend=nf[0])).reshape(n, -1)
    out["ev"] = float((dn > 0).mean())
    r, d = (pint - pint.mean(0)).ravel(), (dn - dn.mean(0)).ravel()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out["rho"] = float(stats.spearmanr(r, d)[0]) if d.std() > 0 and r.std() > 0 else float("nan")
    return out


def steady(d, k):
    """Steady-tail rule (7.54): the last KEEP_LAST measured cycles of stage k have no event cycle."""
    st = json.loads((d / "summary.json").read_text())["stages"][k]
    skip = 2 if "h_amp_initial" in st else 1                   # as analyze.load_run (drive corrected after cycle 0)
    return len(analyze.steady_tail(st["cycles"][skip:], AD.BAL_EVENT)) >= KEEP_LAST


def run_measures(d, n_amp):
    files = sorted(d.glob("samples_*_amp*.csv"))
    assert len(files) == n_amp and all(f.name.endswith("_amp%02d.csv" % k) for k, f in enumerate(files)), d
    out = []
    for k, f in enumerate(files):
        s = series(f)
        m = measures(*s) if s is not None else None
        if m is not None:
            m["steady"] = steady(d, k)
        out.append(m)
    return out


def peak_lnB(rf):
    """ln B at the largest RF; parabola in ln B through the largest value and its neighbours; at an end: the end."""
    x = np.log(AD.AMPS_MT)
    i = int(np.nanargmax(rf))
    if i in (0, len(rf) - 1):
        return x[i], True
    if not np.all(np.isfinite(rf[i - 1:i + 2])):
        return x[i], False
    c = np.polyfit(x[i - 1:i + 2], rf[i - 1:i + 2], 2)
    return float(np.clip(-c[1] / (2 * c[0]), x[i - 1], x[i + 1])) if c[0] < 0 else x[i], False


def med_ci(v, rng):
    v = np.asarray([x for x in v if np.isfinite(x)])
    b = [np.median(rng.choice(v, len(v))) for _ in range(N_BOOT)]
    return float(np.median(v)), *np.percentile(b, [2.5, 97.5]), len(v)


def a_pjm(r):
    return 10.0 * ((1000.0 / math.exp(r["lnd"])) / LEX_A10) ** 2


def main():
    import statsmodels.api as sm
    rng = np.random.default_rng(20261010)
    recs, _ = AD.load_runs(R)
    amps = AD.AMPS_MT
    M_all = {r["name"]: run_measures(R / r["name"], len(amps)) for r in recs}
    M = {n: [x if x is not None and x["steady"] else None for x in v] for n, v in M_all.items()}
    rnp = {9: [], 50: [], 150: []}
    for r in recs:
        for st in json.loads((R / r["name"] / "phasemap.json").read_text())["stages"]:
            rnp[min(rnp, key=lambda x: abs(st["B_peak_mT"] / x - 1))].append(st["R_np"])

    print("# Stage 1: jumps in the dissipation time series (exploratory)\n")
    print("n = %d α 0.1 runs (Sobol, centre, corners), the last %d cycles per amplitude if they are steady (7.54), "
          "256 samples per cycle. Medians over the runs, 95 %% CI of the median by bootstrap over the runs (%d "
          "resamples).\n" % (len(recs), KEEP_LAST, N_BOOT))
    print("## A. Measures against the amplitude\n")
    print("| B (mT) | RF non-smooth | RN non-repeating | RH sharp repeating | dS5 burst excess | rho event coupling "
          "(runs > 0) | samples with an n60 change | R_np snapshots | left out, not steady (runs) |")
    print("|---|---|---|---|---|---|---|---|---|")
    curves = {key: [] for key in KEYS}
    for k, m in enumerate(amps):
        v = {key: [M[r["name"]][k][key] for r in recs if M[r["name"]][k] is not None] for key in KEYS}
        for key in KEYS:
            curves[key].append(med_ci(v[key], rng))
        nbad = sum(M[r["name"]][k] is None for r in recs)
        rho = np.array([x for x in v["rho"] if np.isfinite(x)])
        c = lambda key: "%.3f (%.3f … %.3f)" % curves[key][-1][:3]
        print("| %g | %s | %s | %s | %.3f | %.2f (%.2f … %.2f), %d/%d | %.0f %% | %s | %d |" % (
            m, c("RF"), c("RN"), c("RH"), curves["dS5"][-1][0], *curves["rho"][-1][:3], (rho > 0).sum(), len(rho),
            100 * curves["ev"][-1][0], ("%.3f" % np.median(rnp[m])) if m in rnp else "–", nbad))
    print("\nCheck: RF median with the steady rule against all run-amplitudes (the last 3 cycles, also when not steady):\n")
    print("| B (mT) | RF median, steady | RF median, all |")
    print("|---|---|---|")
    for k, m in enumerate(amps):
        a_ = [M_all[r["name"]][k]["RF"] for r in recs if M_all[r["name"]][k] is not None]
        b_ = [M[r["name"]][k]["RF"] for r in recs if M[r["name"]][k] is not None]
        print("| %g | %.3f (n %d) | %.3f (n %d) |" % (m, np.median(b_), len(b_), np.median(a_), len(a_)))

    # B. amplitude of the largest RF per run, against the factors
    pk, ends = [], 0
    for r in recs:
        rf = np.array([x["RF"] if x is not None else np.nan for x in M[r["name"]]])
        lb, end = peak_lnB(rf)
        pk.append(lb)
        ends += end
    pk = np.array(pk)
    X = AD.coded(recs)
    fit = sm.OLS(pk, sm.add_constant(X)).fit(cov_type="HC3")
    print("\n## B. Amplitude B_pk of the largest RF per run\n")
    print("All runs: median B_pk %.0f mT (95 %% CI %.0f … %.0f), n %d. B_pk at 9 or 150 mT (end of the range): %d runs.\n"
          % (*np.exp(med_ci(pk, rng)[:3]), len(pk), ends))
    print("Linear model ln B_pk = c0 + Σ c_j x_j (x_j coded −1 … +1, OLS, HC3 standard errors). Range effect = factor "
          "on B_pk from the low end to the high end of the factor, with 95 % CI:\n")
    print("| factor (low → high end) | factor on B_pk (95 % CI) | p |")
    print("|---|---|---|")
    names = {"lnQ": "L_eff 30 → 12 l_ex (K_eff up)", "r_p": "r_p 0 → 3", "phi": "φ 0.55 → 0.69",
             "lnd": "d/l_ex 212 → 300 (A 20 → 10 pJ/m)"}
    for j, f in enumerate(AD.FACTORS):
        c, se = fit.params[j + 1], fit.bse[j + 1]
        print("| %s | %.2f (%.2f … %.2f) | %.3f |" % (names[f], math.exp(2 * c), math.exp(2 * (c - 1.96 * se)),
                                                     math.exp(2 * (c + 1.96 * se)), fit.pvalues[j + 1]))
    print("\nR² = %.2f. Groups (median B_pk, 95 %% CI):\n" % fit.rsquared)
    print("| group | n | median B_pk (mT) |")
    print("|---|---|---|")
    A = np.array([a_pjm(r) for r in recs])
    q = np.quantile(A, [1 / 3, 2 / 3])
    L = np.array([round(math.sqrt(math.exp(-r["lnQ"]))) for r in recs])
    groups = [("A low third (%.1f … %.1f pJ/m)" % (A.min(), q[0]), A <= q[0]),
              ("A mid third", (A > q[0]) & (A <= q[1])),
              ("A high third (%.1f … %.1f pJ/m)" % (q[1], A.max()), A > q[1])]
    groups += [("L_eff %d l_ex" % v, L == v) for v in sorted(set(L))]
    for name, sel in groups:
        mc = med_ci(pk[sel], rng)
        print("| %s | %d | %.0f (%.0f … %.0f) |" % (name, sel.sum(), *np.exp(mc[:3])))

    # C. mesh pairs
    print("\n## C. Mesh check: the same measures at dx 3 and dx 1.5 l_ex (paired seeds, one job per amplitude)\n")
    print("Jobs with 5 cycles, the last 3 used. A pair counts only if both jobs are steady in these 3 cycles (7.54). "
          "Medians over the seeds. n60 max = largest number of pairs > 60° in the used cycles. Sign test (two-sided, "
          "exact) over the seeds of both points.\n")
    print("| B (mT) | point | n pairs | RF dx 3 | RF dx 1.5 | seeds with RF(dx 3) > RF(dx 1.5) | RN dx 3 | RN dx 1.5 | "
          "n60 max dx 3 | n60 max dx 1.5 | pairs left out (not steady) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    mesh = {}
    for b in MESH_AMPS:
        n_gt, n_all = 0, 0
        for lab, (fmt, seeds) in MESH.items():
            pairs = []
            for s in seeds:
                d3, d15 = R / (fmt % ("3", s, b)), R / (fmt % ("15", s, b))
                if (d3 / "DONE").exists() and (d15 / "DONE").exists():
                    pairs.append([run_measures(d3, 1)[0], run_measures(d15, 1)[0]])
                    for i, d in enumerate((d3, d15)):
                        pairs[-1][i]["n60max"] = int(series(d / "samples_00_amp00.csv")[1].max())
            mesh[(lab, b)] = pairs
            used = [p for p in pairs if p[0]["steady"] and p[1]["steady"]]
            g = lambda i, key: np.array([p[i][key] for p in used])
            gt = int((g(0, "RF") > g(1, "RF")).sum())
            n_gt, n_all = n_gt + gt, n_all + len(used)
            print("| %d | %s | %d | %.3f | %.3f | %d/%d | %.3f | %.3f | %.0f | %.0f | %d |" % (
                b, lab, len(used), np.median(g(0, "RF")), np.median(g(1, "RF")), gt, len(used),
                np.median(g(0, "RN")), np.median(g(1, "RN")), np.median(g(0, "n60max")), np.median(g(1, "n60max")),
                len(pairs) - len(used)))
        print("| %d | both, sign test | %d | | | %d/%d, p = %.3f | | | | | |" % (
            b, n_all, n_gt, n_all, stats.binomtest(n_gt, n_all, 0.5).pvalue))

    fig_out(amps, curves, recs, pk, A, L, mesh)


def fig_out(amps, curves, recs, pk, A, L, mesh):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for f in pathlib.Path("/usr/share/fonts/truetype/crosextra").glob("Carlito*.ttf"):
        font_manager.fontManager.addfont(str(f))
    plt.rcParams.update({"font.family": "Carlito", "font.size": 12, "axes.spines.top": False,
                         "axes.spines.right": False, "legend.frameon": False})
    RED, BLUE, GREY = "#9B111E", "#1F4E79", "#666666"
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
    ax = axs[0]
    for key, col, lab in (("RF", RED, "non-smooth share RF"), ("RN", BLUE, "non-repeating share RN"),
                          ("RH", GREY, "sharp repeating share RH")):
        c = np.array(curves[key])
        ax.plot(amps, c[:, 0], "o-", color=col, lw=2, label=lab)
        ax.fill_between(amps, c[:, 1], c[:, 2], color=col, alpha=0.15)
    ax.set_xscale("log")
    ax.set_xticks(amps)
    ax.set_xticklabels(["%g" % a for a in amps])
    ax.set_xlabel("peak flux density B (mT)")
    ax.set_ylabel("share of the dissipated energy")
    ax.set_title("Jump shares, median of %d runs, 95 %% CI" % len(recs))
    ax.legend(fontsize=10)
    ax = axs[1]
    for v, col in zip(sorted(set(L)), ("#F4B183", "#E06666", RED, "#5B0A12")):
        s = L == v
        ax.scatter(A[s], np.exp(pk[s]), color=col, s=24, alpha=0.8, label="L_eff %d l_ex" % v)
    ax.set_yscale("log")
    ax.set_yticks(amps)
    ax.set_yticklabels(["%g" % a for a in amps])
    ax.set_xlabel("exchange constant A (pJ/m) at d = 1 µm")
    ax.set_ylabel("B_pk: amplitude of the largest RF (mT)")
    ax.set_title("Where the jumps peak, one dot per run")
    ax.legend(fontsize=9, ncol=2)
    ax = axs[2]
    for j, (lab, col) in enumerate(zip(MESH, (RED, BLUE))):
        for i, (dx, mk) in enumerate((("3", "o"), ("1.5", "s"))):
            xs = np.arange(len(MESH_AMPS)) + (2 * j + i - 1.5) * 0.18
            for x, b in zip(xs, MESH_AMPS):
                for p in mesh[(lab, b)]:
                    ax.scatter(x, p[i]["RF"], marker=mk, s=26, color=col if p[i]["steady"] else "none", edgecolors=col,
                               alpha=0.8 if i == 0 else 0.5)
                ax.hlines(np.median([p[i]["RF"] for p in mesh[(lab, b)] if p[0]["steady"] and p[1]["steady"]]), x - 0.08, x + 0.08, color=col, lw=2.5)
            ax.scatter([], [], marker=mk, color=col, label="%s, dx %s l_ex" % (lab.split(" (")[0], dx))
    ax.set_xticks(range(len(MESH_AMPS)))
    ax.set_xticklabels(["%d mT" % b for b in MESH_AMPS])
    ax.set_ylabel("non-smooth share RF")
    ax.set_title("Mesh check: dots = seeds, bars = medians")
    ax.scatter([], [], marker="o", color="none", edgecolors=GREY, label="open: not steady (left out)")
    ax.legend(fontsize=9, loc="upper right")
    fig.tight_layout()
    fig.savefig(HERE / "results" / "jumps_time.png", dpi=200)


if __name__ == "__main__":
    main()
