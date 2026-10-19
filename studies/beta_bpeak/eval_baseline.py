#!/usr/bin/env python
"""
Evaluation of the repeats with the converged virgin state (PLAN 3e, PROTOKOLL 7.17):
  1. baseline (B300_s921 ... s928): d/l_ex 300, dx 3, base point, alpha 0.1, 7 amplitudes;
     quality per run, realisation scatter of beta and of the loss density;
  2. alpha 0.01 vs 0.1: clean runs (A010c_* vs B300_*) and pooled with the spliced runs
     (eval_relaxfix: 9 mT from X*, 20 ... 150 mT from F_base_f_* / A010_*);
  3. number of Sobol points for SE <= 0.1 of the range effects on beta (PLAN 4.3).

    python eval_baseline.py > results/baseline_eval.txt
"""
import json
import math
import pathlib

import numpy as np

import analyze
import eval_beta_segments as G
import eval_beta_windows as W
import eval_relaxfix as E

RUNS = pathlib.Path(__file__).resolve().parent / "runs"
BASE = ["B300_s%d" % s for s in range(921, 929)]
A010C = ["A010c_s%d" % s for s in (904, 905, 906)]
KEYS = ["beta_pow", "beta@9-35mT", "beta@35-70mT", "beta@70-150mT"]     # analyze.features (PLAN 4.3)
KD, FREQ = 895e3, 30e6
SE_PER_SIGMA = 0.31        # PLAN 4.3: review Monte Carlo, N = 128 + 8: SE = 0.11 ... 0.14 at sigma 0.45
N_REF = 136


def mean_se(v):
    v = np.asarray(v, float)
    return v.mean(), v.std(ddof=1), v.std(ddof=1) / math.sqrt(len(v)), len(v)


def welch(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = a.mean() - b.mean()
    se = math.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
    return d, se


def clean(name):
    cfg, rows = analyze.load_run(RUNS / name)                 # gate: converged virgin state
    init = json.loads((RUNS / name / "summary.json").read_text())["init"]
    return cfg, rows, init


def main():
    print("Repeats with the converged virgin state. d/l_ex 300 (1 um), dx 3, L_eff 12 l_ex, r_p 0, phi 0.65,")
    print("f = 30 MHz. beta: analyze.features: power law over 3 neighbouring amplitudes (9-20-35, 35-50-70,")
    print("70-100-150 mT) and over all amplitudes (beta_pow); P from ln w interpolated at 10 / 50 / 100 mT.")
    print("Errors: sd = realisation scatter (one run), SE = sd / sqrt(n).\n")

    print("1. baseline, alpha 0.1, n = %d" % len(BASE))
    print("%-10s %6s %5s %-38s %6s %7s %6s | %5s %5s %5s %5s" %
          ("run", "relax", "iter", "w_dis/w_loop per amplitude", "closure", "offset", "b hit",
           "pow", "9-35", "35-70", "70-150"))
    F, lnP = {k: [] for k in KEYS}, {10: [], 50: [], 100: []}
    for name in BASE:
        cfg, rows, init = clean(name)
        f = analyze.features(rows)[0]
        tgt = sorted(cfg["b_list"])
        hit = max(abs(r["b_peak"] / t - 1.0) for r, t in zip(rows, tgt))
        print("%-10s %6s %5d %-38s %6.4f %7.3f %5.1f%% | %5.2f %5.2f %5.2f %5.2f" %
              (name, init["relax_converged"], init["relax_iterations"],
               " ".join("%.2f" % r["balance"] for r in rows), max(r["closure"] for r in rows),
               max(abs(r["offset"]) for r in rows), 100 * hit, *[f[k] for k in KEYS]))
        for k in KEYS:
            F[k].append(f[k])
        for mT in lnP:
            lnP[mT].append(f["lnw@%gmT" % mT] + math.log(KD * FREQ / 1e6))
    print("\n   %-10s %8s %8s %8s" % ("", "mean", "sd", "SE"))
    for k in KEYS:
        m, sd, se, n = mean_se(F[k])
        print("   %-10s %8.2f %8.2f %8.2f" % (k, m, sd, se))
    for mT, v in lnP.items():
        m, sd, se, n = mean_se(v)
        print("   P@%-3d mT  %8.3g W/cm^3  scatter factor x/ %.2f (sd of ln P %.2f), mean x/ %.2f"
              % (mT, math.exp(m), math.exp(sd), sd, math.exp(se)))

    print("\n2. alpha 0.01 - 0.1 (Welch): 3-point power-law fits on the common amplitudes 9/20/50/100/150 mT")
    for name in A010C:
        cfg, rows, init = clean(name)
        print("   %s: relax %s %d it, balance %s" %
              (name, init["relax_converged"], init["relax_iterations"], " ".join("%.2f" % r["balance"] for r in rows)))
    labels = ["9-20-50", "20-50-100", "50-100-150"]
    res = {}
    for alpha in (0.1, 0.01):
        runs = G.group(alpha)
        res[alpha] = (np.array([[x[0] for x in W.window_fits(W.pick(r, G.COMMON))] for r in runs]),
                      np.array([analyze.features(r)[0]["beta_pow"] for r in runs]),
                      [n.startswith("B300") for n in ["B300"] * 8 + ["s"] * 3] if alpha == 0.1 else
                      [True] * 3 + [False] * 3)
    for i, lab in enumerate(labels + ["all (beta_pow)"]):
        a01 = res[0.01][0][:, i] if i < 3 else res[0.01][1]
        a1 = res[0.1][0][:, i] if i < 3 else res[0.1][1]
        ca, c1 = np.array(res[0.01][2]), np.array(res[0.1][2])
        d1, s1 = welch(a01[ca], a1[c1])
        d2, s2 = welch(a01, a1)
        print("   %-16s clean (3 vs 8) %+.2f ± %.2f (%.1f SE)   all (6 vs 11) %+.2f ± %.2f (%.1f SE)"
              % (lab, d1, s1, abs(d1) / s1, d2, s2, abs(d2) / s2))

    print("\n3. Sobol points for SE <= 0.1 of the range effect on beta (PLAN 4.3)")
    print("   Scaling of the review Monte Carlo: SE = %.2f sigma sqrt(%d / N) (worst factor)." % (SE_PER_SIGMA, N_REF))
    for k in KEYS:
        sd = mean_se(F[k])[1]
        N = N_REF * (SE_PER_SIGMA * sd / 0.1) ** 2
        print("   %-10s sigma %.2f  ->  N = %4.0f  (N = 64: SE %.3f, N = 128: SE %.3f)" %
              (k, sd, N, SE_PER_SIGMA * sd * math.sqrt(N_REF / 64), SE_PER_SIGMA * sd * math.sqrt(N_REF / 128)))
    print("   Note: sigma of one point of the design also contains the parameter dependence that the")
    print("   15-term regression does not describe; the 8 centre replicates check this (lack of fit).")


if __name__ == "__main__":
    main()
