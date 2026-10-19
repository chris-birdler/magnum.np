#!/usr/bin/env python
"""
Harmonics of the smooth (reversible) response (PROTOKOLL 7.22): which part of the reversible
response becomes less lossy at large amplitude?

Per kept cycle (from cycle 2 on), clean runs (B300_*, A010c_*), from samples_*.csv:
  H1, M1     fundamental of H(t) and M_par(t) (complex), loss angle delta = phase lag of M1 behind H1;
             for the fundamental w = pi mu0 |H1| |M1| sin(delta); with mu = const,
             beta = 2 + d ln sin(delta) / d ln B
  M3/M1      third harmonic of M_par (non-linearity of the reversible M(H))
  p2/p0      second harmonic / mean of the dissipation power (a linear viscous response with a
             sinusoidal dM/dt gives p ~ cos^2: p2/p0 = 1)
  eta        p0 / <(dM_par/dt)^2>: dissipation per unit of the macroscopic magnetization rate.
             eta = const for a linear viscous response; a falling eta means that internal motion of
             m (not seen in <M>: ripple, vortex cores, precession) carries a part of the dissipation
             that grows slower than the rate of <M>.
Mean ± SE over the runs (each run: mean over its kept cycles).

    python eval_harmonics.py > results/harmonics.md
"""
import math
import pathlib

import numpy as np

import eval_jumps as J

AMPS = J.AMPS


def harm(x, k):
    n = len(x)
    ph = 2 * math.pi * np.arange(n) / n
    return 2.0 / n * np.sum(x * np.exp(-1j * k * ph))


def cycle(a):
    n = len(a) - 1
    t, H, M, p = a[:n, 1], a[:n, 2], a[:n, 3], a[:n, 8]
    T = (a[n, 1] - a[0, 1])
    w = 2 * math.pi / T
    H1, M1, M3 = harm(H, 1), harm(M, 1), harm(M, 3)
    delta = np.angle(H1) - np.angle(M1)                        # lag of M behind H
    delta = (delta + math.pi) % (2 * math.pi) - math.pi
    p0 = p.mean()
    p2 = abs(harm(p, 2))
    rate2 = 0.5 * sum((k * w) ** 2 * abs(harm(M, k)) ** 2 for k in range(1, 16))   # <(dM/dt)^2>
    return {"sin_delta": math.sin(delta), "M3/M1": abs(M3) / abs(M1), "p2/p0": p2 / p0, "eta": p0 / rate2,
            "mu_r": 1.0 + abs(M1) / abs(H1)}


def main():
    keys = ["sin_delta", "M3/M1", "p2/p0", "eta", "mu_r"]
    data = {}
    for alpha, runs in J.GROUPS.items():
        for run in runs:
            for f in J.stage_files(run):
                a = np.loadtxt(f, delimiter=",", comments="#")
                cyc = a[:, 0].astype(int)
                per = []
                Bk = 0.0
                for k in sorted(set(cyc))[J.SKIP:]:
                    sel = a[cyc == k]
                    per.append(cycle(sel))
                    Bk = 1000.0 * np.abs(sel[:, 7]).max()
                mT = min(AMPS, key=lambda x: abs(x - Bk))
                data.setdefault((alpha, mT), []).append({q: np.mean([c[q] for c in per]) for q in keys})
    print("| α | B_peak | n | sin δ (fundamental) | μ_r of the core (fundamental) | M3/M1 | p2/p0 | η / η(9 mT) |")
    print("|---|---|---|---|---|---|---|---|")
    for alpha in J.GROUPS:
        eta0 = None
        for mT in AMPS:
            d = data.get((alpha, mT))
            if not d:
                continue
            v = {q: np.array([x[q] for x in d]) for q in keys}
            eta0 = eta0 or v["eta"].mean()
            se = lambda x: x.std(ddof=1) / math.sqrt(len(x))
            print("| %g | %d mT | %d | %.3f ± %.3f | %.2f ± %.2f | %.3f ± %.3f | %.2f ± %.2f | %.2f ± %.2f |" %
                  (alpha, mT, len(d), v["sin_delta"].mean(), se(v["sin_delta"]), v["mu_r"].mean(), se(v["mu_r"]),
                   v["M3/M1"].mean(), se(v["M3/M1"]), v["p2/p0"].mean(), se(v["p2/p0"]),
                   v["eta"].mean() / eta0, se(v["eta"]) / eta0))
    print("\nOver the range: β = 2 + d ln sin δ / d ln B (for μ = const); the fit of ln sin δ over ln B per group:")
    for alpha in J.GROUPS:
        B = [mT for mT in AMPS if (alpha, mT) in data]
        sd = [np.mean([x["sin_delta"] for x in data[(alpha, mT)]]) for mT in B]
        s = np.polyfit(np.log(B), np.log(sd), 1)[0]
        print("  α %g: d ln sin δ / d ln B = %.2f → β ≈ %.2f (if μ = const)" % (alpha, s, 2 + s))
    print("\n**Linear cross-check** (not a mechanism): for a linear, viscously damped response below resonance "
          "tan δ ≈ α ω / ω_eff, thus f_eff = α f / tan δ, and tan δ should scale with α. Ratio tan δ(α 0.01) / "
          "tan δ(α 0.1) = 0.1 if linear.")
    print("| B_peak | f_eff α 0.1 [MHz] | f_eff α 0.01 [MHz] | tan δ(0.01) / tan δ(0.1) |")
    print("|---|---|---|---|")
    F = 30.0
    for mT in AMPS:
        if (0.1, mT) not in data or (0.01, mT) not in data:
            continue
        t1 = np.mean([math.tan(math.asin(x["sin_delta"])) for x in data[(0.1, mT)]])
        t2 = np.mean([math.tan(math.asin(x["sin_delta"])) for x in data[(0.01, mT)]])
        print("| %d mT | %.0f | %.0f | %.2f |" % (mT, 0.1 * F / t1, 0.01 * F / t2, t2 / t1))
    print("\nα 0.1: B300_s921 … s928; α 0.01: A010c_s904 … s906; cycles from 2 on. μ_r of the core = 1 + |M1|/|H1| "
          "(M averaged over the unit cell).")


if __name__ == "__main__":
    main()
