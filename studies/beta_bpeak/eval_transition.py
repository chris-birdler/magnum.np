#!/usr/bin/env python
"""
What changes near 30 … 70 mT? Exploratory (PROTOKOLL 7.50), α 0.1 Sobol runs. Per amplitude: median share R1 of the loss
in the local fundamental (accumulators), loss shares core / p top 10 %, jump flags, energy-balance events, local β
between neighbour amplitudes; the amplitude where R1 falls below 0.6, for thirds of each parameter (and the L_eff levels
with the anisotropy field H_K = 2 K_eff / (μ0 Ms)); the drive field H = B / (μ0 μ_r) with μ_r = 8.

    python eval_transition.py > results/transition.md
"""
import json
import math
import pathlib

import numpy as np

import analyze
import analyze_design as AD

R = pathlib.Path(__file__).resolve().parent / "runs"
MU0, MS, MU_R = 4e-7 * math.pi, 1.5 / (4e-7 * math.pi), 8.0


def r1(r, m):
    return min(json.loads((R / r["name"] / "acc.json").read_text())["stages"], key=lambda s: abs(s["B_peak_mT"] / m - 1))


def b_at(R1, amps, level=0.6):
    return next((amps[i] * (amps[i + 1] / amps[i]) ** ((R1[i] - level) / (R1[i] - R1[i + 1])) for i in range(len(amps) - 1)
                 if R1[i] >= level > R1[i + 1]), float("nan"))


def main():
    recs, _ = AD.load_runs(R)
    amps = AD.AMPS_MT
    n = len(recs)
    print("# What changes near 30 … 70 mT (exploratory, n = %d)\n" % n)
    print("| B (mT) | H (kA/m) | R1 median | core share | p top 10 % share | jump flags | balance > 20 % |")
    print("|---|---|---|---|---|---|---|")
    for m in amps:
        st = [r1(r, m) for r in recs]
        rows = [min(analyze.load_run(R / r["name"])[1], key=lambda q: abs(1500 * q["b_peak"] / m - 1)) for r in recs]
        print("| %g | %.1f | %.2f | %.2f | %.2f | %d | %d |" % (
            m, m * 1e-3 / (MU0 * MU_R) / 1e3, np.median([s["R1"] for s in st]), np.median([s["where"]["r<0.5"]["p_share"] for s in st]),
            np.median([s["where"]["p top 10%"]["p_share"] for s in st]), sum("#" in analyze.flags(q) for q in rows),
            sum(abs(q["balance"] - 1) > 0.2 for q in rows)))
    L = np.array([[r["lnw@%gmT" % m] for m in amps] for r in recs])
    lb = np.diff(L, axis=1) / np.diff(np.log(amps))
    print("\nLocal β between neighbour amplitudes (mean ± SE): " + ", ".join(
        "%g–%g mT %.2f ± %.2f" % (amps[i], amps[i + 1], lb[:, i].mean(), lb[:, i].std(ddof=1) / math.sqrt(n)) for i in range(6)))
    print("\n## Amplitude where R1 falls below 0.6\n")
    print("| group | n | R1 at 9 … 150 mT | B (mT) at R1 = 0.6 |")
    print("|---|---|---|---|")
    kd = json.loads((R / recs[0]["name"] / "config.json").read_text())["units"]["Kd"]
    groups = {}
    for r in recs:
        groups.setdefault(round(math.sqrt(math.exp(-r["lnQ"]))), []).append(r)
    for Lv in sorted(groups):
        K = kd / Lv ** 2
        R1 = [np.median([r1(r, m)["R1"] for r in groups[Lv]]) for m in amps]
        print("| L_eff %d l_ex (K_eff %.2f kJ/m³, H_K %.1f kA/m) | %d | %s | %.0f |" % (
            Lv, K / 1e3, 2 * K / (MU0 * MS) / 1e3, len(groups[Lv]), " ".join("%.2f" % x for x in R1), b_at(R1, amps)))
    for lab, f in (("A at 1 µm (pJ/m)", lambda r: 10 * ((1000 / math.exp(r["lnd"])) / 3.342) ** 2),
                   ("r_p", lambda r: r["r_p"]), ("packing fraction", lambda r: r["phi"])):
        v = np.array([f(r) for r in recs])
        q = np.quantile(v, [1 / 3, 2 / 3])
        for name, sel in (("low third", v <= q[0]), ("high third", v > q[1])):
            rr = [r for r, s in zip(recs, sel) if s]
            R1 = [np.median([r1(r, m)["R1"] for r in rr]) for m in amps]
            print("| %s, %s (%.2f … %.2f) | %d | %s | %.0f |" % (lab, name, v[sel].min(), v[sel].max(), len(rr),
                                                                " ".join("%.2f" % x for x in R1), b_at(R1, amps)))


if __name__ == "__main__":
    main()
