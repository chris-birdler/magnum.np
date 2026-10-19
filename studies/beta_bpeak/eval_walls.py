#!/usr/bin/env python
"""
Wall measures of the baseline (B300_s921 ... s928, domains.json from domains_batch.py) against
B_peak and against beta (PROTOKOLL 7.19).
  A: per amplitude, mean ± SE over the 8 seeds: texture A_w/V, A_F90/V, F90, l_C, ripple,
     switched volume per half cycle (H max <-> H min) / V, wall path s = V_sw / A_w.
  B: per seed: virgin-state measures, the change of A_w over the amplitude series, beta.
  C: correlation over the seeds (Pearson r, two-sided p; Spearman rho) of the virgin-state
     texture and of the switched-volume exponent with beta and ln P.
One particle of the base point: V_p = (pi/6) d^3; one 180 deg disc through its centre: pi d^2/4.

    python eval_walls.py > results/walls_eval.md
"""
import json
import math

import numpy as np
from scipy import stats

import analyze

RUNS = analyze.pathlib.Path(__file__).resolve().parent / "runs"
SEEDS = range(921, 929)
FK = ["beta_pow", "beta@9-35mT", "beta@35-70mT", "beta@70-150mT", "lnw@10mT", "lnw@50mT", "lnw@100mT"]


def ms(v, fmt="%.3g"):
    v = np.asarray(v, float)
    return (fmt + " ± " + fmt) % (v.mean(), v.std(ddof=1) / math.sqrt(len(v)))


def load(seed):
    name = "B300_s%d" % seed
    d = json.loads((RUNS / name / "domains.json").read_text())
    cfg, rows = analyze.load_run(RUNS / name)
    return d, analyze.features(rows)[0], cfg


def main():
    D = [load(s) for s in SEEDS]
    V = D[0][0]["init"]["V"]
    d_lex = D[0][2]["d_lex_used"]
    disc = math.pi * d_lex ** 2 / 4.0
    print("**A. Wall measures vs B_peak** (baseline, n = 8, mean ± SE over the seeds; per stage the mean of the 4 "
          "phases of the last cycle; V = %.3g l_ex³ for 4 particles; texture also as 180° discs per particle, "
          "disc = π d²/4 = %.3g l_ex²)\n" % (V, disc))
    print("| B_peak | A_w/V [1/l_ex] | discs per particle | A_F90/V [1/l_ex] | F90 | l_C [l_ex] | ripple [°] "
          "| V_sw half cycle / V | wall path s [l_ex] |")
    print("|---|---|---|---|---|---|---|---|---|")
    keys = [("virgin", None)] + [(None, i) for i in range(len(D[0][0]["stages"]))]
    Bs, vsw_all = [], []
    for lab, i in keys:
        A, AF, F, L, R, VS, P = [], [], [], [], [], [], []
        for d, _, _ in D:
            if i is None:
                st = [d["init"]]
                vs, path = np.nan, np.nan
            else:
                s = d["stages"][i]
                st = s["states"]
                vs = s["half_cycle"]["V_sw"][0] / V
                path = s["half_cycle"]["wall_path"][0]
            A.append(np.mean([x["A_w"] for x in st]) / V)
            AF.append(np.mean([x["A_F90"] for x in st]) / V)
            F.append(np.mean([x["F90"] for x in st]))
            L.append(np.mean([x["l_C"] for x in st]))
            R.append(np.mean([x["ripple_rms_deg"] for x in st]))
            VS.append(vs)
            P.append(path)
        if i is None:
            bl = "virgin"
        else:
            bl = "%.0f mT" % np.mean([d["stages"][i]["B_peak_mT"] for d, _, _ in D])
            Bs.append(np.mean([d["stages"][i]["B_peak_mT"] for d, _, _ in D]))
            vsw_all.append(VS)
        print("| %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            bl, ms(A), ms(np.array(A) * V / 4 / disc, "%.2f"), ms(AF), ms(np.array(F) * 100, "%.1f") + " %",
            ms(L, "%.0f"), ms(R, "%.1f"), "–" if i is None else ms(np.array(VS) * 100, "%.2f") + " %",
            "–" if i is None else ms(P, "%.1f")))
    vsw = np.array(vsw_all).T                                    # [seed, stage]
    expo = [np.polyfit(np.log(Bs), np.log(v), 1)[0] for v in vsw]
    print("\nSwitched volume per half cycle vs B_peak: exponent (fit of ln V_sw over ln B, per seed) = %s"
          % ms(expo, "%.2f"))

    print("\n**B. Per seed**\n")
    print("| seed | A_w/V virgin | A_w/V 150 mT − virgin | F90 virgin | l_C virgin | V_sw exponent | β_pow | β 9–35 "
          "| β 35–70 | β 70–150 |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    X = {"texture": [], "dtexture": [], "F90": [], "l_C": [], "V_sw exponent": []}
    Y = {k: [] for k in FK}
    for (d, f, _), s, e in zip(D, SEEDS, expo):
        a0 = d["init"]["A_w"] / V
        a1 = np.mean([x["A_w"] for x in d["stages"][-1]["states"]]) / V
        X["texture"].append(a0)
        X["dtexture"].append(a1 - a0)
        X["F90"].append(d["init"]["F90"])
        X["l_C"].append(d["init"]["l_C"])
        X["V_sw exponent"].append(e)
        for k in FK:
            Y[k].append(f[k])
        print("| %d | %.4f | %+.4f | %.1f %% | %.0f | %.2f | %.2f | %.2f | %.2f | %.2f |" %
              (s, a0, a1 - a0, 100 * d["init"]["F90"], d["init"]["l_C"], e, f["beta_pow"], f["beta@9-35mT"],
               f["beta@35-70mT"], f["beta@70-150mT"]))

    print("\n**C. Correlation over the seeds (n = 8): Pearson r (two-sided p) / Spearman ρ**\n")
    print("| measure | " + " | ".join(FK) + " |")
    print("|---|" + "---|" * len(FK))
    for xk, xv in X.items():
        cells = []
        for k in FK:
            r, p = stats.pearsonr(xv, Y[k])
            rho = stats.spearmanr(xv, Y[k])[0]
            cells.append("%+.2f (p %.2f) / %+.2f" % (r, p, rho))
        print("| %s | %s |" % (xk, " | ".join(cells)))
    print("\n**C2. The same without seed 921** (the run with the events, PROTOKOLL 7.17; n = 7)\n")
    print("| measure | " + " | ".join(FK) + " |")
    print("|---|" + "---|" * len(FK))
    m = np.array([s != 921 for s in SEEDS])
    for xk, xv in X.items():
        cells = []
        for k in FK:
            xa, ya = np.array(xv)[m], np.array(Y[k])[m]
            r, p = stats.pearsonr(xa, ya)
            cells.append("%+.2f (p %.2f) / %+.2f" % (r, p, stats.spearmanr(xa, ya)[0]))
        print("| %s | %s |" % (xk, " | ".join(cells)))
    print("\nWith n = 8 (7), |r| > 0.71 (0.75) is needed for p < 0.05 (two-sided); 5 measures × 7 outputs = 35 tests "
          "per table, so about 2 false hits at p < 0.05 are expected by chance. A large Pearson r with a small "
          "Spearman ρ means that one point carries the correlation.")


if __name__ == "__main__":
    main()
