#!/usr/bin/env python
"""
All valid runs for beta in one table (PROTOKOLL 7.16, 7.17): d/l_ex 300, dx 3, base point.
  clean   : converged virgin state, all amplitudes from one run (B300_*, A010c_*)
  spliced : 9 mT from X* (converged), 20 ... 150 mT from F_base_f_* / A010_* (eval_relaxfix)
Stages with w_dis/w_loop outside 0.9 ... 1.1 are listed (kept by PLAN 4.3a).

    python eval_beta_table.py > results/beta_all_valid.md
"""
import json
import math

import numpy as np

import analyze
import eval_relaxfix as E

RUNS = E.RUNS
K = ["power", "beta@10", "beta@50", "beta@100", "P@10", "P@50", "P@100"]


def row_clean(name):
    cfg, rows = analyze.load_run(RUNS / name)
    it = json.loads((RUNS / name / "summary.json").read_text())["init"]["relax_iterations"]
    return cfg["alpha"], rows, "%d" % it


def flags(rows):
    out = ["%.0f mT: %.2f" % (r["b_peak"] * 1500, r["balance"]) for r in rows if abs(r["balance"] - 1) > 0.1]
    return ", ".join(out) or "–"


def main():
    groups = [("α 0.1, clean", [("B300_s%d" % s, "clean") for s in range(921, 929)]),
              ("α 0.1, spliced", [(s, "spliced") for s in E.SEEDS]),
              ("α 0.01, clean", [("A010c_s%d" % s, "clean") for s in (904, 905, 906)]),
              ("α 0.01, spliced", [(s, "spliced") for s in E.SEEDS])]
    print("| group | run | relax it. | stages with w_dis/w_loop outside 0.9 … 1.1 | β 9…150 mT | β 10 mT | β 50 mT "
          "| β 100 mT | P 10 mT | P 50 mT | P 100 mT |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    summary = {}
    for gname, members in groups:
        vals = {k: [] for k in K}
        for m, kind in members:
            if kind == "clean":
                alpha, rows, it = row_clean(m)
                label = m
            else:
                alpha = 0.1 if gname.startswith("α 0.1") else 0.01
                rows, _ = E.spliced(alpha, m)
                label = "%s + %s" % tuple(p % m for p in E.SETS[alpha])
                it = json.loads((RUNS / (E.SETS[alpha][1] % m) / "summary.json").read_text())["init"]
                it = "%d" % it["relax_iterations"]
            f = E.fits(rows)
            for k in K:
                vals[k].append(f[k])
            print("| %s | %s | %s | %s | %.2f | %.2f | %.2f | %.2f | %.1f | %.0f | %.0f |" %
                  (gname, label, it, flags(rows), f["power"], f["beta@10"], f["beta@50"], f["beta@100"],
                   f["P@10"], f["P@50"], f["P@100"]))
        summary[gname] = vals
        n = len(members)
        cells = []
        for k in K:
            v = np.array(vals[k])
            if k.startswith("P"):
                lv = np.log(v)
                cells.append("%.3g ×/ %.2f" % (math.exp(lv.mean()), math.exp(lv.std(ddof=1) / math.sqrt(n))))
            else:
                cells.append("**%.2f ± %.2f** (sd %.2f)" % (v.mean(), v.std(ddof=1) / math.sqrt(n), v.std(ddof=1)))
        print("| **%s: mean ± SE** | n = %d | | | %s |" % (gname, n, " | ".join(cells)))
    for a in ("0.1", "0.01"):
        vals = {k: summary["α %s, clean" % a][k] + summary["α %s, spliced" % a][k] for k in K}
        n = len(vals["power"])
        cells = []
        for k in K:
            v = np.array(vals[k])
            if k.startswith("P"):
                lv = np.log(v)
                cells.append("%.3g ×/ %.2f" % (math.exp(lv.mean()), math.exp(lv.std(ddof=1) / math.sqrt(n))))
            else:
                cells.append("**%.2f ± %.2f** (sd %.2f)" % (v.mean(), v.std(ddof=1) / math.sqrt(n), v.std(ddof=1)))
        print("| **α %s, all: mean ± SE** | n = %d | | | %s |" % (a, n, " | ".join(cells)))
    print()
    print("d/l_ex 300 (1 µm), dx 3, L_eff 12 l_ex, r_p 0, φ 0.65, f = 30 MHz. β: power law 9 … 150 mT and "
          "quadratic fit at 10 / 50 / 100 mT, per run. P: loss density of the unit cell [W/cm³] at the reference "
          "field. ± = SE over the runs (sd = scatter of one run). P mean: geometric mean ×/ exp(SE of ln P).")


if __name__ == "__main__":
    main()
