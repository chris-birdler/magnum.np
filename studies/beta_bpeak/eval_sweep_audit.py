#!/usr/bin/env python
"""
Frequency data at 9 mT (C′ and sweep, PROTOKOLL 7.91 … 7.95): what was computed, and how accurate each value is
(question of Chris 2026-10-14). Per run (corner K1, A 20 / A 10, seeds 0 … 2, 8 frequencies): source, samples per
cycle, cycles run, kept cycles (accumulator window), reached B_peak (mean of the kept cycles), w = loss per cycle (mean
of the kept cycles, raw, at the reached B_peak), cycle-to-cycle SD of w, energy balance ε = w_dis / w_loop − 1, and
the model parameters that differ between A 20 and A 10.

    python eval_sweep_audit.py > results/sweep9_audit.md
"""
import json

import numpy as np

import eval_sweep_state as ES

KEYS = ("d_lex_used", "d", "N_used", "Leff_lex", "dx_lex", "alpha", "Js", "r_p", "phi_vox")


def main():
    L = ["# Frequency data at 9 mT: runs, values and accuracy (corner K1)", "",
         "w in 10⁻⁷ (w_loop / K_d, raw at the reached B_peak). Cycle SD: SD of w over the kept cycles / mean.", "",
         "| A | f (MHz) | seed | run | samples / cycle | cycles run | kept cycles | B_peak (mT) | w | cycle SD (%) | ε |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    bp, sd, gm = [], [], {}
    for a in (20, 10):
        for f in ES.F_MHZ:
            ws = []
            for k in ES.SEEDS:
                d = ES.folder(a, f, k)
                cfg = json.loads((d / "config.json").read_text())
                st = json.loads((d / "summary.json").read_text())["stages"][0]
                kept = json.loads((d / "acc.json").read_text())["stages"][0]["harm"]["cycles_kept"]
                cyc = {c["cycle"]: c for c in st["cycles"]}
                w = np.array([cyc[c]["w_loop"] for c in kept])
                b = 1000 * cfg["Js"] * np.mean([cyc[c]["b_peak"] for c in kept])
                s = 100 * w.std(ddof=1) / w.mean() if len(w) > 1 else float("nan")
                e = np.mean([cyc[c]["w_dis"] / cyc[c]["w_loop"] - 1 for c in kept])
                bp.append(b)
                sd.append(s)
                ws.append(w.mean())
                L.append("| %d | %d | %d | %s | %d | %d | %s | %.2f | %.2f | %.1f | %+.3f |" % (
                    a, f, k, d.name, cfg["samples"], len(st["cycles"]), kept, b, 1e7 * w.mean(), s, e))
            gm[(a, f)] = ws
    L += ["", "Reached B_peak: %.2f … %.2f mT (w ∝ B^2.5 gives at most %.0f %% error in w). Cycle SD: median %.1f %%, "
          "largest %.1f %%." % (min(bp), max(bp), 100 * ((max(abs(np.array(bp) / 9 - 1))) * 2.5),
                                 float(np.nanmedian(sd)), float(np.nanmax(sd))), ""]
    a10 = json.loads((ES.folder(10, 30, 0) / "config.json").read_text())
    a20 = json.loads((ES.folder(20, 30, 0) / "config.json").read_text())
    L += ["Model parameters (seed 0; the run code works in units of l_ex at A 10 pJ/m):", ""]
    for kk in KEYS + ("anisotropy",):
        v10, v20 = a10.get(kk), a20.get(kk)
        if kk == "anisotropy":
            for sk in ("K_eff", "L_eff", "K_p", "n_cubes_per_particle"):
                L.append("- anisotropy.%s: A 10 run %s, A 20 run %s" % (sk, v10[sk], v20[sk]))
        else:
            L.append("- %s: A 10 run %s, A 20 run %s" % (kk, v10, v20))
    # loss density P = f K_d w (W/cm³, core volume) at 9 mT, steady rule as in eval_sweep (w9: raw at the reached
    # B_peak); mean ± SD between the seeds (independent microstructures), SE = SD / √n, n seeds with a steady value
    import math
    import eval_sweep as SW
    kd = 1.5 ** 2 / (2 * 4e-7 * math.pi) / 1e6                          # K_d in J/cm³
    L += ["", "## Loss density P at 9 mT (W/cm³), steady rule as in eval_sweep", "",
          "| particle | f (MHz) | n | mean | SD between seeds | SE | per seed |", "|---|---|---|---|---|---|---|"]
    for a, lab in ((10, "large (A 10, d/l_ex 304)"), (20, "small (A 20, d/l_ex 216)")):
        for f in SW.F_MHZ:
            P = np.array([f * 1e6 * kd * SW.w9(SW.folder(a, f, k)) for k in range(3)])
            v = P[np.isfinite(P)]
            sdv = v.std(ddof=1) if len(v) > 1 else float("nan")
            L.append("| %s | %.0f | %d | %.2f | %.2f (%.0f %%) | %.2f | %s |" % (
                lab, f, len(v), v.mean(), sdv, 100 * sdv / v.mean(), sdv / math.sqrt(len(v)),
                " / ".join("%.2f" % x if np.isfinite(x) else "no steady value" for x in P)))
    # the shift by the factor 2 in f at all B (C′ only: 15 / 30 / 60 MHz, ln w at the nominal B of eval_freq). If
    # w20(2 f) = c · w10(f), the two shifted ratios w20(30)/w10(15) and w20(60)/w10(30) are equal.
    import analyze
    import analyze_design as AD
    import eval_freq as EF
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    W = EF.lnw_all()
    L += ["", "## Shift by the factor 2 in f at all B (C′, geometric mean of 3 seeds)", "",
          "| B (mT) | w20/w10 at 15 / 30 / 60 MHz | w20(30)/w10(15) | w20(60)/w10(30) |", "|---|---|---|---|"]
    for i, b in enumerate(EF.AMPS):
        g = lambda a, fr: np.array([W[(a, fr, k)][i] for k in range(3)])
        un = [np.exp(np.mean(g(20, fr) - g(10, fr))) for fr in ("fhalf", "f", "f2")]
        L.append("| %g | %s | %.2f | %.2f |" % (b, " / ".join("%.2f" % x for x in un),
                                               np.exp(np.mean(g(20, "f") - g(10, "fhalf"))),
                                               np.exp(np.mean(g(20, "f2") - g(10, "f")))))
    print("\n".join(L))


if __name__ == "__main__":
    main()
