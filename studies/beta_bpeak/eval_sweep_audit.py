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
    print("\n".join(L))


if __name__ == "__main__":
    main()
