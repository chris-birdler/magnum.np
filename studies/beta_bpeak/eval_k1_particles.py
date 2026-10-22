#!/usr/bin/env python
"""
K1 investigation (PROTOKOLL 7.137): the loss of each particle from the accumulators (acc_reduce.py, key "particles":
the Gilbert dissipation of the particle cells) against the work of the applied field on the particle
(eval_k1_loops.py: f * mu0 * closed integral of H dM_p, PROTOKOLL 7.136). The two are equal only in the sum over
the particles: the stray field moves energy from one particle to another.
  D_p   dissipation density of particle p, W/cm³ of the particle volume: P * share_p * 4 / phi_vox (P: the box value,
        f * W_loop of the original run, mean of the last 2 cycles of the stage)
  F_p   work of the applied field on particle p, W/cm³ of the particle volume (only at the stages of eval_k1_loops)
  c     concentration of the particle loss in its core (r < 0.5) and in the texture top 10 % (acc_reduce.py)
  core  mean m.e in the core of the particle in the prepared state (eval_k1_states.py, key core_mx_mean)

    python eval_k1_particles.py --acc DIR --states k1_states.json [--runs DIR] [--json FILE] RUN ...
DIR holds <RUN>/acc.json (acc_reduce.py --out DIR, or a replay folder) for each RUN (for example K1_4_rep50); the
original run folder is <RUN without the suffix>_sat.
"""
import argparse
import json
import pathlib

import numpy as np

import eval_k1_loops as EL

HERE = pathlib.Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--acc", required=True)
    ap.add_argument("--states", required=True)
    ap.add_argument("--runs", default=str(HERE / "runs"))
    ap.add_argument("--json", default=None)
    ap.add_argument("names", nargs="+")
    a = ap.parse_args()
    runs, states = pathlib.Path(a.runs), json.loads(pathlib.Path(a.states).read_text())
    out = {}
    print("# Loss of each particle: dissipation (accumulators) against the work of the applied field\n")
    print("D_p, F_p in W/cm³ of the particle volume; c in the particle; core: mean m.e in the core of the prepared "
          "state.\n")
    for name in a.names:
        orig = name.rsplit("_", 1)[0] + "_sat" if not name.endswith("_sat") else name
        cfg = json.loads((runs / orig / "config.json").read_text())
        summ = json.loads((runs / orig / "summary.json").read_text())
        acc = json.loads((pathlib.Path(a.acc) / name / "acc.json").read_text())
        core = [p["core_mx_mean"] for p in states[orig]["init.pt"]] if orig in states else [float("nan")] * 4
        out[name] = {"original": orig, "core_m_e": core, "stages": {}}
        for st in acc["stages"]:
            i = int(st["stage"][3:])
            P = float(np.mean([c["W_loop"] for c in summ["stages"][i]["cycles"][-2:]])) * cfg["freq"] / 1e6
            D = [P * q["p_share"] * 4.0 / cfg["phi_vox"] for q in st["particles"]]
            F = None
            if (runs / orig / ("samples_%02d_%s.csv" % (i, st["stage"]))).exists() and i in (0, 3):
                F = [cfg["freq"] / 1e6 * w for w in EL.run_record(runs / orig, "%02d_%s" % (i, st["stage"]))["w_p"]]
            out[name]["stages"][st["stage"]] = {
                "B_mT": st["B_peak_mT"], "P": P, "D_p": D, "F_p": F,
                "R1_p": [q["R1"] for q in st["particles"]],
                "c_core_p": [q["c r<0.5"] for q in st["particles"]],
                "c_texture_p": [q["c texture top 10%"] for q in st["particles"]]}
    for name, r in out.items():
        print("## %s (core m.e: %s)\n" % (name, " / ".join("%.3f" % v for v in r["core_m_e"])))
        print("| stage | B_peak, mT | P | D_0 / D_1 / D_2 / D_3 | F_0 / F_1 / F_2 / F_3 | R1 per particle | "
              "c core | c texture top 10 % |")
        print("|---|---|---|---|---|---|---|---|")
        for stg, s in r["stages"].items():
            f = " / ".join("%.0f" % v for v in s["F_p"]) if s["F_p"] else "-"
            print("| %s | %.1f | %.1f | %s | %s | %s | %s | %s |" % (
                stg, s["B_mT"], s["P"], " / ".join("%.0f" % v for v in s["D_p"]), f,
                " / ".join("%.2f" % v for v in s["R1_p"]), " / ".join("%.2f" % v for v in s["c_core_p"]),
                " / ".join("%.2f" % v for v in s["c_texture_p"])))
        print()
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
