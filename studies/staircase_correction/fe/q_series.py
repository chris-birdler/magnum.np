#!/usr/bin/env python
"""
Series in q with FD (PROTOKOLL 3.25, 3.27): staircase effect of the exact fine staircase (X3) and of the correction
(T3Pc) for each state, from the meta.json files of fe_states.py.

    staircase effect at q = (E_vox − E_q) / S,   S = E_vox − E_X3,q2   (E in K_d V_box, same coarse m)

Estimate of the smooth spheres (an assumption: convergence of first order in 1/q, so each doubling of q halves the
increment): X3 at q 8 plus the increment from q 4 to q 8.

    python q_series.py --data fe_data --out results/q_series   (writes .md and .json)
"""
import argparse
import json
import pathlib

QS = (2, 4, 8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="folder with one sub-folder (meta.json) per state")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out, rows = {}, []
    for d in sorted(pathlib.Path(a.data).iterdir()):
        if not (d / "meta.json").exists():
            continue
        F = json.loads((d / "meta.json").read_text())["energies_KdVbox"]
        if "E_X3_q2" not in F:
            continue
        S = F["E_vox"] - F["E_X3_q2"]
        rec = {"S": S}
        for k in ("X3", "T3Pc"):
            for q in QS:
                if "E_%s_q%d" % (k, q) in F:
                    rec["%s_q%d" % (k, q)] = (F["E_vox"] - F["E_%s_q%d" % (k, q)]) / S
        for q in QS:
            if "X3_q%d" % q in rec and "T3Pc_q%d" % q in rec:
                rec["T3Pc_minus_X3_q%d" % q] = rec["X3_q%d" % q] - rec["T3Pc_q%d" % q]   # energy difference / S
        if all("X3_q%d" % q in rec for q in QS):
            rec["increment_ratio"] = (rec["X3_q8"] - rec["X3_q4"]) / (rec["X3_q4"] - rec["X3_q2"])
            rec["smooth_estimate"] = 2 * rec["X3_q8"] - rec["X3_q4"]
        out[d.name] = rec
        f = lambda key, fmt="%.3f": (fmt % rec[key]) if key in rec else "–"  # noqa: E731
        rows.append("| %s | %.4e | %s / %s / %s | %s / %s / %s | %s / %s / %s | %s | %s |" % (
            d.name, S, f("X3_q2"), f("X3_q4"), f("X3_q8"), f("T3Pc_q2"), f("T3Pc_q4"), f("T3Pc_q8"),
            *(("%+.2f" % (100 * rec["T3Pc_minus_X3_q%d" % q])) if "T3Pc_minus_X3_q%d" % q in rec else "–" for q in QS),
            f("increment_ratio"), f("smooth_estimate")))
    md = ["# Series in q with FD: staircase effect of the q-staircase at the same coarse m (PROTOKOLL 3.27)\n",
          "Staircase effect (E_vox − E_q) / S, with S = E_vox − E_X3,q2 of the state (column S, in K_d V_box). X3: exact "
          "field of the q-staircase (fine FFT on (48 q)³). T3Pc: Newell with M_s v plus the pairwise correction with the "
          "centroid term (p 3). T3Pc − X3: energy difference in % of S. Increment ratio: (X3 q 8 − X3 q 4) / (X3 q 4 − "
          "X3 q 2), 0.5 for convergence of first order in 1/q. Smooth spheres: estimate 2 X3 q 8 − X3 q 4 (assumption: "
          "first order in 1/q).\n",
          "| state | S (K_d V_box) | X3 q 2 / 4 / 8 | T3Pc q 2 / 4 / 8 | T3Pc − X3, q 2 / 4 / 8 (% of S) "
          "| increment ratio | smooth spheres (estimate) |",
          "|---|---|---|---|---|---|---|"] + rows
    pathlib.Path(a.out + ".md").write_text("\n".join(md) + "\n")
    pathlib.Path(a.out + ".json").write_text(json.dumps(out, indent=1))
    print("\n".join(md))


if __name__ == "__main__":
    main()
