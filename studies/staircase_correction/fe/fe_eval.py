#!/usr/bin/env python
"""
FE comparison, preparation (PROTOKOLL 3.24): evaluation of fe_states.py (FD energies) and fe_demag.py (FE energies).

Difference method: the FE energy of the coarse m has a large error from the charge sheets on the inner cell faces, but
this error is nearly the same for two bodies on the same mesh. The energy of the smooth spheres is therefore
    E_s = E_vox(FD, exact) + [E_smooth(FE) − E_vox(FE)].
Check of the method on a pair with an exact FD value (the 5-nm staircase, X3 at q 2):
    ε = {[E_q2(FE) − E_vox(FE)] − [E_X3,q2(FD) − E_vox(FD)]} / S,   S = E_vox(FD) − E_X3,q2(FD)
All differences are in units of S (the staircase effect of the state between the 10-nm and the 5-nm staircase).

    python fe_eval.py --results fe_data/fe_night.jsonl --out results/fe_static   (writes .md and .json)
"""
import argparse
import json
import pathlib
from collections import defaultdict

BASE = pathlib.Path(__file__).resolve().parent.parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    runs = [json.loads(l) for l in pathlib.Path(a.results).read_text().splitlines() if l.strip()]
    fe = defaultdict(dict)                          # (data, mesh) -> {(body, m): run}; the last run of a kind wins
    for r in runs:
        if r.get("degree") != 2:
            continue
        key = (r["body"], r["m"])
        old = fe[(r["data"], r["mesh"])].get(key)
        if old is not None:
            print("note: %s %s %s: run with levels %s/%s replaced by levels %s/%s" % (
                r["data"], r["mesh"], key, old.get("levels"), old.get("levels_cut"), r.get("levels"),
                r.get("levels_cut")))
        fe[(r["data"], r["mesh"])][key] = r
    rows, out = [], {}
    for (data, mesh), d in sorted(fe.items()):
        meta = json.loads((BASE / data / "meta.json").read_text())       # data: relative to the study folder
        F = meta["energies_KdVbox"]
        S = F["E_vox"] - F["E_X3_q2"]
        rec = {"state": pathlib.Path(data).name, "mesh": mesh, "S": S, "fd": F}
        if ("vox", "state") in d:
            ev = d[("vox", "state")]["E_KdVbox"]
            rec["fe_vox_rel"] = ev / F["E_vox"] - 1                       # absolute FE error, for information
            for q in (2, 4):
                if ("q%d" % q, "state") in d:
                    rec["eps_q%d" % q] = ((d[("q%d" % q, "state")]["E_KdVbox"] - ev) - (F["E_X3_q%d" % q] - F["E_vox"])) / S
            if ("smooth", "state") in d:
                Es = F["E_vox"] + d[("smooth", "state")]["E_KdVbox"] - ev
                rec["E_s"] = Es
                rec["levels_cut"] = d[("smooth", "state")].get("levels_cut")
                for q in (2, 4):
                    for k in ("X3", "T3Pc"):
                        key = "E_%s_q%d" % (k, q)
                        if key in F:
                            rec["%s_q%d_minus_smooth" % (k, q)] = (F[key] - Es) / S
            if ("smooth", "interp") in d:
                # raw FE value without a difference correction (no exact FD value exists for the continuous m)
                rec["fe_smooth_interp"] = d[("smooth", "interp")]["E_KdVbox"]
                if "E_s" in rec:
                    rec["interp_raw_minus_smooth"] = (rec["fe_smooth_interp"] - rec["E_s"]) / S
        if ("vox", "100") in d and ("smooth", "100") in d:
            # exact check of the difference method on the smooth body: uniform m, target φ (1 − φ) / 3 (φ = 0.55 of the
            # smooth spheres), in units of S of the state
            phi = 0.55
            Es100 = F["E_vox_100"] + d[("smooth", "100")]["E_KdVbox"] - d[("vox", "100")]["E_KdVbox"]
            rec["eps_smooth_uniform"] = (Es100 - phi * (1 - phi) / 3) / S
            rec["fe_raw_smooth_uniform_err"] = (d[("smooth", "100")]["E_KdVbox"] - phi * (1 - phi) / 3) / S
        if ("full", "state") in d:
            rec["fe_full_rel"] = d[("full", "state")]["E_KdVbox"] / F["E_full"] - 1   # FE error without a surface
            rec["fe_full_abs_S"] = (d[("full", "state")]["E_KdVbox"] - F["E_full"]) / S
        rec["cost"] = {"%s/%s" % k: {"cells": r["cells"], "t_total_s": round(r["t_total_s"]), "rss_GB": r["rss_GB"]}
                       for k, r in d.items()}
        out["%s %s" % (rec["state"], mesh)] = rec
        f = lambda key: ("%+.3f" % rec[key]) if key in rec else "–"  # noqa: E731
        rows.append("| %s | %s | %.2e | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            rec["state"], mesh, S, f("eps_q2"), f("eps_q4"), f("eps_smooth_uniform"), f("X3_q2_minus_smooth"),
            f("X3_q4_minus_smooth"), f("T3Pc_q2_minus_smooth"), f("T3Pc_q4_minus_smooth"), f("interp_raw_minus_smooth"),
            f("fe_raw_smooth_uniform_err")))
    md = ["# Static FE comparison: coarse m in the smooth spheres (PROTOKOLL 3.24)\n",
          "FE: magnum.pi 3.13.3, P2 potential, periodic box mesh of n^3 cubes with 6 tetrahedra each, mean field zero. "
          "FD: magnum.np, Newell with full PBC. Values in units of S = E_vox − E_X3,q2 of the state (column S, in "
          "K_d V_box). E_s: energy of the smooth spheres from the difference method.\n",
          "Checks of the difference method: ε_q2, ε_q4 on the 5-nm and 2.5-nm staircases (exact FD values); ε_u on the "
          "smooth spheres with uniform m (exact target φ (1 − φ) / 3). Continuous m: the raw FE value without a "
          "correction; its FE error is not known; the raw FE error of the smooth spheres with uniform m is given "
          "for scale.\n",
          "| state | mesh | S (K_d V_box) | ε_q2 | ε_q4 | ε_u | X3 q2 − E_s | X3 q4 − E_s | T3Pc q2 − E_s "
          "| T3Pc q4 − E_s | continuous m (raw) − E_s | raw FE error, smooth, uniform m |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"] + rows
    md += ["", "Cost per FE run (cells, wall time in s, peak RAM in GB):", ""]
    for key, rec in out.items():
        md.append("- %s: %s" % (key, "; ".join("%s %d cells, %d s, %.1f GB" % (k, v["cells"], v["t_total_s"],
                                                                             v["rss_GB"]) for k, v in rec["cost"].items())))
    pathlib.Path(a.out + ".md").write_text("\n".join(md) + "\n")
    pathlib.Path(a.out + ".json").write_text(json.dumps(out, indent=1))
    print("\n".join(md))


if __name__ == "__main__":
    main()
