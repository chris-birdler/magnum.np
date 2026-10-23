#!/usr/bin/env python
"""
Criteria of stage A (PLAN.md section 3) from the files static_A_<geometry>[_q<q>].json of check_static.py, and the
smallest eigenvalues of check_sign.py (R3) if present.

    python eval_static.py results/static_A_*.json [--sign results/sign_*.json] > results/static_A.md
"""
import argparse
import glob
import json
import math

TOL_R1 = 5e-4


def crit(d, form, p, EK):
    """Criteria R2 of one state, form and window; returns (dict of numbers, list of used criteria, pass)."""
    E, Eref, S = d["E_%s%d" % (form, p)], d["E_ref"], d["S"]
    out = {"e_E": abs(E - Eref) / abs(S) if S else float("inf"), "dE_EK": abs(E - Eref) / abs(EK),
           "e_h": d["rms_layer_%s%d" % (form, p)] / d["rms_layer_vox"] if d["rms_layer_vox"] else float("inf"),
           "rms_h": d["rms_layer_%s%d" % (form, p)]}
    used, ok = [], True
    if abs(S) >= 0.01 * abs(EK):
        used.append("e_E")
        ok &= out["e_E"] <= 0.10
    used.append("E_K")
    ok &= out["dE_EK"] <= 0.10
    if d["rms_layer_vox"] >= 1e-3:
        used.append("e_h")
        ok &= out["e_h"] <= 0.10
    else:
        used.append("rms_h")
        ok &= out["rms_h"] <= 1e-4
    return out, used, bool(ok)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--sign", nargs="*", default=[])
    a = ap.parse_args()
    sign = {}
    for f in a.sign:
        for g in glob.glob(f):
            sign.update(json.load(open(g)))
    print("# Stage A of PLAN.md: static test of the near-field correction (V1, Donahue & McMichael 2007)\n")
    for f in a.files:
        r = json.load(open(f))
        geo, q = r["geometry"], r["q"]
        print("## Geometry %s, N %d, q %d (%s)\n" % (geo["key"], geo["N"], q, f.split("/")[-1]))
        r0 = r["R0"]
        print("R0: %s. Voxel cells %d, magnetic V1 cells %d, modified %d, new %d; φ_vox %.5f, φ_f %.5f.\n"
              % ("passed" if r0["ok"] else "FAILED", r0["n_voxel"], r0["n_magnetic"], r0["n_modified"], r0["n_new"],
                 r0["phi_vox"], r0["phi_f"]))
        r1 = r["R1"]
        print("R1: %s. Relative deviation from φ(1 − φ)/3: %s.\n" % (
            "passed" if r1["ok"] else "FAILED",
            "; ".join("%s ref %+.1e, voxel %+.1e" % (k, v["ref_rel"], v["vox_rel"]) for k, v in r1["cases"].items())))
        ps = r["p"]
        for name, d in r["fields"].items():
            if name.startswith("uniform"):
                line = "; ".join("P%d %+.1e" % (p, (d["E_P%d" % p] - d["E_ref"]) / d["E_ref"]) for p in ps if "E_P%d" % p in d)
                print("- %s: E_ref %.6f, E_vox %.6f, E_frac %.6f; window error of P (relative): %s" % (
                    name, d["E_ref"], d["E_vox"], d["E_frac"], line))
        print()
        print("| state | E_K | S = E_vox − E_ref | S / |E_K| | E_frac − E_ref | RMS h_vox − h_ref (layer) | RMS h_frac − h_ref |")
        print("|---|---|---|---|---|---|---|")
        states = [n for n in r["fields"] if not n.startswith("uniform")]
        for n in states:
            d, EK = r["fields"][n], r["runs"][n]["E_K"]
            print("| %s | %.3e | %.3e | %.3f | %.3e | %.3e | %.3e |" % (n, EK, d["S"], d["S"] / abs(EK),
                                                                 d["E_frac"] - d["E_ref"], d["rms_layer_vox"],
                                                                 d["rms_layer_frac"]))
        print()
        print("| form | p | per state: e_E / |ΔE|/|E_K| / e_h | criteria in use | R2 in all states | R3 (λ_min / λ_max) |")
        print("|---|---|---|---|---|---|")
        for form in ("P", "L"):
            for p in ps:
                if not all("E_%s%d" % (form, p) in r["fields"][n] for n in states):
                    continue
                cells, use, ok_all = [], set(), True
                for n in states:
                    c, used, ok = crit(r["fields"][n], form, p, r["runs"][n]["E_K"])
                    cells.append("%.3f / %.4f / %.3f" % (c["e_E"], c["dE_EK"], c["e_h"]))
                    use.update(used)
                    ok_all &= ok
                key = "%s_q%d_%s%d" % (geo["key"], q, form, p)
                s = sign.get(key)
                s_txt = ("%.2e (%s)" % (s["ratio"], s["reading"])) if s else "not computed"
                print("| %s | %d | %s | %s | %s | %s |" % (form, p, "; ".join(cells), ", ".join(sorted(use)),
                                                           "yes" if ok_all else "no", s_txt))
        print()


if __name__ == "__main__":
    main()
