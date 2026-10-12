#!/usr/bin/env python
"""
Test of option O2 before FREQ_GO (PROTOKOLL 7.72): two runs of the same job line, one with --acc_harm 7 (h7) and one
without (h0). Checks the pre-registered criteria 1 … 5 and 7 and prints the numbers for criterion 6.

    python o2_check.py <run h7> <run h0>
"""
import filecmp
import json
import pathlib
import subprocess
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent


def main(h7, h0):
    h7, h0 = pathlib.Path(h7), pathlib.Path(h0)
    ok = {}
    ok["1 DONE in both runs"] = (h7 / "DONE").exists() and (h0 / "DONE").exists()
    s7 = sorted(h7.glob("samples_*.csv"))
    ok["2 samples equal"] = bool(s7) and all(filecmp.cmp(f, h0 / f.name, shallow=False) for f in s7)
    for d in (h7, h0):
        r = subprocess.run([sys.executable, str(HERE / "acc_reduce.py"), str(d)], capture_output=True, text=True)
        print("acc_reduce %s: exit %d\n%s" % (d.name, r.returncode, (r.stdout + r.stderr)[-1500:]))
        ok["7 acc_reduce %s" % d.name] = r.returncode == 0
    a7 = {s["stage"]: s for s in json.loads((h7 / "acc.json").read_text())["stages"]}
    a0 = {s["stage"]: s for s in json.loads((h0 / "acc.json").read_text())["stages"]}
    summ = {s["name"]: s for s in json.loads((h7 / "summary.json").read_text())["stages"]}
    summ0 = {s["name"]: s for s in json.loads((h0 / "summary.json").read_text())["stages"]}
    for st, a in sorted(a7.items()):
        ok["3 p_box/p_dis %s" % st] = a["check_ok"]
        hm = a["harm"]
        cyc = [c for c in summ[st]["cycles"] if "acc" in c]
        all_calm = all(c["acc"]["calm"] for c in cyc)
        print("%s: calm cycles %s of %d accumulated; R_k %s; rest %.4f; Bessel max %.4f, cells over 1.001: %.4f"
              % (st, hm["cycles_kept"], hm["cycles_acc"], " ".join("%.4f" % x for x in hm["R_k"]), hm["R_rest"],
                 hm["bessel_max"], hm["bessel_frac_over_1e-3"]))
        ok["5 Bessel %s" % st] = sum(hm["R_k"]) <= 1.0 + 1e-6 and hm["bessel_frac_over_1e-3"] < 0.01
        if all_calm and st in a0:
            z7 = np.load(next(h7.glob("acc_*_%s.npz" % st))) if any(h7.glob("acc_*_%s.npz" % st)) else None
            z0 = np.load(next(h0.glob("acc_*_%s.npz" % st))) if any(h0.glob("acc_*_%s.npz" % st)) else None
            if z7 is not None and z0 is not None:
                rel = {k: float(np.abs(z7[k] - z0[k]).max() / max(np.abs(z0[k]).max(), 1e-300)) for k in ("p", "m0", "m1")}
                print("   %s all calm: max relative difference h7 against h0: %s; R1 %.6f / %.6f"
                      % (st, rel, a["R1"], a0[st]["R1"]))
                ok["4 equal to h0 %s" % st] = max(rel.values()) < 1e-5 and abs(a["R1"] / a0[st]["R1"] - 1) < 1e-5
            else:
                print("   %s: npz deleted already, compare acc.json only: R1 %.6f / %.6f" % (st, a["R1"], a0[st]["R1"]))
                ok["4 equal to h0 %s" % st] = abs(a["R1"] / a0[st]["R1"] - 1) < 1e-5
        else:
            print("   %s: not all accumulated cycles calm (%s): criterion 4 not applicable"
                  % (st, [c["acc"]["calm"] for c in cyc]))
        w7 = np.mean([c["wall_s"] for c in summ[st]["cycles"]])
        w0 = np.mean([c["wall_s"] for c in summ0[st]["cycles"]])
        print("   6 wall per cycle: %.0f s with harmonics, %.0f s without (× %.2f)" % (w7, w0, w7 / w0))
    print()
    for k, v in ok.items():
        print("%-28s %s" % (k, "PASS" if v else "FAIL"))
    print("ALL PASS" if all(ok.values()) else "NOT ALL PASS")
    return 0 if all(ok.values()) else 1


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:3]))
