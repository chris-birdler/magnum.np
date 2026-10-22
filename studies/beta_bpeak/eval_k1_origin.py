#!/usr/bin/env python
"""
K1 investigation (PROTOKOLL 7.136, Chris 2026-10-22): where and in which form the loss of the runs K<c>_<k>_sat
occurs, with the methods of the main study. Local, reads summary.json, acc.json and domains.json of the run folders.
Per run and stage (9, 50, 150 mT, the stages with snapshots):
  ln P      P = f * W_loop (mean of the last 2 cycles), W/cm³ of the simulation box
  channels  (7.70, 7.72; acc.json harm): R1 = share at f, R_harm = share at 2f … Kf, R_rest = the remainder
  location  (7.69; acc.json where): concentration c = loss share / volume share for the core (r < 0.5), the rim
            (r >= 0.8), the 10 % of the cells with the largest texture and the 10 % with the largest loss
  walls     (7.125; domains.json half_cycle): A_w (mean wall area), the wall path s = V_sw / A_w of the two half
            cycles (H max <-> H min, H = 0 falling <-> rising), X = A_w * s_mean² (viscous wall motion: loss ~ X)
Per corner: Pearson r of ln P with ln X over the runs; for K1 also without K1_4 and K1_8.

    python eval_k1_origin.py [--runs DIR] [--json FILE] > report.md
"""
import argparse
import json
import pathlib

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
STAGES = {0: "9 mT", 3: "50 mT", 6: "150 mT"}
REGIONS = [("r<0.5", "core"), ("r>=0.8", "rim"), ("texture top 10%", "texture top 10 %"), ("p top 10%", "loss top 10 %")]
SPECIAL = ("K1_4_sat", "K1_8_sat")


def run_record(run, i):
    cfg = json.loads((run / "config.json").read_text())
    summ = json.loads((run / "summary.json").read_text())
    acc = json.loads((run / "acc.json").read_text())["stages"][i]
    dom = [s for s in json.loads((run / "domains.json").read_text())["stages"] if s["stage"] == "amp%02d" % i][0]
    W = np.mean([c["W_loop"] for c in summ["stages"][i]["cycles"][-2:]])
    rec = {"lnP": float(np.log(W * cfg["freq"] / 1e6)), "R1": acc["harm"]["R_k"][0],
           "R_harm": float(sum(acc["harm"]["R_k"][1:])), "R_rest": acc["harm"]["R_rest"]}
    for key, name in REGIONS:
        w = acc["where"][key]
        rec["c " + name] = w["p_share"] / w["vol_share"]
    h = dom["half_cycle"]
    s = float(np.mean(h["wall_path"]))
    rec.update({"A_w": h["A_w_mean"], "s0": h["wall_path"][0], "s1": h["wall_path"][1],
                "lnX": float(np.log(h["A_w_mean"] * s * s))})
    return rec


def pearson(x, y):
    return float(np.corrcoef(x, y)[0, 1]) if len(x) > 2 else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=str(HERE / "runs"))
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    root = pathlib.Path(a.runs)
    names = sorted((d.name for d in root.iterdir() if d.is_dir() and d.name.startswith("K") and
                    d.name.endswith("_sat") and (d / "DONE").exists() and (d / "domains.json").exists()),
                   key=lambda n: (n[:2], int(n.split("_")[1])))
    out = {n: {STAGES[i]: run_record(root / n, i) for i in STAGES} for n in names}
    keys = list(out[names[0]]["9 mT"].keys())
    print("# Loss origin of the runs K<c>_<k>_sat (methods of 7.69, 7.70, 7.125)\n")
    print("c: loss share / volume share of the region (1: the same density as the mean). A_w in l_ex², s in l_ex. "
          "X = A_w s_mean², s_mean = mean of the 2 half cycles.\n")
    for st in STAGES.values():
        print("## %s\n" % st)
        print("| run | " + " | ".join(keys) + " |")
        print("|---|" + "---|" * len(keys))
        for n in names:
            r = out[n][st]
            print("| %s | %s |" % (n, " | ".join(("%.0f" if k == "A_w" else "%.3f") % r[k] for k in keys)))
        print("\nK1: the other runs (range) against K1_4 and K1_8\n")
        oth = [n for n in names if n.startswith("K1_") and n not in SPECIAL]
        print("| quantity | other K1 runs (n %d) | K1_4 | K1_8 |" % len(oth))
        print("|---|---|---|---|")
        for k in keys:
            v = [out[n][st][k] for n in oth]
            sp = ["%.3f" % out[n][st][k] if n in out else "-" for n in SPECIAL]
            print("| %s | %.3f … %.3f | %s |" % (k, min(v), max(v), " | ".join(sp)))
        print("\nPearson r of ln P with ln X over the runs of a corner\n")
        for c in ("K1", "K2", "K3", "K4"):
            ns = [n for n in names if n.startswith(c + "_")]
            r = pearson([out[n][st]["lnX"] for n in ns], [out[n][st]["lnP"] for n in ns])
            line = "- %s: r %.2f (n %d)" % (c, r, len(ns))
            if c == "K1":
                no = [n for n in ns if n not in SPECIAL]
                line += "; without K1_4 and K1_8: r %.2f (n %d)" % (
                    pearson([out[n][st]["lnX"] for n in no], [out[n][st]["lnP"] for n in no]), len(no))
            print(line)
        print()
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
