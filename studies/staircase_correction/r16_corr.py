#!/usr/bin/env python
"""
PLAN 16 C1 (information): how new are the data of the perturbed states? For each quantity, the correlation over the
16 states of the bias b (and of e) of one tree between the original states and the states with m + 1e-7 r
(perturb_ref.py), and the ratio of the RMS of b and of e. A correlation near 1 means that C1 repeats the original
data for this quantity (the rounding pattern comes from the code and the state, not from the last bits of m); near 0
means a new draw.

    python r16_corr.py ORIGINAL.pt PERTURBED.pt
"""
import math
import statistics
import sys

import torch


def corr(x, y):
    mx, my = statistics.mean(x), statistics.mean(y)
    sx = math.sqrt(sum((a - mx) ** 2 for a in x))
    sy = math.sqrt(sum((b - my) ** 2 for b in y))
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy) if sx > 0 and sy > 0 else float("nan")


def main():
    o = torch.load(sys.argv[1], map_location="cpu", weights_only=False)
    p = torch.load(sys.argv[2], map_location="cpu", weights_only=False)
    ks = sorted(o["states"])
    if sorted(p["states"]) != ks:
        raise SystemExit("r16_corr: other states")
    rms = lambda v: math.sqrt(sum(x * x for x in v) / len(v))
    print("%-34s %8s %8s %10s %10s" % ("quantity", "corr b", "corr e", "RMS b p/o", "RMS e p/o"))
    for q in o["states"][ks[0]]["q"]:
        if q not in p["states"][ks[0]]["q"] or "b" not in o["states"][ks[0]]["q"][q]:
            continue
        bo = [o["states"][k]["q"][q]["b"] for k in ks]
        bp = [p["states"][k]["q"][q]["b"] for k in ks]
        eo = [o["states"][k]["q"][q]["e"] for k in ks]
        ep = [p["states"][k]["q"][q]["e"] for k in ks]
        print("%-34s %8.3f %8.3f %10.3f %10.3f" % (q, corr(bo, bp), corr(eo, ep), rms(bp) / rms(bo) if rms(bo) else
                                                  float("nan"), rms(ep) / rms(eo) if rms(eo) else float("nan")))


if __name__ == "__main__":
    main()
