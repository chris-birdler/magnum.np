#!/usr/bin/env python
"""
Summary of dissmap.json (dissmap_batch.py) over the baseline B300_s921 ... s928 (PROTOKOLL 7.23).

    python eval_dissmap.py > results/dissmap.md
"""
import json
import math
import pathlib

import numpy as np

RUNS = pathlib.Path(__file__).resolve().parent / "runs"
AMPS = (9, 20, 35, 50, 70, 100, 150)
KEYS = [("R", "R = p_fund / p_dis"), ("shell_outer", "share in r/R_p > 0.8 (volume 0.49)"),
        ("shell_inner", "share in r/R_p < 0.5 (volume 0.13)"), ("texture_top10", "share in the 10 % cells of largest texture"),
        ("top10", "share in the 10 % cells of largest p_fund")]


def main():
    D = {}
    for s in range(921, 929):
        for st in json.loads((RUNS / ("B300_s%d" % s) / "dissmap.json").read_text())["stages"]:
            D.setdefault(min(AMPS, key=lambda x: abs(x - st["B_peak_mT"])), []).append(st)
    print("| B_peak | n | " + " | ".join(t for _, t in KEYS) + " |")
    print("|---|---|" + "---|" * len(KEYS))
    for mT in AMPS:
        v = D[mT]
        cells = []
        for k, _ in KEYS:
            x = np.array([r[k] for r in v])
            cells.append("%.2f ± %.2f" % (x.mean(), x.std(ddof=1) / math.sqrt(len(x))))
        print("| %d mT | %d | %s |" % (mT, len(v), " | ".join(cells)))
    print("\np_fund: Gilbert dissipation of the local rotation of m at the drive frequency (from the 4 snapshots of "
          "the last cycle); p_dis: measured dissipation of the same cycle. Shares refer to p_fund. Mean ± SE over the "
          "8 seeds. Base point, α 0.1, d/l_ex 300, dx 3.")


if __name__ == "__main__":
    main()
