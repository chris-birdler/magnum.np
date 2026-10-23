#!/usr/bin/env python
"""
Count of single outlier tries in K1 files (finding of the K1 pilot, PLAN 14a): in each window, go try by try through
the reference and all arms while all files have the same sample interval n and dt (relative 1e-6). At each such try,
take the median of the error norms over the files. An outlier is a try whose norm differs from the median by more
than the factor --factor. The count stops in a window at the first try where the files do not agree in n or dt.
With --pair: each arm against the reference alone (the count of an arm stops at its own first difference of n or dt).
With --reps: the repetitions of the windows in each arm file (step_dump --k1_reps) take the place of the files. Rows
with a second try (step_dump --k1_redo) are listed when the two norms differ by more than --redo_rel.

    python glitch_count.py REF.pt ARM1.pt [ARM2.pt ...] [--pair] [--factor 2] [--json out.json]
    python glitch_count.py --reps ARM.pt [ARM2.pt ...] [--json out.json]
"""
import argparse
import json
import pathlib
import statistics

import torch


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--factor", type=float, default=2.0)
    ap.add_argument("--dt_rel", type=float, default=1e-6)
    ap.add_argument("--pair", action="store_true")
    ap.add_argument("--reps", action="store_true")
    ap.add_argument("--redo_rel", type=float, default=1e-6)
    ap.add_argument("--json")
    a = ap.parse_args()
    if a.reps:
        for f in a.files:
            d = torch.load(f, map_location="cpu", weights_only=False)
            R = d["reps"]
            wins = [{k0: w["rows"] for k0, w in d["win"].items()}]
            wins += [{k0: w["reps"][r - 1]["rows"] for k0, w in d["win"].items()} for r in range(1, R)]
            stem = pathlib.Path(f).stem
            group(a, ["%s_r%d" % (stem, r) for r in range(R)], wins, d["fields"], f)
            dm = [x["dm_end"] for w in d["win"].values() for x in w.get("reps", [])]
            print("%s: %d windows x %d integrations; largest |m - m| at the window end against the first: %.3e"
                  % (stem, len(d["win"]), R, max(dm, default=0.0)))
            if "redo_norm" in d["fields"]:
                redo(a, stem, wins, d["fields"])
    else:
        pairs = [[a.files[0], f] for f in a.files[1:]] if a.pair else [a.files]
        for files in pairs:
            D = [torch.load(f, map_location="cpu", weights_only=False) for f in files]
            group(a, [pathlib.Path(f).stem for f in files],
                  [{k0: w.get("rows", []) for k0, w in d["win"].items()} for d in D], D[1]["fields"], None)


def redo(a, stem, wins, fld):
    """Rows with a second try: the two norms and the place of the largest scaled error."""
    i_no, i_rn, i_im, i_ri = fld.index("norm"), fld.index("redo_norm"), fld.index("imax"), fld.index("redo_imax")
    n_redo, diff = 0, []
    for r, w in enumerate(wins):
        for k0, rows in w.items():
            for i, row in enumerate(rows):
                if row[i_rn] == row[i_rn]:
                    n_redo += 1
                    q = row[i_no] / row[i_rn] if row[i_rn] > 0 else float("inf")
                    if abs(q - 1) > a.redo_rel:
                        diff.append((r, k0, i, row[i_no], row[i_rn], row[i_im], row[i_ri]))
    print("%s: %d tries with a second try (norm above 1.1); %d with a different norm (relative %.0e)"
          % (stem, n_redo, len(diff), a.redo_rel))
    for r, k0, i, n1, n2, i1, i2 in diff:
        print("   rep %d k0 %3d try %4d: norm %.6g, second try %.6g; largest entry at %d, second try at %d"
              % (r, k0, i, n1, n2, i1, i2))


def group(a, names, wins, fld, src):
    i_n, i_dt, i_norm, i_acc = fld.index("n"), fld.index("dt"), fld.index("norm"), fld.index("accept")
    ks = sorted(set.intersection(*[set(w) for w in wins]))
    compared = {nm: 0 for nm in names}
    out = {"files": src or names, "factor": a.factor, "outliers": [], "windows": {}}
    worst = {nm: 1.0 for nm in names}
    for k0 in ks:
        rows = [w[k0] for w in wins]
        i = 0
        while all(i < len(r) for r in rows):
            ns = [r[i][i_n] for r in rows]
            dts = [r[i][i_dt] for r in rows]
            if len(set(ns)) > 1 or max(dts) > min(dts) * (1 + a.dt_rel):
                break
            norms = [r[i][i_norm] for r in rows]
            # two files: the reference is the measure, only the arm is judged
            med = statistics.median(norms) if len(rows) > 2 else norms[0]
            for j, nm in enumerate(names):
                if len(rows) == 2 and j == 0:
                    continue
                compared[nm] += 1
                if med > 0 and norms[j] > 0:
                    q = max(norms[j] / med, med / norms[j])
                    worst[nm] = max(worst[nm], q) if q <= a.factor else worst[nm]
                    if q > a.factor:
                        out["outliers"].append({"file": nm, "k0": k0, "try": i, "n": ns[0], "dt": dts[0],
                                                "norm": norms[j], "median": med,
                                                "accept": bool(rows[j][i][i_acc]),
                                                "accept_others": [bool(r[i][i_acc]) for r in rows]})
            i += 1
        out["windows"][k0] = {"tries_compared": i, "tries": [len(r) for r in rows]}
    out["compared"] = compared
    out["largest_ratio_without_outliers"] = worst
    for o in out["outliers"]:
        print("outlier %-10s k0 %3d try %4d n %3d dt %.3e norm %.4g median %.4g accepted %s" % (
            o["file"], o["k0"], o["try"], o["n"], o["dt"], o["norm"], o["median"], o["accept"]))
    for nm in names:
        print("%-10s tries compared %6d, outliers %d, largest ratio of the other tries %.6f" % (
            nm, compared[nm], sum(1 for o in out["outliers"] if o["file"] == nm), worst[nm]))
    if a.json:
        p = pathlib.Path(a.json)
        p = p.with_name(p.stem + "_" + names[-1] + p.suffix) if a.pair else p
        p.write_text(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
