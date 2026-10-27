#!/usr/bin/env python
"""
Speed result of D7 (PLAN 14b; the speed measurement of PLAN 14) from the json files of prof_step.py --window 0 in one
folder (d7.sh): <conf>_old_r<k>.json, <conf>_new_r<k>.json (k 1 ... 3), kb_<tree>_<N>.json (kernel builds). A run that
failed the halves rule and ran again with 100 intervals of warm-up is <conf>_<arm>_r<k>_rep.json; for the speed it
replaces the run <k>.

For each configuration and arm (n runs): ms per sample interval, ms per RKF45 stage, tries per sample interval and
set-up time, as mean ± standard error (SE = SD / sqrt(n)); the peak GPU memory (largest of the runs); the time of one
drive cycle (ms per sample interval × --samples). For each configuration with an old arm: the gain factor old / new of
the ms per sample interval and of the ms per stage, ± SE (propagated from the two SEs). Halves rule: the ms per stage of
the two halves of a run differ by at most 10 %; a run above it is marked. An arm without a json file is reported with
the last line of its log (for example out of GPU memory).

    python d7_eval.py /root/d7 [--json out.json]
"""
import argparse
import json
import math
import pathlib
import statistics

CONFS = ("off108", "off156", "q4_108", "q4_156", "q8_108", "q8_156")
KEYS = ("ms_per_sample", "ms_per_stage", "tries_per_sample", "setup_s")


def mean_se(v):
    m = statistics.fmean(v)
    return m, (statistics.stdev(v) / math.sqrt(len(v)) if len(v) > 1 else float("nan"))


def run_values(r):
    p = r["plain"]
    return {"ms_per_sample": p["ms_per_sample"], "ms_per_stage": p["ms_per_stage"],
            "tries_per_sample": p["tries_per_sample"], "setup_s": r["setup_s"], "halves": r["halves_rel_diff"],
            "peak_alloc_GiB": r["peak_alloc_GiB"], "peak_reserved_GiB": r["peak_reserved_GiB"], "skip": r["skip"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--runs", type=int, default=3)
    ap.add_argument("--halves", type=float, default=0.10)
    ap.add_argument("--json")
    a = ap.parse_args()
    d = pathlib.Path(a.folder)
    out = {"kernel_builds": {}, "conf": {}}
    for f in sorted(d.glob("kb_*.json")):
        r = json.loads(f.read_text())
        out["kernel_builds"][f.stem] = {"setup_s": r["setup_s"]}
        print("kernel build %-14s set-up %.1f s" % (f.stem[3:], r["setup_s"]))
    for c in CONFS:
        res = out["conf"][c] = {}
        for arm in ("old", "new"):
            runs, marks, missing = [], [], []
            for k in range(1, a.runs + 1):
                f = d / ("%s_%s_r%d.json" % (c, arm, k))
                rep = d / ("%s_%s_r%d_rep.json" % (c, arm, k))
                if rep.exists():
                    f = rep
                if not f.exists():
                    log = d / ("%s_%s_r%d.log" % (c, arm, k))
                    if log.exists():
                        lines = log.read_text().strip().splitlines()
                        missing.append("r%d: %s" % (k, lines[-1][:160] if lines else "empty log"))
                    continue
                r = json.loads(f.read_text())
                v = run_values(r)
                runs.append(v)
                if v["halves"] > a.halves:
                    marks.append("%s: halves differ by %.1f %%" % (f.stem, 100 * v["halves"]))
            if not runs:
                if missing:
                    res[arm] = {"n": 0, "missing": missing}
                    print("%-7s %s: no data; %s" % (c, arm, "; ".join(missing)))
                continue
            samples = int(r["argv"][r["argv"].index("--samples") + 1])
            s = {"n": len(runs), "halves_marks": marks, "missing": missing,
                 "peak_alloc_GiB": max(v["peak_alloc_GiB"] for v in runs),
                 "peak_reserved_GiB": max(v["peak_reserved_GiB"] for v in runs), "samples_per_cycle": samples}
            for key in KEYS:
                s[key] = mean_se([v[key] for v in runs])
            s["s_per_cycle"] = (s["ms_per_sample"][0] * samples / 1e3, s["ms_per_sample"][1] * samples / 1e3)
            res[arm] = s
            print("%-7s %s (n %d): %8.1f ± %5.1f ms per sample interval, %.3f ± %.3f ms per stage, %.2f ± %.2f tries "
                  "per sample interval, set-up %.1f ± %.1f s, peak %.2f GiB allocated (%.2f reserved), %.0f s per "
                  "cycle%s%s" % (c, arm, s["n"], *s["ms_per_sample"], *s["ms_per_stage"], *s["tries_per_sample"],
                                 *s["setup_s"], s["peak_alloc_GiB"], s["peak_reserved_GiB"], s["s_per_cycle"][0],
                                 ("; " + "; ".join(marks)) if marks else "",
                                 ("; missing " + "; ".join(missing)) if missing else ""))
        if res.get("old", {}).get("n") and res.get("new", {}).get("n"):
            for key in ("ms_per_sample", "ms_per_stage"):
                (mo, so), (mn, sn) = res["old"][key], res["new"][key]
                g = mo / mn
                res["gain_" + key] = (g, g * math.sqrt((so / mo) ** 2 + (sn / mn) ** 2))
            print("%-7s gain old / new: %.2f ± %.2f per sample interval, %.2f ± %.2f per stage (n %d old, %d new; SE)"
                  % (c, *res["gain_ms_per_sample"], *res["gain_ms_per_stage"], res["old"]["n"], res["new"]["n"]))
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
