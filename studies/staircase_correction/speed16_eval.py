#!/usr/bin/env python
"""
PLAN 16, final speed measurement (d24.sh, d25.sh): the json files of prof_step.py (--trace) in ROOT/speed and
ROOT/wq. For each configuration and code (the start code M_9e0 / B_c63, ABP, AS5): ms per RKF45 stage, ms per sample
interval and tries per sample interval (mean +- SE, n 3), set-up time, peak GPU memory, the halves rule of section 14
(ms per stage of the two halves within 10 %); the gains ABP / AS5 and start code / AS5 with their SE (error
propagation of the two means); whether the 3 runs of a code have bitwise-equal traces (t, tries, <Ms m_x>, accepted
steps); with --w5 the AS5 traces against the AS3 traces of W5 (PLAN 15) and the ABP traces against those of W5.

    python speed16_eval.py ROOT [--w5 W5_FOLDER] [--json OUT]
"""
import argparse
import json
import math
import pathlib

CONF = ("off108", "off156", "q4_108", "q4_156", "q8_108", "q8_156", "q16_108", "q16_156")
OLD = {"off108": "M_9e0", "off156": "M_9e0", "q4_108": "B_c63", "q4_156": "B_c63"}


def mean_se(xs):
    n = len(xs)
    m = sum(xs) / n
    s = math.sqrt(sum((x - m) ** 2 for x in xs) / (n - 1)) if n > 1 else float("nan")
    return m, s / math.sqrt(n), n


def ratio(a, b):
    (ma, sa, _), (mb, sb, _) = a, b
    g = ma / mb
    return g, g * math.sqrt((sa / ma) ** 2 + (sb / mb) ** 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--w5", default=None)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    root = pathlib.Path(a.root)
    files = {f.stem: f for d in ("speed", "wq") for f in (root / d).glob("*.json") if "_r" in f.stem}
    out = {}
    for c in CONF:
        res = {}
        for t in (OLD.get(c), "ABP", "AS5"):
            if t is None:
                continue
            js = [json.loads(files[k].read_text()) for k in sorted(files) if k.startswith("%s_%s_r" % (c, t))]
            js = [j for j in js if "plain" in j]
            if not js:
                continue
            r = {k: mean_se([j["plain"][k] for j in js]) for k in ("ms_per_stage", "ms_per_sample", "tries_per_sample")}
            r["setup_s"] = [round(j["setup_s"], 1) for j in js]
            r["peak_alloc_GiB"] = max(j["peak_alloc_GiB"] for j in js)
            r["halves_max"] = max(j["halves_rel_diff"] for j in js)
            tr = [j["trace"] for j in js]
            r["traces_equal"] = all(x == tr[0] for x in tr)
            if a.w5 and t in ("ABP", "AS5"):
                ref = "AS3" if t == "AS5" else "ABP"
                w = [pathlib.Path(a.w5) / ("%s_%s_r%d.json" % (c, ref, i)) for i in (1, 2, 3)]
                w = [json.loads(p.read_text())["trace"] for p in w if p.exists()]
                r["traces_equal_W5_" + ref] = (all(x == tr[0] for x in w), len(w)) if w else None
            res[t] = r
        if "AS5" in res:
            for t in (OLD.get(c), "ABP"):
                if t in res:
                    res["gain_%s_AS5_stage" % t] = ratio(res[t]["ms_per_stage"], res["AS5"]["ms_per_stage"])
                    res["gain_%s_AS5_sample" % t] = ratio(res[t]["ms_per_sample"], res["AS5"]["ms_per_sample"])
        out[c] = res
        line = "%-8s" % c
        for t in (OLD.get(c), "ABP", "AS5"):
            if t in res:
                m, s, n = res[t]["ms_per_stage"]
                line += " | %s %.3f +- %.3f ms/stage (n %d, tries %.2f, set-up %s s, halves %.1f %%, traces equal %s%s)" % (
                    t, m, s, n, res[t]["tries_per_sample"][0], res[t]["setup_s"], 100 * res[t]["halves_max"],
                    res[t]["traces_equal"], "" if not res[t].get("traces_equal_W5_" + ("AS3" if t == "AS5" else "ABP"))
                    else ", = W5 %s" % ("AS3" if t == "AS5" else "ABP"))
        for t in (OLD.get(c), "ABP"):
            k = "gain_%s_AS5_stage" % t
            if k in res:
                line += " | gain %s/AS5 %.3f +- %.3f per stage, %.3f +- %.3f per sample interval" % (
                    t, *res[k], *res["gain_%s_AS5_sample" % t])
        print(line)
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
