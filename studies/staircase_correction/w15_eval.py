#!/usr/bin/env python
"""
PLAN 15, evaluation of W4 (dynamics against a double-precision reference), W5 (speed) and T5 (configurations), on box
11 (PROTOKOLL 3.44). Reads the json files of prof_step.py (--trace) in /root/d7, /root/d7c, /root/d10 and /root/w5.

W4, for each configuration off108, off156, q4_108, q4_156:
  - stop checks: 126 samples in each run, t equal at each sample in all runs (1e-12 relative);
  - check of Dref: the largest difference of <Ms m_x> of ABP double (off; q4 from 3.41) and of AS3 double against
    Dref at most 1e-3 A/m, and for AS3 also at most max(10 x that of ABP double, 1e-5 A/m);
  - groups ABP (unperturbed run of W5 and the perturbed runs p1, p2 of d10a.sh) and AS3 (W5 r1, p1, p2): fewer than
    3 distinct traces (t, tries, <Ms m_x>, accepted steps) in a group gives "not determined";
  - e: the largest distance of <Ms m_x> to Dref over the 126 sample intervals and the RMS; e of a group: the largest
    over its realisations; pass: e_AS3 at most 2 e_ABP for both;
  - information: e of the old code (D7 runs of 9e0ea7c or c63e862 and their null runs), the accepted steps per
    sample interval of AS3 against ABP, and whether the 3 unperturbed AS3 runs of W5 are bitwise equal.
W5: ms per sample interval and per stage (mean +- SE, n 3), gain ABP / AS3 +- SE, tries, set-up time, peak memory,
the halves rule of section 14 (ms per stage of the two halves within 10 %).
T5: config.json of AS3 against ABP (W5 r1) differs only in out, magnumnp_path and git_commit.

    python w15_eval.py --json /root/w5/w15_eval.json
"""
import argparse
import json
import math
import pathlib

R = pathlib.Path("/root")
CONFS = ["off108", "off156", "q4_108", "q4_156"]
SPEED = CONFS + ["q8_108", "q8_156"]
DREF = {"off108": "d10/off108_M9e0_double", "off156": "d10/off156_M9e0_double",
        "q4_108": "d7c/q4_108_c63_double", "q4_156": "d7c/q4_156_c63_double"}
ABP_DOUBLE = {"off108": "d10/off108_ABP_double", "off156": "d10/off156_ABP_double",
              "q4_108": "d7c/q4_108_ABP_double", "q4_156": "d7c/q4_156_ABP_double"}
OLD = {"off108": ["d7/off108_old_r1", "d7/off108_null_reverse", "d7/off108_null_rotate"],
       "off156": ["d7/off156_old_r1", "d7/off156_null_reverse", "d7/off156_null_rotate"],
       "q4_108": ["d7/q4_108_old_r%d" % r for r in (1, 2, 3)] + ["d7/q4_108_null_reverse", "d7/q4_108_null_rotate"],
       "q4_156": ["d7/q4_156_old_r%d" % r for r in (1, 2, 3)] + ["d7/q4_156_null_reverse", "d7/q4_156_null_rotate"]}


def load(name):
    return json.loads((R / (name + ".json")).read_text())


def trace(name):
    return load(name)["trace"]


def dist(T, D):
    d = [abs(a[2] - b[2]) for a, b in zip(T, D)]
    return max(d), math.sqrt(sum(x * x for x in d) / len(d))


def mean_se(xs):
    n = len(xs)
    m = sum(xs) / n
    s = math.sqrt(sum((x - m) ** 2 for x in xs) / (n - 1)) if n > 1 else float("nan")
    return m, s / math.sqrt(n)


def w4(c, out):
    D = trace(DREF[c])
    runs = {"ABP": ["w5/%s_ABP_r1" % c, "d10/%s_ABP_p1" % c, "d10/%s_ABP_p2" % c],
            "AS3": ["w5/%s_AS3_r1" % c, "w5/%s_AS3_p1" % c, "w5/%s_AS3_p2" % c]}
    T = {k: [trace(f) for f in v] for k, v in runs.items()}
    Tall = [D, trace(ABP_DOUBLE[c]), trace("d10/%s_AS3_double" % c)] + T["ABP"] + T["AS3"]
    stop = []
    for t in Tall:
        if len(t) != 126:
            stop.append("sample count %d" % len(t))
        elif any(abs(a[0] - b[0]) > 1e-12 * abs(b[0]) for a, b in zip(t, D)):
            stop.append("t differs")
    r = {"stop": stop}
    dA = max(abs(a[2] - b[2]) for a, b in zip(trace(ABP_DOUBLE[c]), D))
    dS = max(abs(a[2] - b[2]) for a, b in zip(trace("d10/%s_AS3_double" % c), D))
    r["dref"] = {"ABP_double": dA, "AS3_double": dS, "bound_AS3": max(10 * dA, 1e-5),
                 "pass": dA <= 1e-3 and dS <= 1e-3 and dS <= max(10 * dA, 1e-5)}
    distinct = {}
    for k, ts in T.items():
        keys = [json.dumps(t) for t in ts]
        distinct[k] = len(set(keys))
    r["distinct"] = distinct
    e = {}
    for k, ts in T.items():
        ds = [dist(t, D) for t in ts]
        e[k] = {"largest": max(x[0] for x in ds), "rms": max(x[1] for x in ds),
                "each": [[round(x[0], 5), round(x[1], 5)] for x in ds]}
    r["e"] = e
    old = [trace(f) for f in OLD[c]]
    ds = [dist(t, D) for t in old]
    r["e_old_information"] = {"n": len(old), "distinct": len(set(json.dumps(t) for t in old)),
                              "largest": max(x[0] for x in ds), "rms": max(x[1] for x in ds)}
    sa = [T["ABP"][0][i][3] - (T["ABP"][0][i - 1][3] if i else 0) for i in range(126)]
    ss = [T["AS3"][0][i][3] - (T["AS3"][0][i - 1][3] if i else 0) for i in range(126)]
    r["steps_information"] = {"largest_diff_per_interval": max(abs(a - b) for a, b in zip(sa, ss)),
                              "sum_ABP": sum(sa), "sum_AS3": sum(ss)}
    w5r = [json.dumps(trace("w5/%s_AS3_r%d" % (c, i))) for i in (1, 2, 3)]
    r["AS3_unperturbed_runs_bitwise_equal"] = len(set(w5r)) == 1
    if stop or not r["dref"]["pass"]:
        r["verdict"] = "STOP (check of Dref or stop check)"
    elif min(distinct.values()) < 3:
        r["verdict"] = "not determined (fewer than 3 distinct realisations)"
    else:
        ok = e["AS3"]["largest"] <= 2 * e["ABP"]["largest"] and e["AS3"]["rms"] <= 2 * e["ABP"]["rms"]
        r["verdict"] = "PASS" if ok else "FAIL"
    r["ratio_AS3_ABP"] = {"largest": e["AS3"]["largest"] / e["ABP"]["largest"], "rms": e["AS3"]["rms"] / e["ABP"]["rms"]}
    out["w4"][c] = r
    print("W4 %-7s %-55s Dref: ABP double %.1e, AS3 double %.1e A/m; distinct %s; e largest ABP %.4f AS3 %.4f (x%.2f), "
          "RMS ABP %.4f AS3 %.4f (x%.2f); old %.4f / %.4f (n %d); steps diff %d; AS3 W5 runs bitwise equal %s"
          % (c, r["verdict"], dA, dS, distinct, e["ABP"]["largest"], e["AS3"]["largest"], r["ratio_AS3_ABP"]["largest"],
             e["ABP"]["rms"], e["AS3"]["rms"], r["ratio_AS3_ABP"]["rms"], r["e_old_information"]["largest"],
             r["e_old_information"]["rms"], r["e_old_information"]["n"],
             r["steps_information"]["largest_diff_per_interval"], r["AS3_unperturbed_runs_bitwise_equal"]))


def w5(c, out):
    J = {t: [load("w5/%s_%s_r%d" % (c, t, i)) for i in (1, 2, 3)] for t in ("ABP", "AS3")}
    r = {}
    for t, js in J.items():
        r[t] = {k: mean_se([j["plain"][k] for j in js]) for k in ("ms_per_sample", "ms_per_stage", "tries_per_sample")}
        r[t]["setup_s"] = [round(j["setup_s"], 1) for j in js]
        r[t]["peak_alloc_GiB"] = max(j["peak_alloc_GiB"] for j in js)
        r[t]["halves_max"] = max(j["halves_rel_diff"] for j in js)
    for k in ("ms_per_sample", "ms_per_stage"):
        (a, sa), (b, sb) = r["ABP"][k], r["AS3"][k]
        g = a / b
        r["gain_" + k] = (g, g * math.sqrt((sa / a) ** 2 + (sb / b) ** 2))
    r["halves_rule_ok"] = max(r["ABP"]["halves_max"], r["AS3"]["halves_max"]) <= 0.10
    out["w5"][c] = r
    print("W5 %-7s ABP %8.1f +- %.1f ms per sample interval (%.3f ms/stage, %.2f tries), AS3 %8.1f +- %.1f (%.3f "
          "ms/stage, %.2f tries); gain %.3f +- %.3f per sample interval, %.3f +- %.3f per stage; set-up ABP %s AS3 %s s; "
          "peak %.2f / %.2f GiB; halves max %.1f %%"
          % (c, *r["ABP"]["ms_per_sample"], r["ABP"]["ms_per_stage"][0], r["ABP"]["tries_per_sample"][0],
             *r["AS3"]["ms_per_sample"], r["AS3"]["ms_per_stage"][0], r["AS3"]["tries_per_sample"][0],
             *r["gain_ms_per_sample"], *r["gain_ms_per_stage"], r["ABP"]["setup_s"], r["AS3"]["setup_s"],
             r["ABP"]["peak_alloc_GiB"], r["AS3"]["peak_alloc_GiB"],
             100 * max(r["ABP"]["halves_max"], r["AS3"]["halves_max"])))


def t5(c, out):
    a = json.loads((R / ("w5/run_%s_ABP_r1/config.json" % c)).read_text())
    b = json.loads((R / ("w5/run_%s_AS3_r1/config.json" % c)).read_text())
    diff = sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))
    ok = set(diff) <= {"out", "magnumnp_path", "git_commit"}
    out["t5"][c] = {"diff": diff, "pass": ok}
    print("T5 %-7s %s %s" % (c, "PASS" if ok else "FAIL", diff))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    out = {"w4": {}, "w5": {}, "t5": {}}
    for c in CONFS:
        w4(c, out)
    for c in SPEED:
        w5(c, out)
        t5(c, out)
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
