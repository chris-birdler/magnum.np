#!/usr/bin/env python
"""
Evaluation of K1 of PLAN 14 (windows over the cycle) from the files of step_dump.py --k1_ref / --k1_arm.

    python k_eval.py --ref ref.pt --pair NAME OLD NEW CTRL [--null NULL] [--align_a4] [--bound 1e-3]
                     [--null_bound 1e-4] [--json out.json]
    python k_eval.py --a4_count ref.pt
    python k_eval.py --same A B            (resync check: tries and m bitwise equal in all windows)

The reference file can also be an arm (OLD): it has the same windows.

Determined try: OLD, NEW and CTRL have the same branch tuple (step_dump.BRANCH) and the same sample interval in all
tries from the checkpoint up to this try. Determined point: the end of interval k0 + p (p in the full-m points of the
file), if all tries up to it are determined and the three arms have the same number of tries up to it.
--align_a4 (OLD without A4, NEW and CTRL with A4): the tries past the end of OLD are left out, and the A4 clamp and
the past-end flag are not compared (the tries are aligned by accepted step).

Criterion 1 at each determined point: rho_all (RMS of |dm| over the cells with |m0| above 0), rho_part (the same over
the particle cells) and dm_max (largest |dm| over the particle cells) of NEW - OLD, each at most --bound times the same
measure of CTRL - OLD. Criterion 2: the RMS over the tries up to the point of (norm_NEW - norm_OLD) / norm_OLD, at most
--bound times the same value of CTRL. Pass: each window has a determined point, and both criteria hold at all
determined points. A window without a determined point is "not determined". The null run (--null) is evaluated in the
same way with NULL in place of NEW and --null_bound (pilot stop rule of PLAN 14).
"""
import argparse
import json
import math
import sys

import torch

F = ("n", "dt", "norm", "dt_opt", "remaining", "accept", "a4_clamp", "past_end", "norm_class", "c_minscale",
     "c_maxscale", "c_one", "c_maxstep", "retry")
I_N, I_NORM, I_PAST = F.index("n"), F.index("norm"), F.index("past_end")
KEY_ALL = [F.index(k) for k in F[5:]]
KEY_A4 = [F.index(k) for k in F[5:] if k not in ("a4_clamp", "past_end")]


def load(p):
    d = torch.load(p, map_location="cpu", weights_only=False)
    if tuple(d["fields"]) != F:
        raise SystemExit("k_eval: %s has other try fields %s" % (p, d["fields"]))
    return d


def rows_of(d, k0, drop_past=False):
    r = d["win"][k0].get("rows", [])
    return [x for x in r if not (drop_past and x[I_PAST])]


def measures(dm, sel_all, sel_part):
    a = dm.double().norm(dim=-1)
    return {"rho_all": float(a[sel_all].pow(2).mean().sqrt()), "rho_part": float(a[sel_part].pow(2).mean().sqrt()),
            "dm_max": float(a[sel_part].max())}


def norm_rms(rows_a, rows_b):
    v = [(b[I_NORM] - a[I_NORM]) / a[I_NORM] for a, b in zip(rows_a, rows_b) if a[I_NORM] > 0]
    return math.sqrt(sum(x * x for x in v) / len(v)) if v else 0.0


def window(old, new, ctrl, k0, key, align_a4, sel_all, sel_part, bound):
    ro, rn, rc = rows_of(old, k0, align_a4), rows_of(new, k0), rows_of(ctrl, k0)
    n_det = 0
    for a, b, c in zip(ro, rn, rc):
        if not (a[I_N] == b[I_N] == c[I_N] and [a[i] for i in key] == [b[i] for i in key] == [c[i] for i in key]):
            break
        n_det += 1
    res = {"tries": [len(ro), len(rn), len(rc)], "determined_tries": n_det, "points": {}}
    for p in sorted(old["win"][k0]["m"]):
        n_end = k0 + p
        cnt = [sum(1 for x in r if x[I_N] <= n_end) for r in (ro, rn, rc)]
        if not (cnt[0] == cnt[1] == cnt[2] and cnt[0] <= n_det and cnt[0] > 0):
            res["points"][p] = {"determined": False, "tries": cnt}
            continue
        mo = old["win"][k0]["m"][p]
        pair = measures(new["win"][k0]["m"][p] - mo, sel_all, sel_part)
        ctl = measures(ctrl["win"][k0]["m"][p] - mo, sel_all, sel_part)
        pair["norm_rel_rms"] = norm_rms(ro[:cnt[0]], rn[:cnt[0]])
        ctl["norm_rel_rms"] = norm_rms(ro[:cnt[0]], rc[:cnt[0]])
        ratio = {k: (pair[k] / ctl[k] if ctl[k] > 0 else (0.0 if pair[k] == 0 else math.inf)) for k in pair}
        res["points"][p] = {"determined": True, "tries": cnt[0], "pair": pair, "control": ctl, "ratio": ratio,
                            "max_ratio": max(ratio.values()), "pass": all(v <= bound for v in ratio.values())}
    det = [v for v in res["points"].values() if v["determined"]]
    res["status"] = "not determined" if not det else ("pass" if all(v["pass"] for v in det) else "fail")
    res["max_ratio"] = max((v["max_ratio"] for v in det), default=None)
    return res


def evaluate(name, ref, old, new, ctrl, align_a4, bound):
    m0 = ref["ck"][min(ref["ck"])]["m"]
    sel_all = m0.double().norm(dim=-1) > 0
    sel_part = ref["particle"] & sel_all if ref.get("particle") is not None else sel_all
    key = KEY_A4 if align_a4 else KEY_ALL
    wins = {k0: window(old, new, ctrl, k0, key, align_a4, sel_all, sel_part, bound) for k0 in sorted(old["win"])}
    st = [w["status"] for w in wins.values()]
    status = "fail" if "fail" in st else ("not determined" if "not determined" in st else "pass")
    mr = [w["max_ratio"] for w in wins.values() if w["max_ratio"] is not None]
    return {"name": name, "bound": bound, "status": status, "windows": wins,
            "not_determined": st.count("not determined"), "max_ratio": max(mr) if mr else None,
            "cells": int(sel_all.sum()), "particle_cells": int(sel_part.sum())}


def show(r):
    print("%-14s %s (bound %.0e): largest ratio %s; windows not determined %d of %d"
          % (r["name"], r["status"].upper(), r["bound"], "%.2e" % r["max_ratio"] if r["max_ratio"] is not None
             else "-", r["not_determined"], len(r["windows"])))
    for k0, w in r["windows"].items():
        pts = " ".join("p%d %s" % (p, ("%.1e" % v["max_ratio"]) if v["determined"] else "nd")
                       for p, v in w["points"].items())
        print("    k0 %3d: %-14s tries %s, determined %d; %s" % (k0, w["status"], w["tries"], w["determined_tries"], pts))


def a4_count(p):
    d = load(p)
    rows = d["rows"]
    past = [x for x in rows if x[I_PAST]]
    case = [x for x in past if x[F.index("retry")] == 1]
    print("A4 case in %s: %d tries past the end, of them %d with dt_opt below the rest of the interval (%d tries in all)"
          % (p, len(past), len(case), len(rows)))
    return {"tries": len(rows), "past_end": len(past), "a4_case": len(case),
            "a4_case_rows": [list(x) for x in case]}


def same(pa, pb):
    a, b = load(pa), load(pb)
    ok = sorted(a["win"]) == sorted(b["win"])
    for k0 in sorted(a["win"]):
        wa, wb = a["win"][k0], b["win"].get(k0, {})
        eq_rows = wa.get("rows") == wb.get("rows")
        eq_m = all(torch.equal(wa["m"][p], wb["m"][p]) for p in wa["m"]) and sorted(wa["m"]) == sorted(wb["m"])
        eq_t = wa["t_end"] == wb["t_end"] and wa["dt_end"] == wb["dt_end"]
        eq_p = (wa["probe_m"] is None and wb["probe_m"] is None) or torch.equal(wa["probe_m"], wb["probe_m"])
        print("k0 %3d: tries equal %s (%d), m equal %s, t and dt equal %s, probe m equal %s"
              % (k0, eq_rows, len(wa.get("rows", [])), eq_m, eq_t, eq_p))
        ok &= eq_rows and eq_m and eq_t and eq_p
    print("RESULT windows bitwise equal: %s" % ok)
    return ok


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref")
    ap.add_argument("--pair", nargs=4, metavar=("NAME", "OLD", "NEW", "CTRL"))
    ap.add_argument("--null")
    ap.add_argument("--align_a4", action="store_true")
    ap.add_argument("--bound", type=float, default=1e-3)
    ap.add_argument("--null_bound", type=float, default=1e-4)
    ap.add_argument("--json")
    ap.add_argument("--a4_count")
    ap.add_argument("--same", nargs=2)
    a = ap.parse_args(argv)
    if a.same:
        sys.exit(0 if same(*a.same) else 1)
    out = {}
    if a.a4_count:
        out["a4"] = a4_count(a.a4_count)
    if a.pair:
        ref = load(a.ref)
        name, po, pn, pc = a.pair
        old, ctrl = load(po), load(pc)
        if a.null:
            out["null"] = evaluate(name + " null", ref, old, load(a.null), ctrl, False, a.null_bound)
            show(out["null"])
        if pn != "-":
            out["pair"] = evaluate(name, ref, old, load(pn), ctrl, a.align_a4, a.bound)
            show(out["pair"])
    if a.json:
        with open(a.json, "w") as f:
            json.dump(out, f, indent=1, default=str)


if __name__ == "__main__":
    main(sys.argv[1:])
