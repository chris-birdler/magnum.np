#!/usr/bin/env python
"""
Evaluation of K1 of PLAN 14 / 14a (windows over the cycle) from the files of step_dump.py --k1_ref / --k1_arm.

    python k_eval.py --ref ref.pt --pair NAME OLD NEW --base BASE --ctrl CTRL --null N1 N2 [--align_a4]
                     [--json out.json]
    python k_eval.py --ref ref.pt --pilot NAME --base BASE --ctrl CTRL --null N1 N2 [--json out.json]
    python k_eval.py --a4_count ref.pt
    python k_eval.py --same A B            (resync check: tries and m bitwise equal in all windows)

A reference file can also serve as an arm (it has the same windows). BASE is the old code of the control CTRL and of
the null runs N1 (reverse order of the field sum) and N2 (rotated order); each of their measures is the difference to
BASE. The measures of the pair are NEW - OLD.

Determined point: the end of interval k0 + p (p in the full-m points of the files), if all arms (OLD, NEW, BASE,
CTRL and the null runs) have the same sample interval and branch tuple (step_dump.BRANCH, with the decision) in all
tries up to it and the same number of tries up to it. --align_a4 (one side without A4): the tries past the end are
left out of all arms, the A4 clamp and the past-end flag are not compared, and in a window with an A4 case (a try past
the end with dt_opt below the rest of the interval, retry branch 1) the points after the case are not judged.

Criteria (PLAN 14a, T4):
  D:  each window has at least one determined point (--align_a4: at least --a4_min_windows windows).
  C1: at each determined point of the windows k0 >= --c1_from, rho_all, rho_part and dm_max of the pair are each at
      most --c1 times the same measure of CTRL; the same for the sample rows of the runs (M_par, Mx, My, Mz, p_dis;
      the diagnostic path), with a tolerance of 2 units of the 10th digit of the CSV.
  C2: for rho_all, rho_part, dm_max and norm_rms (RMS over the tries up to the point of the relative difference of the
      error norms), N = RMS of the null runs N1 and N2, or the kernel null run N3 if it is larger, at least 2.2e-16:
      the geometric mean of pair / N over the values of each window is at most --gm_window, over all windows at most
      --gm. A pair value of 0 counts as a pass and is not in the means.
Pilot (before the patch arms): the null runs are determined in each window; null / control in C1 at most --pilot_c1 at
the points with k0 >= --c1_from; the geometric mean of N2 / N1 in --pilot_gm over all windows and in
--pilot_gm_window in each window.
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
M_KEYS = ("rho_all", "rho_part", "dm_max")
C2_KEYS = M_KEYS + ("norm_rms",)
N_FLOOR = 2.2e-16


def load(p):
    d = torch.load(p, map_location="cpu", weights_only=False)
    if tuple(d["fields"]) != F:
        raise SystemExit("k_eval: %s has other try fields %s" % (p, d["fields"]))
    return d


class Geo:
    """Cell sets of the reference: all cells with |m0| above 0, and the particle cells."""

    def __init__(self, ref):
        m0 = ref["ck"][min(ref["ck"])]["m"]
        self.all = m0.double().norm(dim=-1) > 0
        self.part = ref["particle"] & self.all if ref.get("particle") is not None else self.all


def measures(dm, geo):
    a = dm.double().norm(dim=-1)
    return {"rho_all": float(a[geo.all].pow(2).mean().sqrt()), "rho_part": float(a[geo.part].pow(2).mean().sqrt()),
            "dm_max": float(a[geo.part].max())}


def norm_rms(rows_a, rows_b):
    v = [(b[I_NORM] - a[I_NORM]) / a[I_NORM] for a, b in zip(rows_a, rows_b) if a[I_NORM] > 0]
    return math.sqrt(sum(x * x for x in v) / len(v)) if v else 0.0


def det_tries(rows, key):
    """Number of tries from the checkpoint on with the same interval and branch in all arms."""
    n = 0
    for xs in zip(*rows):
        if len({x[I_N] for x in xs}) != 1 or len({tuple(x[i] for i in key) for x in xs}) != 1:
            break
        n += 1
    return n


def samples(d):
    """Rows of the sample CSV of the run of an arm (row n: after sample interval n), or None."""
    import glob
    import numpy as np
    argv = d.get("argv", [])
    if "--out" not in argv:
        return None
    f = sorted(glob.glob(argv[argv.index("--out") + 1] + "/samples_*.csv"))
    if not f:
        return None
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        v = np.loadtxt(f[0], delimiter=",", comments="#", ndmin=2)
    return v if v.shape[0] > 1 else None


S_COLS = {"M_par": 3, "Mx": 4, "My": 5, "Mz": 6, "p_dis": 8}


def window(arms, k0, key, align, geo, n_case):
    """arms: dict name -> file (old, new, base, ctrl, n1, n2, n3; some may be missing). All measures at the points.
    n_case: first sample interval with an A4 case in this window (points after it are not judged), or None."""
    rows = {a: [x for x in d["win"][k0].get("rows", []) if not (align and x[I_PAST])] for a, d in arms.items()}
    nd = det_tries([rows[a] for a in arms], key)
    res = {"tries": {a: len(r) for a, r in rows.items()}, "det_tries": nd, "a4_case_interval": n_case, "points": {}}
    for p in sorted(arms["base"]["win"][k0]["m"]):
        n_end = k0 + p
        cnt = {a: sum(1 for x in r if x[I_N] <= n_end) for a, r in rows.items()}
        c0 = cnt["base"]
        det = c0 > 0 and all(v == c0 for v in cnt.values()) and c0 <= nd and (n_case is None or n_end < n_case)
        pt = {"det": det, "tries": c0}
        if det:
            mb = arms["base"]["win"][k0]["m"][p]
            for a in [x for x in ("ctrl", "n1", "n2", "n3") if x in arms]:
                pt[a] = measures(arms[a]["win"][k0]["m"][p] - mb, geo)
                pt[a]["norm_rms"] = norm_rms(rows["base"][:c0], rows[a][:c0])
            if "new" in arms:
                pt["pair"] = measures(arms["new"]["win"][k0]["m"][p] - arms["old"]["win"][k0]["m"][p], geo)
                pt["pair"]["norm_rms"] = norm_rms(rows["old"][:c0], rows["new"][:c0])
            nn = {k: math.sqrt(0.5 * (pt["n1"][k] ** 2 + pt["n2"][k] ** 2)) for k in C2_KEYS}
            if "n3" in pt:
                nn = {k: max(nn[k], pt["n3"][k]) for k in C2_KEYS}
            pt["N"] = {k: max(v, N_FLOOR) for k, v in nn.items()}
        res["points"][p] = pt
    return res


def ratio(a, b):
    return a / b if b > 0 else (0.0 if a == 0 else math.inf)


def a4_cases(d, k0):
    """First sample interval of window k0 with an A4 case (a try past the end with dt_opt below the rest)."""
    c = [x[I_N] for x in d["win"][k0].get("rows", []) if x[I_PAST] and x[F.index("retry")] == 1]
    return min(c) if c else None


def gm(v):
    return 10 ** (sum(math.log10(x) for x in v) / len(v)) if v else None


def evaluate(name, ref, arms, align, a):
    geo = Geo(ref)
    key = KEY_A4 if align else KEY_ALL
    wins = {}
    for k0 in sorted(arms["base"]["win"]):
        n_case = None
        if align:
            cs = [c for c in (a4_cases(arms[x], k0) for x in ("old", "new") if x in arms) if c is not None]
            n_case = min(cs) if cs else None
        wins[k0] = window(arms, k0, key, align, geo, n_case)
    out = {"name": name, "windows": wins, "cells": int(geo.all.sum()), "particle_cells": int(geo.part.sum())}
    out["windows_without_point"] = [k0 for k0, w in wins.items() if not any(p["det"] for p in w["points"].values())]
    out["a4_case_windows"] = [k0 for k0, w in wins.items() if w["a4_case_interval"] is not None]
    pts = [(k0, p, v) for k0, w in wins.items() for p, v in w["points"].items() if v["det"]]
    # pilot quantities: the null runs against the control and against each other
    q = {"null_c1": [], "n2_n1_window": {}, "n2_n1_all": []}
    for k0, p, v in pts:
        if k0 >= a.c1_from:
            q["null_c1"].append(max(ratio(v[n][k], v["ctrl"][k]) for n in ("n1", "n2", "n3") if n in v for k in M_KEYS))
        lr = [v["n2"][k] / v["n1"][k] for k in C2_KEYS if v["n1"][k] > 0 and v["n2"][k] > 0]
        q["n2_n1_window"].setdefault(k0, []).extend(lr)
        q["n2_n1_all"].extend(lr)
    out["pilot"] = {"windows_without_point": out["windows_without_point"],
                    "null_c1_max": max(q["null_c1"]) if q["null_c1"] else None,
                    "n2_n1_gm": gm(q["n2_n1_all"]),
                    "n2_n1_gm_window": {k0: gm(v) for k0, v in q["n2_n1_window"].items()}}
    pq = out["pilot"]
    wg = [g for g in pq["n2_n1_gm_window"].values() if g is not None]
    pq["pass"] = (not pq["windows_without_point"] and pq["null_c1_max"] is not None and pq["null_c1_max"] <= a.pilot_c1
                  and pq["n2_n1_gm"] is not None and a.pilot_gm[0] <= pq["n2_n1_gm"] <= a.pilot_gm[1]
                  and all(a.pilot_gm_window[0] <= g <= a.pilot_gm_window[1] for g in wg))
    if "new" not in arms:
        return out
    # D: each window has a determined point (A4: at least --a4_min_windows windows with judged points)
    judged = [k0 for k0, w in wins.items() if any(p["det"] for p in w["points"].values())]
    out["D"] = (len(judged) >= a.a4_min_windows) if align else not out["windows_without_point"]
    # C1: physical scale, at the points of the windows k0 >= c1_from; also the sample rows (diagnostic path)
    c1, c1s = [], []
    sm = {x: samples(arms[x]) for x in ("old", "new", "base", "ctrl")}
    for k0, p, v in pts:
        if k0 < a.c1_from:
            continue
        r = {k: ratio(v["pair"][k], v["ctrl"][k]) for k in M_KEYS}
        v["c1_ratio"] = r
        c1.append(max(r.values()))
        if all(s is not None for s in sm.values()):
            n = k0 + p
            for nm, col in S_COLS.items():
                dp = abs(sm["new"][n, col] - sm["old"][n, col])
                dc = abs(sm["ctrl"][n, col] - sm["base"][n, col])
                tol = max(a.c1 * dc, 2e-10 * abs(sm["old"][n, col]))      # CSV values have 10 digits
                c1s.append((k0, p, nm, dp, dc, dp <= tol))
    out["C1_max_ratio"] = max(c1) if c1 else None
    out["C1_points"] = len(c1)
    out["C1_samples"] = {"checked": len(c1s), "fails": [x for x in c1s if not x[5]],
                         "available": all(s is not None for s in sm.values())}
    out["C1"] = bool(c1) and max(c1) <= a.c1 and not out["C1_samples"]["fails"]
    # C2: round-off scale, geometric means of pair / N per window and over all windows (pair 0 counts as a pass)
    per_w, allv, zero = {}, [], 0
    for k0, p, v in pts:
        for k in C2_KEYS:
            rr = v["pair"][k] / v["N"][k]
            if rr == 0:
                zero += 1
                continue
            per_w.setdefault(k0, []).append(rr)
            allv.append(rr)
    out["C2_gm_window"] = {k0: gm(v) for k0, v in per_w.items()}
    out["C2_gm_all"] = gm(allv)
    out["C2_zero_values"] = zero
    out["C2"] = (all(g <= a.gm_window for g in out["C2_gm_window"].values())
                 and (out["C2_gm_all"] is None or out["C2_gm_all"] <= a.gm))
    out["pass"] = out["D"] and out["C1"] and out["C2"]
    return out


def show(r, a):
    print("== %s: cells %d, particle cells %d" % (r["name"], r["cells"], r["particle_cells"]))
    for k0, w in r["windows"].items():
        s = []
        for p, v in w["points"].items():
            if not v["det"]:
                s.append("p%d nd" % p)
                continue
            t = "p%d" % p
            if "c1_ratio" in v:
                t += " C1 %.1e" % max(v["c1_ratio"].values())
            if "pair" in v:
                t += " C2 %.1e" % max(v["pair"][k] / v["N"][k] for k in C2_KEYS)
            s.append(t)
        g = r.get("C2_gm_window", {}).get(k0)
        print("   k0 %3d: tries %s, determined %d%s; %s%s" % (
            k0, w["tries"], w["det_tries"], "" if w["a4_case_interval"] is None else
            ", A4 case in interval %d" % w["a4_case_interval"], "; ".join(s), "" if g is None else "; C2 GM %.2f" % g))
    q = r["pilot"]
    print("   pilot: windows without a point %s; null/control in C1 at most %s (bound %.0e); GM null2/null1 %s (bound "
          "%s), per window %s .. %s (bound %s) -> %s"
          % (q["windows_without_point"], "%.1e" % q["null_c1_max"] if q["null_c1_max"] is not None else "-",
             a.pilot_c1, "%.2f" % q["n2_n1_gm"] if q["n2_n1_gm"] else "-", a.pilot_gm,
             "%.2f" % min(g for g in q["n2_n1_gm_window"].values() if g) if q["n2_n1_gm_window"] else "-",
             "%.2f" % max(g for g in q["n2_n1_gm_window"].values() if g) if q["n2_n1_gm_window"] else "-",
             a.pilot_gm_window, "PASS" if q["pass"] else "STOP"))
    if "pass" in r:
        print("   D %s (windows without a point %s, A4 case windows %s); C1 %s (largest %s, bound %.0e, %d points; "
              "sample rows %d checked, %d fail); C2 %s (GM all %s, bound %.0f; largest window GM %s, bound %.0f; "
              "%d values 0)" % (r["D"], r["windows_without_point"], r["a4_case_windows"], r["C1"],
                                "%.2e" % r["C1_max_ratio"] if r["C1_max_ratio"] is not None else "-", a.c1,
                                r["C1_points"], r["C1_samples"]["checked"], len(r["C1_samples"]["fails"]), r["C2"],
                                "%.2f" % r["C2_gm_all"] if r["C2_gm_all"] is not None else "-", a.gm,
                                "%.2f" % max(r["C2_gm_window"].values()) if r["C2_gm_window"] else "-", a.gm_window,
                                r["C2_zero_values"]))
        print("   RESULT %s: %s" % (r["name"], "PASS" if r["pass"] else "FAIL"))


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


def t1(name, old_dir, new_dir, a):
    """T1 of PLAN 14a: operators in double precision, new tree against old tree, as equivalent fields in cu."""
    import pathlib
    od, nd = pathlib.Path(old_dir), pathlib.Path(new_dir)
    mo, mn = torch.load(od / "meta.pt", map_location="cpu"), torch.load(nd / "meta.pt", map_location="cpu")
    if mo["names"] != mn["names"] or mo["ck"] != mn["ck"]:
        raise SystemExit("k_eval --t1: terms or states differ: %s %s / %s %s" % (mo["names"], mo["ck"], mn["names"],
                                                                               mn["ck"]))
    cu = 1e-6 * mn["H_amp"]
    al = mn["alpha"]
    gp = mn["gamma"] / (1.0 + (al ** 2 if torch.is_tensor(al) else al * al))   # gamma' of each cell (or scalar)
    Ms = mn["Ms"]
    full = Ms >= mn["Ms_full"] * (1 - 1e-9) if torch.is_tensor(Ms) else None
    void = Ms <= mn["Ms_full"] * 1e-6 if torch.is_tensor(Ms) else None
    res = {"name": name, "cu_A_per_m": cu, "states": {}, "worst": {}}
    for k0 in mn["ck"]:
        qo = torch.load(od / ("state_%03d.pt" % k0), map_location="cpu")
        qn = torch.load(nd / ("state_%03d.pt" % k0), map_location="cpu")
        dt = qn["dt"]
        st = {}
        for k in qn["q"]:
            if k not in qo["q"]:
                continue                                     # e.g. the B1 correction without last_h (T6)
            d = (qn["q"][k].double() - qo["q"][k].double()).norm(dim=-1)
            if k.startswith("h"):
                f = d                                        # field [A/m]
            elif k in ("rhs",) or k.startswith("stage"):
                f = d / gp                                   # right-hand side -> equivalent field
            else:
                f = d / (dt * gp)                            # increment, error estimate -> equivalent field
            v = {"rms_cu": float(f.pow(2).mean().sqrt()) / cu, "max_cu": float(f.max()) / cu}
            if k == "h_eff" and full is not None:
                cut = ~(full | void)
                for nm, msk in (("cut", cut), ("full", full), ("void", void)):
                    v["rms_cu_" + nm] = float(f[msk].pow(2).mean().sqrt()) / cu if bool(msk.any()) else 0.0
                    v["sq_sum_" + nm], v["n_" + nm] = float(f[msk].pow(2).sum()) / cu ** 2, int(msk.sum())
            st[k] = v
        e_cu = 4e-7 * math.pi * mn["Ms_full"] * cu * mn["volume"]   # energy of 1 cu on the full magnetization [J]
        for k in qn["sc"]:
            if k.startswith("E") and k in qo["sc"]:
                st[k] = {"rms_cu": abs(qn["sc"][k] - qo["sc"][k]) / e_cu, "max_cu": abs(qn["sc"][k] - qo["sc"][k]) / e_cu}
        if qn["sc"].get("dt_try") != qo["sc"].get("dt_try"):
            raise SystemExit("k_eval --t1: the trees used another dt for the try at k0 %d" % k0)
        no, nn = qo["sc"]["norm"], qn["sc"]["norm"]
        cls = lambda x: 2 if x > 1.1 else (1 if x >= 0.5 else 0)
        st["norm"] = {"rel": abs(nn / no - 1.0) if no > 0 else (0.0 if nn == 0 else math.inf), "old": no, "new": nn,
                      "same_class": cls(no) == cls(nn), "used_last_new": qn.get("used_last")}
        res["states"][k0] = st
    fields = [k for k in res["states"][mn["ck"][0]] if k != "norm"]
    worst_rms = max(res["states"][k0][k]["rms_cu"] for k0 in mn["ck"] for k in fields)
    worst_max = max(res["states"][k0][k]["max_cu"] for k0 in mn["ck"] for k in fields)
    worst_norm = max(res["states"][k0]["norm"]["rel"] for k0 in mn["ck"])
    same = all(res["states"][k0]["norm"]["same_class"] for k0 in mn["ck"])
    res["worst"] = {"rms_cu": worst_rms, "max_cu": worst_max, "norm_rel": worst_norm, "same_class": same}
    res["pass"] = worst_rms <= a.t1_rms and worst_max <= a.t1_max and worst_norm <= a.t1_norm and same
    if "h_eff" in res["states"][mn["ck"][0]] and "sq_sum_cut" in res["states"][mn["ck"][0]]["h_eff"]:
        pooled = {}
        for nm in ("cut", "full", "void"):
            ss = sum(res["states"][k0]["h_eff"]["sq_sum_" + nm] for k0 in mn["ck"])
            nn_ = sum(res["states"][k0]["h_eff"]["n_" + nm] for k0 in mn["ck"])
            pooled[nm] = math.sqrt(ss / nn_) if nn_ else 0.0
        res["h_eff_pooled_rms_cu"] = pooled
    print("T1 %-10s worst over %d states: RMS %.2e cu, max %.2e cu (bounds %.0e, %.0e); norm %.1e (bound %.0e); same "
          "branch %s -> %s" % (name, len(mn["ck"]), worst_rms, worst_max, a.t1_rms, a.t1_max, worst_norm, a.t1_norm,
                               same, "PASS" if res["pass"] else "FAIL"))
    for k in fields:
        print("     %-28s RMS %.2e cu (largest state), max %.2e cu" % (
            k, max(res["states"][k0][k]["rms_cu"] for k0 in mn["ck"]),
            max(res["states"][k0][k]["max_cu"] for k0 in mn["ck"])))
    if "h_eff_pooled_rms_cu" in res:
        print("     h_eff pooled over the states: cut %.2e cu, full %.2e cu, void %.2e cu" % tuple(
            res["h_eff_pooled_rms_cu"][n] for n in ("cut", "full", "void")))
    return res


NOISE = ("err", "err_rms", "norm")


def t2(name, f_new, f_old, f_start, f_null, a):
    """T2 of PLAN 14a: single-precision errors of the new tree against those of the old tree and of the chain start."""
    R = {}
    for k, f in (("new", f_new), ("old", f_old), ("start", f_start), ("null", f_null)):
        R[k] = torch.load(f, map_location="cpu", weights_only=False) if f else None
    cu = 1e-6 * R["new"]["meta"]["H_amp"]
    ks = sorted(R["new"]["states"])

    def rms(r, q, fld):
        # a quantity of the B1 trees only (the correction without the stored stray field): in the other trees it is
        # the correction itself
        if q not in r["states"][ks[0]]["q"] and q.endswith("_noh1"):
            q = q[:-len("_noh1")]
        v = [r["states"][k0]["q"][q].get(fld) for k0 in ks]
        return math.sqrt(sum(x * x for x in v) / len(v)) if all(x is not None for x in v) else None

    out, ok = {}, True
    for q in R["new"]["states"][ks[0]]["q"]:
        ref_small = q.startswith("h") and max(R["new"]["states"][k0]["q"][q].get("ref", 1e300) for k0 in ks) < cu
        if ref_small:
            vn, vo = rms(R["new"], q, "abs") / cu, rms(R["old"], q, "abs") / cu
            r = {"abs_cu_new": vn, "abs_cu_old": vo, "pass": vn <= vo + 1e-3}
        elif q in NOISE:
            en, eo, es = rms(R["new"], q, "e"), rms(R["old"], q, "e"), rms(R["start"], q, "e")
            enull = rms(R["null"], q, "e") if R["null"] else 0.0
            bound = 2 * max(eo, enull)
            r = {"e_new": en, "e_old": eo, "e_null": enull, "e_start": es, "bound": bound,
                 "pass": en <= bound and en <= 2 * max(es, enull)}
        else:
            en, eo, es = rms(R["new"], q, "e"), rms(R["old"], q, "e"), rms(R["start"], q, "e")
            bn, bo, bs = rms(R["new"], q, "b"), rms(R["old"], q, "b"), rms(R["start"], q, "b")
            fl = 2.0 ** -24
            r = {"e_new": en, "e_old": eo, "e_start": es, "b_new": bn, "b_old": bo, "b_start": bs,
                 "pass": en <= 2 * eo and en <= 2 * es and bn <= 2 * max(bo, fl) and bn <= 2 * max(bs, fl)}
        out[q] = r
        ok &= r["pass"]
    print("T2 %-10s %s" % (name, "PASS" if ok else "FAIL"))
    for q, r in out.items():
        print("     %-28s %s %s" % (q, " ".join("%s %.2e" % (k, v) for k, v in r.items() if k != "pass"),
                                    "" if r["pass"] else "<- FAIL"))
    return {"name": name, "quantities": out, "pass": ok}


def t6(old_dir, new_dir, a):
    """T6 of PLAN 14a (B1): (a) the correction of the new tree without the stored stray field equals the correction
    of the old tree (T1 bounds); (b) the change of the effective field in control units: cut cells pooled over the
    states in the band, full and void cells at most --t1_rms."""
    import pathlib
    od, nd = pathlib.Path(old_dir), pathlib.Path(new_dir)
    mn = torch.load(nd / "meta.pt", map_location="cpu")
    cu = 1e-6 * mn["H_amp"]
    Ms = mn["Ms"]
    full = Ms >= mn["Ms_full"] * (1 - 1e-9)
    void = Ms <= mn["Ms_full"] * 1e-6
    cut = ~(full | void)
    wr, wm, ss = 0.0, 0.0, {"cut": [0.0, 0], "full": [0.0, 0], "void": [0.0, 0]}
    for k0 in mn["ck"]:
        qo = torch.load(od / ("state_%03d.pt" % k0), map_location="cpu")["q"]
        qn = torch.load(nd / ("state_%03d.pt" % k0), map_location="cpu")["q"]
        kn = [k for k in qn if k.endswith("_noh1")]
        if len(kn) != 1:
            raise SystemExit("k_eval --t6: the new tree has no correction without last_h")
        ko = kn[0][:-len("_noh1")]
        d = (qn[kn[0]].double() - qo[ko].double()).norm(dim=-1)
        wr, wm = max(wr, float(d.pow(2).mean().sqrt()) / cu), max(wm, float(d.max()) / cu)
        dh = (qn["h_eff"].double() - qo["h_eff"].double()).norm(dim=-1)
        for nm, msk in (("cut", cut), ("full", full), ("void", void)):
            ss[nm][0] += float(dh[msk].pow(2).sum())
            ss[nm][1] += int(msk.sum())
    pooled = {nm: (math.sqrt(v[0] / v[1]) / cu if v[1] else 0.0) for nm, v in ss.items()}
    ok_a = wr <= a.t1_rms and wm <= a.t1_max
    ok_b = a.t6_band[0] <= pooled["cut"] <= a.t6_band[1] and pooled["full"] <= a.t1_rms and pooled["void"] <= a.t1_rms
    print("T6 (a) correction without last_h against the old correction: RMS %.2e cu, max %.2e cu (bounds %.0e, %.0e) "
          "-> %s" % (wr, wm, a.t1_rms, a.t1_max, "PASS" if ok_a else "FAIL"))
    print("T6 (b) change of the effective field, pooled over %d states: cut cells %.2e cu (band %s), full %.2e cu, "
          "void %.2e cu (bound %.0e) -> %s" % (len(mn["ck"]), pooled["cut"], a.t6_band, pooled["full"], pooled["void"],
                                                a.t1_rms, "PASS" if ok_b else "FAIL"))
    return {"a_rms_cu": wr, "a_max_cu": wm, "pooled_cu": pooled, "pass": ok_a and ok_b}


def g_check(old, new, a):
    """G of PLAN 14a: traces of the speed runs (prof_step.py --trace), old against new, scaled by old against old."""
    R = {"old": [json.load(open(f)) for f in old], "new": [json.load(open(f)) for f in new]}
    tr = {k: [r["trace"] for r in v] for k, v in R.items()}
    n = min(len(t) for v in tr.values() for t in v)
    t_equal = all(abs(x[i][0] - tr["old"][0][i][0]) <= 1e-12 * abs(tr["old"][0][i][0]) + 1e-30
                  for v in tr.values() for x in v for i in range(n))

    def rms(a_, b_, j, diff):
        if diff:   # tries per sample interval
            da = [a_[i][j] - a_[i - 1][j] for i in range(1, n)]
            db = [b_[i][j] - b_[i - 1][j] for i in range(1, n)]
        else:
            da = [a_[i][j] for i in range(n)]
            db = [b_[i][j] for i in range(n)]
        return math.sqrt(sum((x - y) ** 2 for x, y in zip(da, db)) / len(da))

    out = {"samples": n, "t_equal": t_equal}
    for nm, j, diff in (("tries", 1, True), ("Mx", 2, False)):
        oo = [rms(tr["old"][i], tr["old"][k], j, diff) for i in range(len(old)) for k in range(i + 1, len(old))]
        on = [rms(x, y, j, diff) for x in tr["old"] for y in tr["new"]]
        s_oo = math.sqrt(sum(v * v for v in oo) / len(oo))
        s_on = math.sqrt(sum(v * v for v in on) / len(on))
        out[nm] = {"old_old": s_oo, "old_new": s_on, "pass": s_on <= a.g_factor * s_oo}
    out["pass"] = t_equal and all(out[k]["pass"] for k in ("tries", "Mx"))
    print("G: %d samples, t equal %s; tries per sample interval RMS old/new %.3g, old/old %.3g; <Ms m_x> RMS old/new "
          "%.3g, old/old %.3g (bound %.0f x) -> %s" % (n, t_equal, out["tries"]["old_new"], out["tries"]["old_old"],
                                                      out["Mx"]["old_new"], out["Mx"]["old_old"], a.g_factor,
                                                      "PASS" if out["pass"] else "FAIL"))
    return out


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref")
    ap.add_argument("--pair", nargs=3, metavar=("NAME", "OLD", "NEW"))
    ap.add_argument("--pilot", metavar="NAME")
    ap.add_argument("--base")
    ap.add_argument("--ctrl")
    ap.add_argument("--null", nargs="+", metavar="N", help="N1 N2 [N3]: reverse, rotated, kernel null runs")
    ap.add_argument("--align_a4", action="store_true")
    ap.add_argument("--a4_min_windows", type=int, default=12)
    ap.add_argument("--c1", type=float, default=1e-2)
    ap.add_argument("--c1_from", type=int, default=8)
    ap.add_argument("--gm", type=float, default=10.0)
    ap.add_argument("--gm_window", type=float, default=30.0)
    ap.add_argument("--pilot_c1", type=float, default=1e-3)
    ap.add_argument("--pilot_gm", type=float, nargs=2, default=[0.5, 2.0])
    ap.add_argument("--pilot_gm_window", type=float, nargs=2, default=[0.1, 10.0])
    ap.add_argument("--json")
    ap.add_argument("--a4_count")
    ap.add_argument("--t1", nargs=3, metavar=("NAME", "OLD_DIR", "NEW_DIR"), help="T1: double-precision folders")
    ap.add_argument("--t1_rms", type=float, default=1e-3)
    ap.add_argument("--t1_max", type=float, default=0.1)
    ap.add_argument("--t1_norm", type=float, default=1e-6)
    ap.add_argument("--t2", nargs="+", metavar="X", help="T2: NAME NEW OLD START [NULL] (files of s1_states --against)")
    ap.add_argument("--t6", nargs=2, metavar=("OLD_DIR", "NEW_DIR"), help="T6: double-precision folders (B1)")
    ap.add_argument("--g_old", nargs="+", help="G: json files of prof_step.py --trace of the old code")
    ap.add_argument("--g_new", nargs="+", help="G: json files of prof_step.py --trace of the new code")
    ap.add_argument("--g_factor", type=float, default=3.0)
    ap.add_argument("--t6_band", type=float, nargs=2, default=[0.1, 10.0])
    ap.add_argument("--same", nargs=2)
    a = ap.parse_args(argv)
    if a.same:
        sys.exit(0 if same(*a.same) else 1)
    out = {}
    if a.a4_count:
        out["a4"] = a4_count(a.a4_count)
    if a.t1:
        out["t1"] = t1(*a.t1, a)
    if a.t6:
        out["t6"] = t6(*a.t6, a)
    if a.g_old:
        out["g"] = g_check(a.g_old, a.g_new, a)
    if a.t2:
        out["t2"] = t2(a.t2[0], a.t2[1], a.t2[2], a.t2[3], a.t2[4] if len(a.t2) > 4 else None, a)
    if a.pair or a.pilot:
        if not (a.ref and a.base and a.ctrl and a.null):
            raise SystemExit("k_eval: --ref, --base, --ctrl and --null are necessary")
        cache = {}

        def get(p):                                      # one copy of a file that has two roles (e.g. OLD = BASE)
            if p not in cache:
                cache[p] = load(p)
            return cache[p]

        ref = get(a.ref)
        arms = {"base": get(a.base), "ctrl": get(a.ctrl), "n1": get(a.null[0]), "n2": get(a.null[1])}
        if len(a.null) > 2:
            arms["n3"] = get(a.null[2])
        if a.pair:
            arms["old"], arms["new"] = get(a.pair[1]), get(a.pair[2])
        r = evaluate(a.pair[0] if a.pair else a.pilot, ref, arms, a.align_a4, a)
        show(r, a)
        out["result"] = r
    if a.json:
        with open(a.json, "w") as f:
            json.dump(out, f, indent=1, default=str)
    ok = True
    if "result" in out:
        r = out["result"]
        ok &= r["pass"] if "pass" in r else r["pilot"]["pass"]
    for k in ("t1", "t2", "t6", "g"):
        if k in out:
            ok &= out[k]["pass"]
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main(sys.argv[1:])
