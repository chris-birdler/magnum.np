#!/usr/bin/env python
"""
Long runs at 60 MHz after the C′ history (decision Chris 2026-10-15; reading fixed in PROTOKOLL 7.103 before the runs).
Runs (jobs/long60.txt, pair_jobs.py --long60): cond_K1A10_k_f2_b9 (large particle, 9 mT) and cond_K1_k_f2_b20 / _b35
(small particle, 20 / 35 mT), k = 0 … 2: the C′ line of K1A10_k_f2 / K1_k_f2 with the C′ history up to the tested
amplitude, then exactly 40 cycles at it. Reference: the same stage of the C′ run.

Control of the history (fixed before the runs): the cycles of the earlier stages and the first n cycles of the tested
stage (n = the cycles of that stage in C′) must repeat C′ (7.59: the same start path gives the same cycles). Largest
relative w_loop difference: ≤ 1e-9 "bit-identical"; ≤ 1e-4 "same within the solver tolerance"; above, or other cycle
counts, "differs". A restart inside the run ("[resume]" or "[init] m from" in stdout.log) is flagged. A run with
"differs" or a restart still enters the readings; a group with this in at least 2 of 3 seeds is marked "history not
reproduced: the comparison with C′ is statistical only".
Completeness: a group reading needs all 3 runs of the group, the A-effect reading needs all 9 runs; else "not
determined (runs missing)".
Secondary value (reported, not part of the readings): the event cycles (|w_dis / w_loop − 1| > 0.2 or w_loop ≤ 0)
among the last 10 cycles and P of the calm cycles among them; the C′ value excludes event cycles (steady tail).
Per run, the tested (last) stage from cycle 2 on: c_settle and class from eval_atol60.windows (window of 3 cycles: RN
≤ 0.05, closure ≤ 0.001 and w within ±5 % of the mean of the last 10 cycles; "periodic" if c_settle exists, "stays
irregular" if the last 10 cycles have median RN > 0.10 or median closure > 0.002, else "marginal").
Readings per group (large 9 mT, small 20 mT, small 35 mT; 3 seeds):
  motion        "periodic" if periodic in ≥ 2 of 3 seeds; "stays irregular" if stays irregular in ≥ 2 of 3; else
                "not determined";
  loss value    P = mean of the last 10 cycles; "the C′ value holds" if within ±10 % of the C′ value in ≥ 2 of 3 seeds;
                "the C′ value is off" if outside ±10 % in ≥ 2 of 3; else "not determined".
Consequence for the A effect on β at 60 MHz (fixed now): per seed β over 9 … 150 mT (main-fit definition) of K1_k_f2
(small) and K1A10_k_f2 (large), with the C′ values and again with ln w of the tested stages replaced by the long-run mean
(last 10 cycles). A effect = β(small) − β(large), paired by seed, mean ± SE (n 3), 90 % CI with t(2). Reading: "the
sign change holds" if the corrected mean < 0 and the CI excludes 0; "the sign change does not hold" if the corrected
mean ≥ 0; else "not determined".

    python eval_long60.py > runs/long60_report.md      # also results/long60.json
"""
import json
import math
import pathlib

import numpy as np
from scipy import stats

import analyze
import analyze_design as AD
import eval_atol60 as EA

R = pathlib.Path(__file__).resolve().parent / "runs"
SEEDS = range(3)
K_D = 1.5 ** 2 / (2 * 4e-7 * math.pi) / 1e6          # J/cm³
F_HZ = 60e6
GROUPS = (("large particle, 9 mT", "K1A10_%d_f2", 9), ("small particle, 20 mT", "K1_%d_f2", 20),
          ("small particle, 35 mT", "K1_%d_f2", 35))
STAGE = {9: 0, 20: 1, 35: 2}


def motion_reading(classes):
    if sum(c == "periodic" for c in classes) >= 2:
        return "periodic"
    if sum(c == "stays irregular" for c in classes) >= 2:
        return "stays irregular"
    return "not determined"


def value_reading(p, ref):
    ok = [x is not None and r is not None for x, r in zip(p, ref)]
    hold = sum(o and abs(x / r - 1) <= 0.10 for o, x, r in zip(ok, p, ref))
    off = sum(o and abs(x / r - 1) > 0.10 for o, x, r in zip(ok, p, ref))
    return "the C′ value holds" if hold >= 2 else ("the C′ value is off" if off >= 2 else "not determined")


def sign_reading(a, n_need=3):
    a = np.array([x for x in a if x is not None and math.isfinite(x)])
    if len(a) < n_need:
        return "not determined (runs missing)", float("nan"), float("nan"), (float("nan"), float("nan"))
    m, se = a.mean(), a.std(ddof=1) / math.sqrt(len(a))
    t = stats.t.ppf(0.95, len(a) - 1)
    ci = (m - t * se, m + t * se)
    if m < 0 and ci[1] < 0:
        return "the sign change holds", m, se, ci
    if m >= 0:
        return "the sign change does not hold", m, se, ci
    return "not determined", m, se, ci


def stage_series(d, si):
    """Per-window RN and closure and per-cycle w of stage si from cycle 2 on (as eval_atol60.run_series)."""
    st = [x for x in json.loads((d / "summary.json").read_text())["stages"] if x.get("measure")][si]
    cyc = {c["cycle"]: c for c in st["cycles"]}
    a = np.loadtxt(sorted(d.glob("samples_*_amp*.csv"))[si], delimiter=",", comments="#")
    P = {c: a[a[:, 0] == c][1:, 8] for c in cyc}
    n = int(np.bincount([len(v) for v in P.values()]).argmax())
    rn, clo, w = [], [], []
    for c in sorted(cyc):
        if c < EA.FIRST:
            continue
        win = [c - 2, c - 1, c]
        if not all(x in P and len(P[x]) == n for x in win):
            continue
        X = np.array([P[x] for x in win])
        rn.append(float(np.abs(X - X.mean(0)).sum() / X.sum()))
        clo.append(float(np.mean([cyc[x]["closure"] for x in win])))
        w.append(cyc[c]["w_loop"])
    return rn, clo, w


def history_class(diff, same_counts, restart):
    if not same_counts or diff > 1e-4:
        c = "differs"
    else:
        c = "bit-identical" if diff <= 1e-9 else "same within the solver tolerance"
    return c + (", restart" if restart else "")


def history_check(d, ref, si):
    """The earlier stages and the first n cycles of the tested stage si against the C′ run ref (n = cycles of stage si
    in C′): (largest relative w_loop difference, same cycle counts of the earlier stages, restart flag)."""
    a = [x for x in json.loads((d / "summary.json").read_text())["stages"] if x.get("measure")]
    b = [x for x in json.loads((ref / "summary.json").read_text())["stages"] if x.get("measure")]
    diffs = [abs(ca["w_loop"] / cb["w_loop"] - 1) for sa, sb in zip(a[:si + 1], b[:si + 1])
             for ca, cb in zip(sa["cycles"], sb["cycles"])]
    same = all(len(sa["cycles"]) == len(sb["cycles"]) for sa, sb in zip(a[:si], b[:si]))
    log = (d / "stdout.log").read_text() if (d / "stdout.log").exists() else ""
    restart = "[resume]" in log or "[init] m from" in log
    return (max(diffs) if diffs else float("nan")), same, restart


def secondary(d, si):
    """Event cycles among the last 10 cycles of stage si and P (W/cm³) of the calm cycles among them."""
    st = [x for x in json.loads((d / "summary.json").read_text())["stages"] if x.get("measure")][si]
    last = st["cycles"][-EA.LAST:]
    ev = [c for c in last if c["w_loop"] <= 0 or abs(c["w_dis"] / c["w_loop"] - 1) > AD.BAL_EVENT]
    calm = [c["w_loop"] for c in last if c not in ev]
    return len(ev), (F_HZ * K_D * float(np.mean(calm)) if calm else float("nan"))


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    res = {"groups": {}, "runs": {}}
    L = ["# Long runs at 60 MHz after the C′ history (reading fixed in PROTOKOLL 7.103)", "",
         "| run | history control (largest rel. w_loop difference to C′) | class | c_settle | P last 10 cycles (W/cm³, "
         "mean ± SD) | C′ value (W/cm³) | P / C′ | events in the last 10 / P of their calm cycles / C′ |",
         "|---|---|---|---|---|---|---|---|"]
    lnw_long = {}
    for title, pat, mT in GROUPS:
        cls, pv, ref, hist = [], [], [], []
        for k in SEEDS:
            name = "cond_%s_b%d" % (pat % k, mT)
            d, c = R / name, R / (pat % k)
            if not (d / "DONE").exists():
                L.append("| %s | missing | | | | | | |" % name)
                continue
            si = STAGE[mT]
            rn, clo, w = stage_series(d, si)
            s, cl = EA.windows(rn, clo, w)
            p = F_HZ * K_D * np.array(w[-EA.LAST:])
            rows = analyze.load_run(c)[1]
            pref = F_HZ * K_D * rows[si]["w_loop"]
            hd, same, rs = history_check(d, c, si)
            hc = history_class(hd, same, rs)
            nev, pcalm = secondary(d, si)
            cls.append(cl)
            pv.append(float(p.mean()))
            ref.append(pref)
            hist.append(hc)
            lnw_long[(pat % k, mT)] = math.log(float(np.mean(w[-EA.LAST:])))
            res["runs"][name] = {"class": cl, "c_settle": s, "P_mean": float(p.mean()), "P_sd": float(p.std(ddof=1)),
                                 "P_ref": pref, "history": hc, "history_diff": hd, "events_last10": nev,
                                 "P_calm_last10": pcalm}
            L.append("| %s | %s (%.1e) | %s | %s | %.1f ± %.1f | %.1f | %.2f | %d / %.1f / %.2f |" % (
                name, hc, hd, cl, s if s is not None else "–", p.mean(), p.std(ddof=1), pref, p.mean() / pref, nev,
                pcalm, pcalm / pref))
        if len(cls) < len(SEEDS):
            mr = vr = "not determined (runs missing)"
        else:
            mr, vr = motion_reading(cls), value_reading(pv, ref)
        bad = sum(h.startswith("differs") or "restart" in h for h in hist)
        note = " History not reproduced: the comparison with C′ is statistical only." if bad >= 2 else ""
        res["groups"][title] = {"motion": mr, "value": vr, "history_note": note.strip()}
        L += ["", "**%s: motion %s; loss value: %s.**%s" % (title, mr, vr, note), ""]
    # consequence: A effect on β at 60 MHz
    A0, A1 = [], []
    rows = []
    complete = len(lnw_long) == len(SEEDS) * len(GROUPS)
    for k in SEEDS:
        rs, rl = AD.run_record(R / ("K1_%d_f2" % k)), AD.run_record(R / ("K1A10_%d_f2" % k))
        if rs is None or rl is None:
            A0.append(None)
            A1.append(None)
            continue
        ls = np.array([rs["lnw@%gmT" % m] for m in AD.AMPS_MT])
        ll = np.array([rl["lnw@%gmT" % m] for m in AD.AMPS_MT])
        bs0, bl0 = AD.per_run_outputs(ls)["beta_all"], AD.per_run_outputs(ll)["beta_all"]
        ls1, ll1 = ls.copy(), ll.copy()
        for m in (20, 35):
            if ("K1_%d_f2" % k, m) in lnw_long:
                ls1[AD.AMPS_MT.index(m)] = lnw_long[("K1_%d_f2" % k, m)]
        if ("K1A10_%d_f2" % k, 9) in lnw_long:
            ll1[0] = lnw_long[("K1A10_%d_f2" % k, 9)]
        bs1, bl1 = AD.per_run_outputs(ls1)["beta_all"], AD.per_run_outputs(ll1)["beta_all"]
        A0.append(bs0 - bl0)
        A1.append(bs1 - bl1 if complete else None)
        rows.append("| %d | %.3f / %.3f | %.3f / %.3f | %+.3f | %+.3f |" % (k, bs0, bs1, bl0, bl1, bs0 - bl0, bs1 - bl1))
    r0, m0, se0, ci0 = sign_reading(A0)
    r1, m1, se1, ci1 = sign_reading(A1)
    res["A_effect"] = {"cprime": [m0, se0, list(ci0)], "corrected": [m1, se1, list(ci1)], "reading": r1}
    L += ["## A effect on β at 60 MHz (β over 9 … 150 mT, main-fit definition)", "",
          "| seed | β small (C′ / corrected) | β large (C′ / corrected) | A effect C′ | A effect corrected |",
          "|---|---|---|---|---|"] + rows
    corr = ("corrected: %+.3f ± %.3f (90 %% CI %+.3f … %+.3f), n 3" % (m1, se1, ci1[0], ci1[1]) if complete
            else "corrected: not computed (runs missing)")
    L += ["", "C′: %+.3f ± %.3f (90 %% CI %+.3f … %+.3f); %s." % (m0, se0, ci0[0], ci0[1], corr), "",
          "**Reading: %s.**" % r1]
    print("\n".join(L))
    (pathlib.Path(__file__).resolve().parent / "results" / "long60.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
