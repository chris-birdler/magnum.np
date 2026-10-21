#!/usr/bin/env python
"""
Dynamic series in q (PLAN section 9, PROTOKOLL 3.29): the readings of phase A (q 4) and phase B (q 8), fixed before
the data. Loss density P at 9 mT, d/l_ex 96, φ 0.55, dx 3 l_ex; the value of a run is ln P of the calm tail
(eval_conv.run_record, as in eval_dynamic.py).

Arms: Coarse (T3), Coarse m corrected at q 2 / 4 / 8 (T3Pc, T3Pc4, T3Pc8), Coarse m exact at q 4 (X3q4).
Pairs: only the seeds in which both runs of an effect are free of an event (rule of 3.28; --events: the JSON of
check_events.py with all these arms). The event numbers of each arm are always reported.

    N3(q 4) = ln P(T3Pc4) − ln P(X3q4)            reading (K4)
    Δ_24    = ln P(T3Pc4) − ln P(T3Pc)            step q 2 → 4 (information)
    2 Δ_24                                        step q 2 → smooth spheres with first order in 1/q: reading (K4)
    Δ_c(4)  = ln P(T3Pc4) − ln P(T3)              information
    Δ_48    = ln P(T3Pc8) − ln P(T3Pc4)           step q 4 → 8 (information)
    Δ_24 + 2 Δ_48                                 step q 2 → smooth spheres with phase B: reading (K4)
    Δ_48 − Δ_24 / 2                               check of first order: "consistent" if the interval includes 0

Outcomes of a reading: "within ±0.15" (the 95 % t-interval lies inside the band), "outside ±0.15" (it lies outside),
else "not determined"; with fewer than 3 pairs always "not determined". Before Δ_24 and Δ_48 are read, the overlap of
the relaxed states (init.pt) of the two arms is given per particle (common magnetic cells; design B had 0.997 … 1.000).

    python eval_qseries.py --runs ROOT1 ROOT2 ... --events results/events_Q.json --json results/qseries_dyn.json
"""
import argparse
import json
import math
import pathlib
import sys

import numpy as np
import torch

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "beta_bpeak"))
sys.path.insert(0, str(HERE / "fe"))
import eval_conv as EC  # noqa: E402
from eval_dynamic import stats, interval  # noqa: E402
from check_events import magnetic_ids  # noqa: E402
from anisotropy_model import box_edge  # noqa: E402
from units import units  # noqa: E402

ARMS = ("T3", "T3Pc", "T3Pc4", "X3q4", "T3Pc8")
BAND = 0.15
# (key, meaning, function of the ln P of one seed, arms of the pair, reading or information)
EFFECTS = [
    ("N3(q 4)", "correction against exact field at q 4: ln P(T3Pc4) − ln P(X3q4)",
     lambda L: L["T3Pc4"] - L["X3q4"], ("T3Pc4", "X3q4"), True),
    ("Δ_24", "step q 2 → 4: ln P(T3Pc4) − ln P(T3Pc)", lambda L: L["T3Pc4"] - L["T3Pc"], ("T3Pc4", "T3Pc"), False),
    ("2 Δ_24", "step q 2 → smooth spheres (first order in 1/q), phase A", lambda L: 2 * (L["T3Pc4"] - L["T3Pc"]),
     ("T3Pc4", "T3Pc"), True),
    ("Δ_c(4)", "correction at q 4 against Coarse: ln P(T3Pc4) − ln P(T3)", lambda L: L["T3Pc4"] - L["T3"],
     ("T3Pc4", "T3"), False),
    ("Δ_48", "step q 4 → 8: ln P(T3Pc8) − ln P(T3Pc4)", lambda L: L["T3Pc8"] - L["T3Pc4"], ("T3Pc8", "T3Pc4"), False),
    ("Δ_24 + 2 Δ_48", "step q 2 → smooth spheres (first order in 1/q), phase B",
     lambda L: (L["T3Pc4"] - L["T3Pc"]) + 2 * (L["T3Pc8"] - L["T3Pc4"]), ("T3Pc", "T3Pc4", "T3Pc8"), True),
    ("Δ_48 − Δ_24 / 2", "check of first order in 1/q",
     lambda L: (L["T3Pc8"] - L["T3Pc4"]) - 0.5 * (L["T3Pc4"] - L["T3Pc"]), ("T3Pc", "T3Pc4", "T3Pc8"), False),
]


def reading(st, ci):
    if ci is None or st["n"] < 3:
        return "not determined (n %d)" % st["n"]
    if -BAND <= ci[0] and ci[1] <= BAND:
        return "within ±0.15"
    if ci[0] > BAND or ci[1] < -BAND:
        return "outside ±0.15"
    return "not determined"


def folder(roots, name):
    return next((r / name for r in roots if (r / name / "DONE").exists()),
                next((r / name for r in roots if (r / name).exists()), roots[0] / name))


def state_overlap(roots, a, b, s, cache):
    """Mean of m_a · m_b over the common magnetic cells of each particle, relaxed states (init.pt); None if missing."""
    da, db = folder(roots, "%s_s%d" % (a, s)), folder(roots, "%s_s%d" % (b, s))
    if not ((da / "init.pt").exists() and (db / "init.pt").exists()):
        return None
    cfg = json.loads((da / "config.json").read_text())
    lex = units(cfg["Js"], cfg["A"])["l_ex"]
    a_lex, d_lex = box_edge(cfg["d_lex"], cfg["phi"], cfg["a_quant"])
    N = int(round(a_lex / cfg["dx_lex"]))
    for arm in (a, b):
        if (arm, N) not in cache:
            cache[(arm, N)] = magnetic_ids(arm, N, d_lex * lex, cfg["dx_lex"] * lex)
    ia, ib = cache[(a, N)], cache[(b, N)]
    ma = torch.load(str(da / "init.pt"), map_location="cpu", weights_only=False)["m"].double().numpy()
    mb = torch.load(str(db / "init.pt"), map_location="cpu", weights_only=False)["m"].double().numpy()
    dot = (ma * mb).sum(-1)
    return [float(dot[(ia == k) & (ib == k)].mean()) for k in range(4)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True, help="run roots, searched in this order")
    ap.add_argument("--seeds", type=int, nargs="+", default=list(range(941, 953)))
    ap.add_argument("--events", required=True, help="JSON of check_events.py with the arms " + " ".join(ARMS))
    ap.add_argument("--all_runs", action="store_true",
                    help="information: use all runs with an event record, also the runs with an event (PROTOKOLL 3.31)")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    roots = [pathlib.Path(r) for r in a.runs]
    ev = json.loads(pathlib.Path(a.events).read_text())
    lnP, event, status = {}, {}, {}
    for arm in ARMS:
        for s in a.seeds:
            name = "%s_s%d" % (arm, s)
            rec = EC.run_record(folder(roots, name))
            status[name] = rec["status"]
            event[name] = bool(ev[name]["changed"]) if name in ev else None
            use = event[name] is not None if a.all_runs else event[name] is False   # no event record: not used
            lnP[(arm, s)] = rec["lnP"] if use else float("nan")
    out = {"lnP": {"%s_s%d" % k: v for k, v in lnP.items()}, "event": event, "status": status}
    print("# Dynamic series in q: readings of PLAN 9 (%s)\n" % (
        "all runs, also the runs with an event; information only" if a.all_runs
        else "pairs without an event, rule of 3.28"))
    print("| arm | runs done | runs with an event | runs without an event record | geometric mean of P (W/cm³), "
          "%s |" % ("runs used" if a.all_runs else "event-free runs"))
    print("|---|---|---|---|---|")
    for arm in ARMS:
        names = ["%s_s%d" % (arm, s) for s in a.seeds]
        done = sum((folder(roots, n) / "DONE").exists() for n in names)
        ok = [lnP[(arm, s)] for s in a.seeds if np.isfinite(lnP[(arm, s)])]
        print("| %s | %d | %d | %d | %s |" % (arm, done, sum(event[n] is True for n in names),
                                             sum(event[n] is None for n in names),
                                             "%.3g (n %d)" % (math.exp(np.mean(ok)), len(ok)) if ok else "–"))
    print("\n## Overlap of the relaxed states (init.pt), per particle, common magnetic cells\n")
    print("| pair | seed | particle 0 / 1 / 2 / 3 |")
    print("|---|---|---|")
    cache, out["overlap"] = {}, {}
    for x, y in (("T3Pc4", "T3Pc"), ("T3Pc8", "T3Pc4"), ("T3Pc4", "X3q4")):
        for s in a.seeds:
            ov = state_overlap(roots, x, y, s, cache)
            if ov is not None:
                out["overlap"]["%s-%s_s%d" % (x, y, s)] = ov
                print("| %s − %s | %d | %s |" % (x, y, s, " / ".join("%.4f" % v for v in ov)))
    print("\n## Effects (ln P; mean ± SE; 95 %% t-interval; band ±%.2f)\n" % BAND)
    print("| effect | meaning | per seed %s | mean ± SE (n) | 95 %% interval | reading |" % " / ".join(map(str, a.seeds)))
    print("|---|---|---|---|---|---|")
    out["effects"] = {}
    for key, txt, f, arms, is_reading in EFFECTS:
        per = []
        for s in a.seeds:
            L = {arm: lnP[(arm, s)] for arm in ARMS}
            per.append(f(L) if all(np.isfinite(L[x]) for x in arms) else float("nan"))
        st = stats([v for v in per if np.isfinite(v)])
        ci = interval(st)
        if key.startswith("Δ_48 − "):
            r = "not determined (n %d)" % st["n"] if ci is None or st["n"] < 3 else (
                "consistent with first order" if ci[0] <= 0.0 <= ci[1] else "not consistent with first order")
        else:
            r = reading(st, ci) + ("" if is_reading else " (information)")
        ms = ("%+.4f ± %.4f (SE, n %d)" % (st["mean"], st["se"], st["n"]) if st["n"] >= 2 and np.isfinite(st["se"])
              else "n %d" % st["n"])
        print("| %s | %s | %s | %s | %s | %s |" % (key, txt, " / ".join("%+.4f" % v if np.isfinite(v) else "–"
                                                                       for v in per),
                                                 ms, "%+.4f … %+.4f" % tuple(ci) if ci else "–", r))
        out["effects"][key] = dict(st, per_seed=per, ci95=ci, reading=r)
    if a.all_runs:
        print("\nStart rule of phase B (PLAN 9.3): not applied here (it uses the event-free pairs).")
        out["start_phase_B"] = None
    else:
        n3, a2 = out["effects"]["N3(q 4)"]["reading"], out["effects"]["2 Δ_24"]["reading"]
        start_b = n3 == "within ±0.15" and a2 != "within ±0.15"
        print("\nStart rule of phase B (PLAN 9.3): N3(q 4) %s; 2 Δ_24 %s; phase B starts: %s." % (
            n3, a2, "yes" if start_b else "no"))
        out["start_phase_B"] = start_b
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
