#!/usr/bin/env python
"""
Surface roughness with the converged virgin state (PROTOKOLL 7.29). d/l_ex 300, base point, seeds 921 ... 924,
one job per amplitude (9 / 50 / 100 mT), paired by seed.
  alpha 0.1 : A = staircase 10 nm (M300_dx3_*), B = staircase 5 nm (M300_dx15_*), C = smoothed 10 nm (R300f_*)
  alpha 0.01: staircase 10 nm (R300v01_*) vs smoothed 10 nm (R300f01_*), 9 cycles
Δ ln w (paired, mean ± SE) and β of the segments 9→50, 50→100, 9→100 mT.

    python eval_rough.py > results/rough.md

--hf: the high-field test of the staircase correction (PLAN 4.16, readings fixed before the data). Arms HFA (voxel),
HFC (fraction surface), HFD (fraction surface + correction q 8), all from the voxel state HFP_s<seed>. Per job: w =
mean w_loop of the steady tail of the corner evaluation (the cycles after the last cycle with |w_dis / w_loop − 1| >
BAL_EVENT 0.20, at least 2), b_peak measured. A job without a value (not done, relaxation not converged, no steady
tail) removes its seed from all three comparisons at that amplitude (one common seed set). Paired by seed: mean ± SE
and the 95 % interval of the t distribution for Δ ln w and Δβ, C − A, D − A, D − C; the class at 100 mT. Bloch
points: a seed whose three arms do not have the same Bloch flags (a Bloch point in the start state of the run, a Bloch
point at one sample or more) is reported; the results are given with and without such seeds.

    python eval_rough.py --hf [--runs DIR] > results/hf416.md
"""
import argparse
import json
import math
import pathlib

import numpy as np

import analyze
import analyze_design as AD
import eval_corners_sat as ECS

RUNS = pathlib.Path(__file__).resolve().parent / "runs"
SEEDS = (921, 922, 923, 924)
AMPS = (9, 50, 100)
SETS = {0.1: {"A staircase 10 nm": "M300_dx3_s%d_b%d", "B staircase 5 nm": "M300_dx15_s%d_b%d",
              "C smoothed 10 nm": "R300f_s%d_b%d"},
        0.01: {"staircase 10 nm": "R300v01_s%d_b%d", "smoothed 10 nm": "R300f01_s%d_b%d"}}
HF_ARMS = {"A": "HFA_s%d_b%d", "C": "HFC_s%d_b%d", "D": "HFD_s%d_b%d"}
HF_LABEL = {"A": "voxel", "C": "fraction", "D": "fraction + correction q 8"}
HF_PAIRS = (("C", "A"), ("D", "A"), ("D", "C"))
SEGMENTS = ((9, 50), (50, 100), (9, 100))


def one(name):
    _, rows = analyze.load_run(RUNS / name)
    return rows[0]["w_loop"], rows[0]["b_peak"], rows[0]["balance"]


def ms(v):
    v = np.asarray(v, float)
    return v.mean(), v.std(ddof=1) / math.sqrt(len(v))


def beta(W, s, lo, hi):
    return math.log(W[(s, hi)][0] / W[(s, lo)][0]) / math.log(W[(s, hi)][1] / W[(s, lo)][1])


def hf_job(d):
    """PLAN 4.16 reading of one job folder: {"w", "b", "bloch"} or {"fail": reason}. The steady tail of the corner
    evaluation (analyze.STEADY_TAIL = BAL_EVENT, STEADY_MIN 2); the module state of analyze is restored."""
    if not (d / "DONE").exists():
        return {"fail": "not done"}
    old = analyze.STEADY_TAIL, analyze.STEADY_MIN
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    try:
        _, rows = analyze.load_run(d)
    except analyze.UnconvergedRun:
        return {"fail": "relaxation not converged"}
    finally:
        analyze.STEADY_TAIL, analyze.STEADY_MIN = old
    if not rows:
        return {"fail": "no complete stage"}
    if rows[0]["missing"] or not rows[0]["w_loop"] > 0:
        return {"fail": "no steady tail of 2 cycles"}
    init = json.loads((d / "summary.json").read_text()).get("init", {})
    bloch = (sum(init.get("bloch_virgin") or [0]) > 0, ECS.stage_flags_b(d)[0])     # None: no count file
    return {"w": rows[0]["w_loop"], "b": rows[0]["b_peak"], "bloch": bloch}


def hf_class(ca, da, dc):
    """PLAN 4.16 class at 100 mT from the 95 % intervals (ECS.cstat tuples: mean, SE, SD, low, high, n) of Δ ln w for
    C − A, D − A and D − C."""
    if ca is None or da is None or dc is None:
        return 5, "not determined (fewer than 2 seeds)"
    if ca[3] <= 0:
        return 1, "excess not reproduced"
    below = dc[4] < 0
    if below and da[3] <= 0 <= da[4]:
        return 2, "the correction removes the excess"
    if not below and da[3] > 0:
        return 3, "the correction does not remove the excess"
    if below and da[3] > 0:
        return 4, "partly removed"
    return 5, "not determined"


def hf_tables(J, seeds):
    """Δ ln w per amplitude and Δβ per segment for HF_PAIRS over the seeds with a value in all three arms. Returns
    ({(pair, key): cstat}, {key: seeds used})."""
    out, used = {}, {}
    for mT in AMPS:
        ss = [s for s in seeds if all("w" in J[(a, s, mT)] for a in HF_ARMS)]
        used[mT] = ss
        for x, y in HF_PAIRS:
            out[((x, y), mT)] = ECS.cstat([math.log(J[(x, s, mT)]["w"] / J[(y, s, mT)]["w"]) for s in ss])
    for lo, hi in SEGMENTS:
        ss = [s for s in seeds if all("w" in J[(a, s, m)] for a in HF_ARMS for m in (lo, hi))]
        used[(lo, hi)] = ss

        def bt(a, s):
            return (math.log(J[(a, s, hi)]["w"] / J[(a, s, lo)]["w"]) /
                    math.log(J[(a, s, hi)]["b"] / J[(a, s, lo)]["b"]))
        for x, y in HF_PAIRS:
            out[((x, y), (lo, hi))] = ECS.cstat([bt(x, s) - bt(y, s) for s in ss])
    return out, used


def fmt_ci(c, f="%+.3f"):
    if c is None:
        return "–"
    return ("%s ± %s [%s, %s] n %d" % (f % c[0], (f % c[1]).lstrip("+"), f % c[3], f % c[4], c[5]))


def hf_report(runs=None):
    """The PLAN 4.16 report (markdown) and the record {"jobs", "fails", "bloch_unequal", "all", "no_bloch"}."""
    runs = pathlib.Path(runs) if runs else RUNS
    J = {(a, s, m): hf_job(runs / (p % (s, m))) for a, p in HF_ARMS.items() for s in SEEDS for m in AMPS}
    lines = ["**PLAN 4.16, high-field test** (arms: %s; paired by seed; mean ± SE [95 %% interval, t distribution] "
             "n seeds)\n" % ", ".join("%s %s" % (a, HF_LABEL[a]) for a in HF_ARMS)]
    fails = sorted((HF_ARMS[a] % (s, m), r["fail"]) for (a, s, m), r in J.items() if "fail" in r)
    lines.append("Jobs without a value: %d of %d%s\n" % (len(fails), len(J), "".join(
        "\n- %s: %s" % f for f in fails)))
    unequal = sorted({(s, m) for s in SEEDS for m in AMPS
                      if len({J[(a, s, m)].get("bloch") for a in HF_ARMS if "w" in J[(a, s, m)]}) > 1})
    lines.append("Seeds with unequal Bloch flags (start state, samples) between the arms: %s\n" % (
        ", ".join("s%d at %d mT %s" % (s, m, {a: J[(a, s, m)].get("bloch") for a in HF_ARMS}) for s, m in unequal)
        or "none"))
    rec = {"jobs": {HF_ARMS[a] % (s, m): r for (a, s, m), r in J.items()}, "fails": fails, "bloch_unequal": unequal}
    for tag, title in (("all", "all seeds"), ("no_bloch", "without the seeds with unequal Bloch flags")):
        Jt = J if tag == "all" else {k: (v if (k[1], k[2]) not in unequal else {"fail": "Bloch"}) for k, v in J.items()}
        out, used = hf_tables(Jt, SEEDS)
        lines.append("\n*%s*\n" % title)
        lines.append("| comparison | Δ ln w 9 mT | Δ ln w 50 mT | Δ ln w 100 mT |")
        lines.append("|---|---|---|---|")
        for x, y in HF_PAIRS:
            lines.append("| %s − %s | %s |" % (x, y, " | ".join(fmt_ci(out[((x, y), m)]) for m in AMPS)))
        lines.append("\n| comparison | Δβ 9→50 | Δβ 50→100 | Δβ 9→100 |")
        lines.append("|---|---|---|---|")
        for x, y in HF_PAIRS:
            lines.append("| %s − %s | %s |" % (x, y, " | ".join(fmt_ci(out[((x, y), g)], "%+.2f") for g in SEGMENTS)))
        k, text = hf_class(*(out[(p, 100)] for p in HF_PAIRS))
        lines.append("\nClass at 100 mT: %d, %s. Seeds used: %s\n" % (k, text, "; ".join(
            "%s: %s" % (("%d mT" % g) if isinstance(g, int) else ("%d→%d" % g), ",".join(map(str, used[g])) or "–")
            for g in list(AMPS) + list(SEGMENTS))))
        rec[tag] = {"class": k, "class_text": text, "seeds": {str(g): v for g, v in used.items()},
                    "stats": {"%s-%s %s" % (p[0], p[1], g): v for (p, g), v in out.items()}}
    return "\n".join(lines), rec


def main():
    for alpha, sets in SETS.items():
        if not all((RUNS / (p % (s, m)) / "DONE").exists() for p in sets.values() for s in SEEDS for m in AMPS):
            print("**α %g: not complete yet**\n" % alpha)
            continue
        D = {k: {(s, m): one(p % (s, m)) for s in SEEDS for m in AMPS} for k, p in sets.items()}
        keys = list(sets)
        ref = keys[0]
        print("**α %g** (reference: %s; paired over seeds %s, mean ± SE)\n" % (alpha, ref, ", ".join(map(str, SEEDS))))
        print("| comparison | Δ ln w 9 mT | Δ ln w 50 mT | Δ ln w 100 mT | Δβ 9→50 | Δβ 50→100 | Δβ 9→100 |")
        print("|---|---|---|---|---|---|---|")
        for k in keys[1:]:
            cells = []
            for m in AMPS:
                cells.append("%+.3f ± %.3f" % ms([math.log(D[k][(s, m)][0] / D[ref][(s, m)][0]) for s in SEEDS]))
            for lo, hi in ((9, 50), (50, 100), (9, 100)):
                cells.append("%+.2f ± %.2f" % ms([beta(D[k], s, lo, hi) - beta(D[ref], s, lo, hi) for s in SEEDS]))
            print("| %s − %s | %s |" % (k, ref, " | ".join(cells)))
        print("\n| set | β 9→50 | β 50→100 | β 9→100 | balance range |")
        print("|---|---|---|---|---|")
        for k in keys:
            b = [ms([beta(D[k], s, lo, hi) for s in SEEDS]) for lo, hi in ((9, 50), (50, 100), (9, 100))]
            bal = [D[k][(s, m)][2] for s in SEEDS for m in AMPS]
            print("| %s | %.2f ± %.2f | %.2f ± %.2f | %.2f ± %.2f | %.2f … %.2f |" %
                  (k, *b[0], *b[1], *b[2], min(bal), max(bal)))
        print()
    print("d/l_ex 300 (1 µm), base point, converged virgin state; α 0.1: 5 cycles, α 0.01: 9 cycles; cycles kept "
          "from 2 on. The α 0.01 runs start from the relaxed states of the α 0.1 runs.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--hf", action="store_true", help="the PLAN 4.16 report")
    ap.add_argument("--runs", default=None)
    ap.add_argument("--json", default=None, help="--hf: write the record to this file")
    a = ap.parse_args()
    if a.hf:
        text, rec = hf_report(a.runs)
        print(text)
        if a.json:
            pathlib.Path(a.json).write_text(json.dumps(rec, indent=1))
    else:
        main()
