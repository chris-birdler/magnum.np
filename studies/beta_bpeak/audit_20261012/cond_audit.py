"""Task 7: conditioning, determinism and slow drifts (audit 2026-10-12).
(1) recompute the conditioning test (PLAN 4.5 D) independently; add a linear trend of ln P over the calm cycles 5 ... 29
    and the closure trend;
(2) slow drift in the design runs: slope of ln w_loop against the cycle number inside the calm tails (stages with >= 4
    kept cycles), per amplitude, mean ± SE over runs;
(3) determinism: compare the cond runs (earlier stages) and the reruns with their originals (bitwise).
Output: cond_audit.md
"""
import math

import numpy as np
from scipy import stats

import audit_lib as al
import design_data as dd


def tested(name):
    _, cfg, summ = al.read(name)
    return [s for s in summ["stages"] if s.get("measure")][-1], cfg


def main():
    L = ["# Audit: conditioning test, determinism, slow drifts", ""]
    L += ["## (1) Conditioning test recomputed (centre point, C0 … C5)", "",
          "| run | design P (W/cm³, calm tail of cycles 2 … 6) | late P (cycles 28, 29) | d = ln late/design | events (cycle) | "
          "trend of ln P over the calm cycles ≥ 5 (per 100 cycles, ± SE) | closure first / last calm cycle |",
          "|---|---|---|---|---|---|---|"]
    D = {9: [], 20: [], 35: []}
    T = {9: [], 20: [], 35: []}
    for b in (9, 20, 35):
        for s in range(6):
            st, cfg = tested("cond_C%d_b%d" % (s, b))
            pf = al.p_factor(cfg)
            cy = st["cycles"]
            tail = al.tail(cy[2:7])
            design = np.mean([c["w_loop"] for c in tail]) if len(tail) >= 2 else np.nan
            late = np.mean([c["w_loop"] for c in cy if c["cycle"] in (28, 29)])
            ev = [c["cycle"] for c in cy[2:] if al.is_event(c)]
            last_ev = max(ev) if ev else 1
            calm = [c for c in cy if c["cycle"] >= max(5, last_ev + 1)]
            k = np.array([c["cycle"] for c in calm], float)
            y = np.log([c["w_loop"] for c in calm])
            sl = stats.linregress(k, y)
            d = math.log(late / design)
            D[b].append(d)
            T[b].append(sl.slope * 100)
            L.append("| cond_C%d_b%d | %.3g | %.3g | %+.3f | %s | %+.3f ± %.3f | %.1e / %.1e |" % (
                s, b, design * pf, late * pf, d, ",".join(map(str, ev)) or "-", 100 * sl.slope, 100 * sl.stderr,
                calm[0]["closure"], calm[-1]["closure"]))
    L += ["", "| B (mT) | mean d ± SE (n 6) | 90 % t-CI | mean trend per 100 cycles ± SE | 95 % t-CI of the trend |", "|---|---|---|---|---|"]
    for b in (9, 20, 35):
        v = np.array(D[b])
        m, se = v.mean(), v.std(ddof=1) / math.sqrt(6)
        t90, t95 = stats.t.ppf(0.95, 5), stats.t.ppf(0.975, 5)
        tr = np.array(T[b])
        mt, st_ = tr.mean(), tr.std(ddof=1) / math.sqrt(6)
        L.append("| %d | %+.4f ± %.4f | %+.3f … %+.3f | %+.3f ± %.3f | %+.3f … %+.3f |" % (b, m, se, m - t90 * se, m + t90 * se,
                                                                                      mt, st_, mt - t95 * st_, mt + t95 * st_))
    # (2) slow drift in the design runs
    L += ["", "## (2) Drift inside the calm tails of the 116 design runs (main rule, reruns; stages with ≥ 4 kept cycles)", "",
          "| B (mT) | n stages | mean slope of ln w per cycle ± SE | stages with negative slope | median closure of the kept cycles |",
          "|---|---|---|---|---|"]
    for k, a in enumerate(al.AMPS):
        sl, clo = [], []
        for n in dd.design_names():
            _, _, summ = al.read(al.use_name(n))
            stg = [s for s in summ["stages"] if s.get("measure")]
            stg.sort(key=lambda s: s["b_target"])
            cyc = al.tail(al.measured(stg[k]))
            if len(cyc) >= 4:
                kk = np.arange(len(cyc), dtype=float)
                sl.append(np.polyfit(kk, np.log([c["w_loop"] for c in cyc]), 1)[0])
            clo += [c["closure"] for c in cyc]
        sl = np.array(sl)
        L.append("| %g | %d | %+.4f ± %.4f | %d | %.1e |" % (a, len(sl), sl.mean(), sl.std(ddof=1) / math.sqrt(len(sl)),
                                                            (sl < 0).sum(), np.median(clo)))
    # (3) determinism
    L += ["", "## (3) Determinism", ""]
    for s in range(6):
        _, _, c0 = al.read("C%d" % s)
        for b, nst in ((20, 1), (35, 2)):
            _, _, cc = al.read("cond_C%d_b%d" % (s, b))
            eq = all(a["w_loop"] == x["w_loop"] and a["w_dis"] == x["w_dis"] for i in range(nst)
                     for a, x in zip(c0["stages"][i]["cycles"], cc["stages"][i]["cycles"]))
            first = all(a["w_loop"] == x["w_loop"] for a, x in zip(c0["stages"][nst]["cycles"][:7], cc["stages"][nst]["cycles"][:7]))
            if not (eq and first):
                L.append("- C%d vs cond_C%d_b%d: NOT bit-identical" % (s, s, b))
    L.append("- cond runs: the earlier stages and the first 7 cycles of the tested stage are bit-identical to C0 … C5 "
             "(all 12 checks), unless listed above.")
    (al.HERE / "cond_audit.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
