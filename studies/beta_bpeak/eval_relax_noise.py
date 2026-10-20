#!/usr/bin/env python
"""
Noise floor of the energy in the converge relaxations (Chris 2026-10-20; local, no GPU).

At alpha 1 without a field the LLG dynamics only lowers E. A rise of E from one chunk to the next is therefore numerical
noise (fp32 rounding, integration error). Per relaxation (E per chunk from summary.json, K_d V units):
  rises      the positive dE per chunk -> median and 90 % quantile = noise level
  tail       median |dE| of the last quarter of the chunks (at least 4)
  q_share    share of the chunks from the first rise on with |dE| <= threshold (chance of a quiet chunk on the plateau)
Duplicates (the same history in nx/ad/pd/h16 runs) are counted once.

Option 3 (trend criterion, offline check; Chris 2026-10-20): the slope of a straight line fitted to E over the last W
chunks, per chunk; the relaxation would stop at the first chunk with |slope| <= threshold. Per relaxation: the chunk of
that stop against the real stop (the recorded history ends at the real stop), and the slope over the last W chunks at the
real stop (the residual drift that the present criterion let through). --plateau FILE: a long flat history (lines of
run_loops with "E/K_d <value>") to test whether fp32 rounding gives a trend: slopes of all windows and of the whole file.

    python eval_relax_noise.py <folder with run folders (summary.json, config.json)> [--threshold 3e-10]
                               [--windows 6 8 12] [--plateau relax_sat.txt] > report.md
"""
import argparse
import json
import pathlib

import numpy as np


def window_slopes(E, W):
    """Slope (per chunk) and its standard error of a straight line fitted to each window of W chunks; entry j belongs to
    the window that ends at chunk j + W - 1."""
    E = np.asarray(E, float)
    x = np.arange(W) - (W - 1) / 2.0
    sl, se = [], []
    for j in range(E.size - W + 1):
        y = E[j:j + W]
        b = float((x * (y - y.mean())).sum() / (x * x).sum())
        r = y - y.mean() - b * x
        sl.append(b)
        se.append(float(np.sqrt((r * r).sum() / (W - 2) / (x * x).sum())) if W > 2 else float("nan"))
    return np.array(sl), np.array(se)


def trend_stop(E, W, th):
    """Index of the first chunk at which the trend criterion stops (|slope over the last W chunks| <= th), else None."""
    sl, _ = window_slopes(E, W)
    k = np.nonzero(np.abs(sl) <= th)[0]
    return int(k[0]) + W - 1 if k.size else None


def noise_stats(E, th):
    """Noise measures of one relaxation from its energies per chunk E (K_d V units)."""
    dE = np.diff(np.asarray(E, float))
    rises = dE[dE > 0]
    k0 = int(np.argmax(dE > 0)) if rises.size else dE.size
    plateau = np.abs(dE[k0:])
    tail = np.abs(dE[-max(4, dE.size // 4):])
    return {"n_rise": int(rises.size), "rise_med": float(np.median(rises)) if rises.size else float("nan"),
            "rise_p90": float(np.quantile(rises, 0.9)) if rises.size else float("nan"),
            "tail_med": float(np.median(tail)), "q_share": float((plateau <= th).mean()) if plateau.size else float("nan")}


def relaxations(root):
    seen = set()
    for d in sorted(pathlib.Path(root).iterdir()):
        try:
            s = json.loads((d / "summary.json").read_text())
            c = json.loads((d / "config.json").read_text())
        except (OSError, ValueError):
            continue
        i = s.get("init") or {}
        rel = [("virgin", i.get("relax_history") or [])]
        rel += [("relax2_p%d" % p["pass"], p.get("relax_history") or []) for p in i.get("predemag_passes") or []]
        for kind, h in rel:
            if len(h) < 3:
                continue
            key = (len(h), round(h[-1]["E_over_Kd"], 14))
            if key in seen:
                continue
            seen.add(key)
            nb = h[-1].get("n_bloch")
            yield {"run": d.name, "kind": kind, "it": h[-1]["it"], "E": h[-1]["E_over_Kd"], "n_cells": c.get("n_cells"),
                   "bp": sum(nb) if nb else None, "E_hist": [x["E_over_Kd"] for x in h]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--threshold", type=float, default=3e-10)
    ap.add_argument("--windows", type=int, nargs="+", default=[6, 8, 12])
    ap.add_argument("--plateau", default=None)
    a = ap.parse_args()
    th = a.threshold
    rows = []
    for r in relaxations(a.root):
        E = r.pop("E_hist")
        r.update(noise_stats(E, th))
        r["chunks"] = len(E)
        r["it_chunk"] = r["it"] / len(E)
        for W in a.windows:
            r["stop%d" % W] = trend_stop(E, W, th)
            sl, se = window_slopes(E, W) if len(E) >= W else (np.array([np.nan]), np.array([np.nan]))
            r["end%d" % W], r["end_se%d" % W] = float(sl[-1]), float(se[-1])
        rows.append(r)
    L = ["# Noise floor of the energy in the converge relaxations (eval_relax_noise.py)", "",
         "Threshold %.0e K_d V per chunk (PLAN 4.13 S7, D1 b). A rise of E is noise (alpha 1, no field)." % th, "",
         "| run | relaxation | iterations | E/(K_d V) | rises | rise median | rise 90 % | tail median \\|dE\\| | "
         "quiet share on plateau | cells | BPs end |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append("| %s | %s | %d | %.3e | %d | %.1e | %.1e | %.1e | %.2f | %.2e | %s |" % (
            r["run"], r["kind"], r["it"], r["E"], r["n_rise"], r["rise_med"], r["rise_p90"], r["tail_med"],
            r["q_share"], r["n_cells"], r["bp"]))
    ok = [r for r in rows if r["n_rise"] >= 3]
    rm = np.array([r["rise_med"] for r in ok])
    p9 = np.array([r["rise_p90"] for r in ok])
    qs = np.array([r["q_share"] for r in ok])
    L += ["", "Relaxations: %d (with at least 3 rises: %d)." % (len(rows), len(ok)),
          "Rise median: min %.1e, median %.1e, max %.1e; rise 90 %%: median %.1e, max %.1e." % (
              rm.min(), np.median(rm), rm.max(), np.median(p9), p9.max()),
          "Rise median / threshold: median %.2f, max %.2f; relaxations with the rise median above the threshold: %d."
          % (np.median(rm) / th, rm.max() / th, int((rm > th).sum())),
          "Quiet share on the plateau: min %.2f, median %.2f." % (qs.min(), np.median(qs))]
    for lab, key in (("|E|", "E"), ("cells", "n_cells"), ("iterations", "it")):
        x = np.log(np.array([abs(r[key]) for r in ok], float))
        L.append("Correlation of log rise median with log %s: %.2f." % (lab, np.corrcoef(x, np.log(rm))[0, 1]))
    bp = [(r["bp"], r["rise_med"]) for r in ok if r["bp"] is not None]
    wb, nb = [v for b, v in bp if b > 0], [v for b, v in bp if b == 0]
    if wb and nb:
        L.append("With BPs at the end: n %d, rise median %.1e; without: n %d, rise median %.1e." % (
            len(wb), np.median(wb), len(nb), np.median(nb)))
    W0 = a.windows[len(a.windows) // 2]
    L += ["", "## Option 3: trend criterion (straight line over the last W chunks, |slope| <= threshold)", "",
          "Iterations at the stop: real (present criterion) and with the trend criterion; '>' = not within the recorded "
          "history (it ends at the real stop). Slope at the real stop: over the last %d chunks, ± standard error." % W0, "",
          "| run | relaxation | real stop | " + " | ".join("trend W %d" % W for W in a.windows)
          + " | slope at the real stop (W %d) |" % W0, "|---|---|---|" + "---|" * len(a.windows) + "---|"]
    for r in rows:
        cells = []
        for W in a.windows:
            k = r["stop%d" % W]
            cells.append("%d" % ((k + 1) * r["it_chunk"]) if k is not None else "> %d" % r["it"])
        L.append("| %s | %s | %d | %s | %+.1e ± %.1e |" % (r["run"], r["kind"], r["it"], " | ".join(cells),
                                                          r["end%d" % W0], r["end_se%d" % W0]))
    for W in a.windows:
        st = [r for r in rows if r["stop%d" % W] is not None]
        ratio = [(r["stop%d" % W] + 1) * r["it_chunk"] / r["it"] for r in st]
        L.append("")
        L.append("W %d: stops within the history in %d of %d relaxations; trend stop / real stop: median %.2f, "
                 "min %.2f, max %.2f." % (W, len(st), len(rows), np.median(ratio) if ratio else np.nan,
                                          min(ratio) if ratio else np.nan, max(ratio) if ratio else np.nan))
    ends = np.array([r["end%d" % W0] for r in rows if np.isfinite(r["end%d" % W0])])
    L.append("Slope at the real stop (W %d): |slope| median %.1e, max %.1e; above the threshold in %d of %d." % (
        W0, np.median(np.abs(ends)), np.abs(ends).max(), int((np.abs(ends) > th).sum()), ends.size))
    if a.plateau:
        import re
        E = np.array([float(m.group(1)) for m in re.finditer(r"E/K_d (\S+)", pathlib.Path(a.plateau).read_text())])
        L += ["", "## Long flat history: %s (%d chunks)" % (pathlib.Path(a.plateau).name, E.size), ""]
        x = np.arange(E.size) - (E.size - 1) / 2.0
        b = float((x * (E - E.mean())).sum() / (x * x).sum())
        res = E - E.mean() - b * x
        se = float(np.sqrt((res * res).sum() / (E.size - 2) / (x * x).sum()))
        dE = np.diff(E)
        L.append("Whole file: slope %+.2e ± %.1e per chunk; |dE| median %.1e; SD of E %.1e." % (
            b, se, np.median(np.abs(dE)), E.std()))
        for W in a.windows:
            sl, _ = window_slopes(E, W)
            L.append("W %d: %d windows, |slope| median %.1e, 90 %% %.1e, max %.1e; share <= threshold %.2f; first "
                     "stop at chunk %s." % (W, sl.size, np.median(np.abs(sl)), np.quantile(np.abs(sl), 0.9),
                                            np.abs(sl).max(), float((np.abs(sl) <= th).mean()), trend_stop(E, W, th)))
    print("\n".join(L))


if __name__ == "__main__":
    main()
