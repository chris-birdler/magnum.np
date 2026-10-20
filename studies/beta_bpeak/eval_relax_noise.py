#!/usr/bin/env python
"""
Noise floor of the energy in the converge relaxations (Chris 2026-10-20; local, no GPU).

At alpha 1 without a field the LLG dynamics only lowers E. A rise of E from one chunk to the next is therefore numerical
noise (fp32 rounding, integration error). Per relaxation (E per chunk from summary.json, K_d V units):
  rises      the positive dE per chunk -> median and 90 % quantile = noise level
  tail       median |dE| of the last quarter of the chunks (at least 4)
  q_share    share of the chunks from the first rise on with |dE| <= threshold (chance of a quiet chunk on the plateau)
Duplicates (the same history in nx/ad/pd/h16 runs) are counted once.

    python eval_relax_noise.py <folder with run folders (summary.json, config.json)> [--threshold 3e-10] > report.md
"""
import argparse
import json
import pathlib

import numpy as np


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
    a = ap.parse_args()
    th = a.threshold
    rows = [dict(r, **noise_stats(r.pop("E_hist"), th)) for r in relaxations(a.root)]
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
    print("\n".join(L))


if __name__ == "__main__":
    main()
