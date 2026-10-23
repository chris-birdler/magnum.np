#!/usr/bin/env python
"""
Design B (PROTOKOLL 3.28) and the series in q (--set Q, PROTOKOLL 3.31): the effects of eval_dynamic.py or
eval_qseries.py with all runs and without the runs with an event, on one axis. Seed pairs with an event in one of the two arms are open circles; the means and 95 % intervals (t distribution)
of both sets are next to each other in each row. Information only: the exclusion rule was set after the data.

    python plot_events_compare.py --all results/dynamic_B.json --noevents results/dynamic_B_noevents.json \
        --plot results/proofs_B_compare.png
    python plot_events_compare.py --set Q --all results/dynamic_B.json results/qseries_dyn_all.json \
        --noevents results/dynamic_B_noevents.json results/qseries_dyn.json --plot results/proofs_Q_compare.png
With --set Q the rows are the effects that have the same meaning at q 2 and q 4 (N3, Δ_c), and the step Δ_24; N1 and
Δ_g do not depend on q, m and total need a fine reference with the roughness of q 4 (not computed).
"""
import argparse
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

SETS = {"B": (["N1", "N3", "Δ_g", "Δ_c", "m", "total"],
               [r"$N_1$", r"$N_3$", r"$\Delta_\mathrm{g}$", r"$\Delta_\mathrm{c}$", r"$m$", "total"]),
        "Q": (["N3", "N3(q 4)", "Δ_c", "Δ_c(4)", "Δ_24"],
              [r"$N_3(q = 2)$", r"$N_3(q = 4)$", r"$\Delta_\mathrm{c}(q = 2)$", r"$\Delta_\mathrm{c}(q = 4)$",
               r"$\Delta_{24}$"])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", nargs="+", required=True, help="JSON of eval_dynamic.py (or eval_qseries.py) with all runs")
    ap.add_argument("--noevents", nargs="+", required=True, help="JSON of eval_dynamic.py --exclude_events (or "
                    "eval_qseries.py), in the same order")
    ap.add_argument("--plot", required=True)
    ap.add_argument("--set", default="B", choices=sorted(SETS), help="B: design B (eval_dynamic.py); Q: series in q "
                    "(eval_qseries.py)")
    a = ap.parse_args()
    KEYS, LABELS = SETS[a.set]
    full, free = {}, {}
    for dst, files in ((full, a.all), (free, a.noevents)):
        for f in files:
            d = json.loads(pathlib.Path(f).read_text())
            dst.update(d.get("proofs", {}))                    # eval_dynamic.py
            dst.update(d.get("effects", {}))                   # eval_qseries.py
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.axvspan(-0.15, 0.15, color="0.9", zorder=0)
    ax.axvline(0.0, color="0.5", lw=0.8)
    for i, k in enumerate(KEYS):
        for v, w in zip(full[k]["per_seed"], free[k]["per_seed"]):
            kept = w == w                                          # NaN in the event-free set: pair with an event
            ax.plot([v], [i], "o", ms=4, zorder=2, color="0.45", mfc="0.45" if kept else "none")
        for r, dy, c in ((full[k], -0.17, "k"), (free[k], +0.17, "tab:blue")):
            if r["ci95"] is not None:
                ax.plot(r["ci95"], [i + dy] * 2, color=c, lw=2.0, zorder=3)
            ax.plot([r["mean"]], [i + dy], "D", color=c, ms=5, zorder=4)
    ax.set_yticks(range(len(KEYS)))
    ax.set_yticklabels(["%s  ($n = %d$ / $%d$)" % (t, full[k]["n"], free[k]["n"]) for t, k in zip(LABELS, KEYS)])
    ax.invert_yaxis()
    ax.set_xlabel(r"difference of $\ln P$")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=4, color="0.45", label="seed pair without event"),
                       Line2D([], [], marker="o", ls="", ms=4, color="0.45", mfc="none", label="seed pair with event"),
                       Line2D([], [], marker="D", ms=5, color="k", lw=2.0, label="mean, 95 % interval: all runs"),
                       Line2D([], [], marker="D", ms=5, color="tab:blue", lw=2.0,
                              label="mean, 95 % interval: without events"),
                       Patch(color="0.9", label="±0.15")],
              loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(a.plot, dpi=160)


if __name__ == "__main__":
    main()
