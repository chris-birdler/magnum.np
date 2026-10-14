#!/usr/bin/env python
"""
Where the loss channels sit (request of Chris 2026-10-13, item 2; A 10 pJ/m, corner K1, K1A10_0 … 2; descriptive, no
test). Channels of the cell loss p (O2, PROTOKOLL 7.72, 7.83): f (p_1 = p − harmonics − rest, the cycle mean of
|m_1|²), harmonics 2 … 7 f, rest (motion above 7 f or not repeating).
1. Zones in 3D (acc.json "where", calm cycles): centre r < 0.5 R, middle, rim r ≥ 0.8 R (R = particle radius), the
   10 % of the cells where the loop-mean m changes most between cells (texture top 10 %), the 10 % of the cells with
   the highest loss (p top 10 %). Per channel and B: concentration = share of the channel in the zone / volume share
   of the zone (1 = spread evenly); mean ± SD over the seeds. The f column uses the coherent m_1 (acc_reduce R1).
2. Maps in the cut z = 0 (acc_slices.npz) for the seed with the median A lever (as slide 10): each channel divided by
   the mean total loss of the particle cells in the cut (white = that mean), log scale; negative rest values (Bessel
   excess, 7.83) are shown at the bottom of the scale.

    python eval_channel_maps.py > results/channel_maps_K1A10.md      # also .png and .json (deck)
"""
import json
import pathlib

import numpy as np

R = pathlib.Path(__file__).resolve().parent / "runs"
RES = pathlib.Path(__file__).resolve().parent / "results"
AMPS = [9, 20, 35, 50, 70, 100, 150]
RUNS = ["K1A10_%d" % k for k in range(3)]
RUNS20 = ["K1_%d_f1" % k for k in range(3)]       # A 20 pJ/m with the same protocol (C′ runs at f), added 2026-10-14
ZONES = (("r<0.5", "centre r < 0.5 R"), ("0.5<=r<0.8", "middle"), ("r>=0.8", "rim r ≥ 0.8 R"),
         ("texture top 10%", "texture top 10 %"), ("p top 10%", "p top 10 %"))
CH = (("p_share", "total"), ("p_fund_share", "f"), ("p_harm_share", "2 … 7 f"), ("p_rest_share", "rest"))


def stages(run):
    out = {}
    for s in json.loads((R / run / "acc.json").read_text())["stages"]:
        out[min(AMPS, key=lambda x: abs(x / s["B_peak_mT"] - 1))] = s
    return out


def map_seed():
    pair = json.loads((RES / "pair_K1A10.json").read_text())
    lev = [b20 - b10 for b20, b10 in zip(pair["beta"]["9 … 150"]["20"], pair["beta"]["9 … 150"]["10"])]
    return int(np.argsort(lev)[len(lev) // 2])


def channel_maps(run, stage):
    """(f, harmonics, rest) in the cut, each / mean total loss of the particle cells; mask of the particle cells."""
    z = np.load(R / run / "acc_slices.npz")
    mask = z["ids"] >= 0
    p = z[stage + "_p_rel"].astype(np.float64)                  # already / the mean total loss in the cut
    pm = float(z[stage + "_p_mean"])
    h = z[stage + "_harm_rel"].astype(np.float64) * float(z[stage + "_harm_mean"]) / pm
    r = z[stage + "_rest_rel"].astype(np.float64) * float(z[stage + "_rest_mean"]) / pm
    return p - h - r, h, r, mask


def main():
    S = {run: stages(run) for run in RUNS}
    L = ["# Where the loss channels sit (A 10 pJ/m, corner K1, K1A10_0 … 2)", "",
         "Concentration = share of the channel in the zone / volume share of the zone (1 = spread evenly), mean ± SD "
         "over 3 seeds, 3D, calm cycles. Total = all of the cell loss. At 9 mT the rest is ≈ 1 % of the loss and its "
         "zone shares are not meaningful (cells with a Bessel excess have a negative rest).", ""]
    res = {"zones": {}, "amps": AMPS}
    for zk, zl in ZONES:
        res["zones"][zk] = {"label": zl, "vol": float(np.mean([S[r][9]["where"][zk]["vol_share"] for r in RUNS]))}
        vol = np.mean([S[r][9]["where"][zk]["vol_share"] for r in RUNS])
        L += ["## %s (volume share %.2f)" % (zl, vol), "", "| B (mT) | " + " | ".join(c for _, c in CH) + " |",
              "|---|" + "---|" * len(CH)]
        for b in AMPS:
            cells = []
            for ck, _ in CH:
                v = np.array([S[r][b]["where"][zk][ck] / S[r][b]["where"][zk]["vol_share"] for r in RUNS])
                res["zones"][zk].setdefault(ck, []).append([float(v.mean()), float(v.std(ddof=1))])
                cells.append("%.2f ± %.2f" % (v.mean(), v.std(ddof=1)))
            L.append("| %g | %s |" % (b, " | ".join(cells)))
        L.append("")
    seed = map_seed()
    L.append("Maps: %s (the seed with the median A lever of the pair)." % RUNS[seed])
    res["map_run"] = RUNS[seed]
    res["map_run20"] = RUNS20[seed]
    # A 20 against A 10: concentration of the total loss and of each group per zone (both protocols at f)
    S20 = {run: stages(run) for run in RUNS20 if (R / run / "acc.json").exists()}
    if len(S20) == 3:
        L.append("")
        L += ["## A 20 pJ/m (K1_k_f1) against A 10 pJ/m (K1A10_k): concentration c, mean ± SD over 3 seeds", "",
              "| zone | B (mT) | total A 20 / A 10 | f A 20 / A 10 | 2 … 7 f A 20 / A 10 | rest A 20 / A 10 |",
              "|---|---|---|---|---|---|"]
        res["zones20"] = {}
        for zk, zl in ZONES:
            res["zones20"][zk] = {}
            for b in AMPS:
                cells = []
                for ck, _ in CH:
                    v20 = np.array([S20[r][b]["where"][zk][ck] / S20[r][b]["where"][zk]["vol_share"] for r in RUNS20])
                    v10 = np.array([S[r][b]["where"][zk][ck] / S[r][b]["where"][zk]["vol_share"] for r in RUNS])
                    res["zones20"][zk].setdefault(ck, []).append([float(v20.mean()), float(v20.std(ddof=1))])
                    cells.append("%.2f / %.2f" % (v20.mean(), v10.mean()))
                if b in (9, 50, 150):
                    L.append("| %s | %g | %s |" % (zl, b, " | ".join(cells)))
        print("\n".join(L[-(3 + 3 * len(ZONES)):]))
    (RES / "channel_maps_K1A10.json").write_text(json.dumps(res, indent=1))
    print("\n".join(L))
    plot(RUNS[seed], S[RUNS[seed]])


def plot(run, st):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    fig, axs = plt.subplots(3, 3, figsize=(9.6, 9.6))
    norm = LogNorm(1e-2, 1e2)
    for c, b in enumerate((9, 50, 150)):
        f, h, r, mask = channel_maps(run, st[b]["stage"])
        for row, (v, lab) in enumerate(((f, "f"), (h, "2 … 7 f"), (r, "rest"))):
            ax = axs[row, c]
            img = np.where(mask, np.clip(v, 1e-2, None), np.nan)
            im = ax.imshow(img.T, origin="lower", cmap="coolwarm", norm=norm)
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            if row == 0:
                ax.set_title("%d mT" % b, fontsize=12)
            if c == 0:
                ax.set_ylabel(lab, fontsize=12)
            ax.text(0.02, 0.02, "%.0f %%" % (100 * np.sum(v[mask]) / np.sum((f + h + r)[mask])), transform=ax.transAxes,
                    fontsize=9, color="black", va="bottom")
    cb = fig.colorbar(im, ax=axs, shrink=0.6, pad=0.02)
    cb.set_label("loss density / mean total loss density in the cut", fontsize=10)
    fig.savefig(RES / "channel_maps_K1A10.png", dpi=130, bbox_inches="tight")


if __name__ == "__main__":
    main()
