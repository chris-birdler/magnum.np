#!/usr/bin/env python
"""
Jumps in the dissipation (PROTOKOLL 7.21): direct test of the jump part of the loss.

Per kept cycle (from cycle 2 on) of each stage: the dissipation power p_dis(t) (256 samples per cycle,
samples_*.csv) is split into a smooth part (Fourier series of the drive phase up to harmonic KH; the
reversible linear response gives p ~ (dM/dt)^2, mainly harmonics 0 and 2) and a residual. A jump is
a sample with residual > NSIG robust sigma (1.4826 MAD). Outputs per stage:
  jumps per cycle, jump energy share = sum of the positive residual of the jump samples / total
  dissipation, non-smooth share = all positive residual / total dissipation (no threshold).
Sensitivity: KH 4 / 8 and NSIG 3 / 5. Note: a jump shorter than the sample interval (period / 256 =
130 ps at 30 MHz) is seen only as one sample or missed (sampling limit).

    python eval_jumps.py > results/jumps.md
"""
import math
import pathlib

import numpy as np

RUNS = pathlib.Path(__file__).resolve().parent / "runs"
GROUPS = {0.1: ["B300_s%d" % s for s in range(921, 929)], 0.01: ["A010c_s%d" % s for s in (904, 905, 906)]}
AMPS = (9, 20, 35, 50, 70, 100, 150)
SKIP = 2


def stage_files(run):
    return sorted((RUNS / run).glob("samples_*.csv"))


def analyse_cycle(t, p, KH, NSIG):
    n = len(p) - 1                                     # 257 samples: the last one closes the cycle
    t, p = t[:n], p[:n]
    ph = 2 * math.pi * np.arange(n) / n
    X = [np.ones(n)]
    for k in range(1, KH + 1):
        X += [np.cos(k * ph), np.sin(k * ph)]
    X = np.array(X).T
    c, *_ = np.linalg.lstsq(X, p, rcond=None)
    r = p - X @ c
    sig = 1.4826 * np.median(np.abs(r - np.median(r)))
    jump = r > NSIG * sig
    tot = p.sum()
    return int(jump.sum()), float(r[jump].sum() / tot), float(np.clip(r, 0, None).sum() / tot), float(tot)


def collect(KH, NSIG):
    out = {}
    for alpha, runs in GROUPS.items():
        for run in runs:
            for f in stage_files(run):
                a = np.loadtxt(f, delimiter=",", comments="#")
                cyc = a[:, 0].astype(int)
                B = None
                for k in sorted(set(cyc))[SKIP:]:
                    sel = cyc == k
                    Bk = 1000.0 * np.abs(a[sel, 7]).max()
                    nj, ej, ens, tot = analyse_cycle(a[sel, 1], a[sel, 8], KH, NSIG)
                    mT = min(AMPS, key=lambda x: abs(x - Bk))
                    out.setdefault((alpha, mT), {}).setdefault(run, []).append((nj, ej, ens, tot))
    return out


def ms(v):
    v = np.asarray(v, float)
    return v.mean(), (v.std(ddof=1) / math.sqrt(len(v)) if len(v) > 1 else float("nan"))


def main():
    base = collect(8, 4.0)
    print("**Jumps in the dissipation per stage** (KH = 8, threshold 4 σ; mean ± SE over the runs, each run "
          "averaged over its kept cycles)\n")
    print("| α | B_peak | runs | jumps per cycle | jump energy share | non-smooth share | jump energy per cycle "
          "[rel. to 9 mT] |")
    print("|---|---|---|---|---|---|---|")
    expo = {}
    for alpha in GROUPS:
        ref, EJ, BB = None, [], []
        for mT in AMPS:
            d = base.get((alpha, mT))
            if not d:
                continue
            nj = [np.mean([c[0] for c in v]) for v in d.values()]
            ej = [np.mean([c[1] for c in v]) for v in d.values()]
            en = [np.mean([c[2] for c in v]) for v in d.values()]
            eabs = np.mean([np.mean([c[1] * c[3] for c in v]) for v in d.values()])
            ref = ref or eabs
            EJ.append(eabs)
            BB.append(mT)
            print("| %g | %d mT | %d | %.1f ± %.1f | %.1f ± %.1f %% | %.1f ± %.1f %% | %.3g |" %
                  (alpha, mT, len(d), *ms(nj), *[100 * x for x in ms(ej)], *[100 * x for x in ms(en)], eabs / ref))
        ok = np.array(EJ) > 0
        expo[alpha] = np.polyfit(np.log(np.array(BB)[ok]), np.log(np.array(EJ)[ok]), 1)[0] if ok.sum() > 2 else np.nan
    print("\nGrowth of the jump energy per cycle with B_peak (fit of ln E_jump over ln B, group means): " +
          ", ".join("α %g: exponent %.2f" % kv for kv in expo.items()))
    print("\n**Sensitivity** (jump energy share, mean over the runs, α 0.1 / α 0.01)\n")
    print("| KH, threshold | " + " | ".join("%d mT" % m for m in AMPS) + " |")
    print("|---|" + "---|" * len(AMPS))
    for KH, NS in ((4, 4.0), (8, 3.0), (8, 4.0), (8, 5.0), (12, 4.0)):
        d = collect(KH, NS)
        cells = []
        for mT in AMPS:
            v = []
            for alpha in GROUPS:
                x = d.get((alpha, mT))
                v.append("%.1f" % (100 * np.mean([np.mean([c[1] for c in r]) for r in x.values()])) if x else "–")
            cells.append(" / ".join(v) + " %")
        print("| %d, %g σ | %s |" % (KH, NS, " | ".join(cells)))
    print("\nα 0.1: B300_s921 … s928 (7 amplitudes); α 0.01: A010c_s904 … s906 (9 / 20 / 50 / 100 / 150 mT). "
          "Cycles from 2 on. p_dis: LLG dissipation power of the unit cell per sample.")


if __name__ == "__main__":
    main()
