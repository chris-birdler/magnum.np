#!/usr/bin/env python
"""
Frequency test at the corner K3 (decision Chris 2026-10-15, PROTOKOLL 7.106; rules fixed before the runs).

Question: does the onset of the non-repeating motion at 9 mT move with K_eff? At K1 (K_eff 1.0 kJ/m³, small particle,
d/l_ex 212) the motion at 9 mT repeats up to 42 MHz (RN 0.037) and does not repeat at 60 MHz (RN 0.201; 7.99). The time
scale of motions that the anisotropy field restores scales with K_eff (7.94: f_K ∝ K); K3 has 2.76 × K_eff, thus the
prediction for an anisotropy-controlled onset is 2.76 × (42 … 60) = 116 … 166 MHz.

Runs (pair_jobs.py --k3freq [newell]): K3_k_fNNN (FD stray field) or K3N_k_fNNN (reference solver), k = 0, 1, 2,
f = 30, 42, 60, 85, 120, 170 MHz.
Valid run: DONE, converged virgin relaxation (summary init relax_converged), a value at all 7 amplitudes (steady tail
≥ 2 calm cycles). Other runs are listed and left out. A frequency enters the readings with ≥ 2 valid seeds.
Regularity at 9 mT (as 7.99, eval_sweep_regime.py): RN = non-repeating share of the box dissipation over the last 3
cycles of the 9 mT stage, mean over the valid seeds. Class: repeating RN ≤ 0.05; non-repeating RN > 0.10; marginal
otherwise. Onset = the lowest frequency with the class non-repeating.
Reading (pre-registered):
  shifted with K_eff     onset at ≥ 120 MHz (all lower frequencies repeating or marginal);
  partly shifted         onset at 85 MHz;
  not shifted            onset at ≤ 60 MHz (as K1);
  no onset up to 170 MHz no frequency non-repeating;
  not determined         fewer than 2 valid seeds at a frequency below the onset or at the onset.
Information (no reading): β over 9 … 150 mT per frequency (main-fit definition, eval_freq.beta), mean ± SE; the paired
change β(f) − β(30 MHz) with 90 % CI (t); RN over B; the excluded runs.

    python eval_k3freq.py [--prefix K3N] > runs/k3freq_report.md     (also results/k3freq_<prefix>.json)
"""
import json
import math
import pathlib
import sys

import numpy as np
from scipy import stats

import analyze_design as AD
import eval_freq as EF
import eval_jumps_time as EJ
from pair_jobs import K3_FREQ_MHZ

R = pathlib.Path(__file__).resolve().parent / "runs"
SEEDS = range(3)
AMPS = list(AD.AMPS_MT)
RN_REP, RN_NON = 0.05, 0.10
MIN_SEEDS = 2


def name(prefix, k, f):
    return "%s_%d_f%03d" % (prefix, k, round(f))


def rn_stages(d):
    """RN of the last 3 cycles per stage (eval_sweep_regime.run_stats), None where the samples are missing."""
    out = []
    for f in sorted(d.glob("samples_*_amp*.csv")):
        s = EJ.series(f)
        out.append(EJ.measures(*s)["RN"] if s is not None else None)
    return out


def run_values(d):
    """(valid, reason, lnw [7], RN per stage)."""
    if not (d / "DONE").exists():
        return False, "no DONE", None, []
    summ = json.loads((d / "summary.json").read_text())
    if not summ.get("init", {}).get("relax_converged"):
        return False, "relaxation not converged", None, []
    r = AD.run_record(d)
    lnw = np.array([r["lnw@%gmT" % b] for b in AMPS]) if r else np.full(len(AMPS), np.nan)
    if not np.all(np.isfinite(lnw)):
        return False, "no steady tail at %d amplitudes" % int(np.sum(~np.isfinite(lnw))), lnw, []
    return True, "", lnw, rn_stages(d)


def klass(rn):
    if rn is None or not np.isfinite(rn):
        return "–"
    return "repeating" if rn <= RN_REP else ("non-repeating" if rn > RN_NON else "marginal")


def reading(rn9, nvalid, freqs=K3_FREQ_MHZ):
    """rn9: {f: mean RN at 9 mT}; nvalid: {f: valid seeds}. Returns (key, text, onset)."""
    onset = None
    for f in freqs:
        if nvalid.get(f, 0) < MIN_SEEDS:
            return "undetermined", "not determined: %d valid seeds at %g MHz (need %d)" % (nvalid.get(f, 0), f,
                                                                                          MIN_SEEDS), None
        if klass(rn9[f]) == "non-repeating":
            onset = f
            break
    if onset is None:
        return "none", "no onset up to %g MHz" % freqs[-1], None
    if onset >= 120.0:
        return "shifted", "shifted with K_eff: onset at %g MHz (prediction 116 … 166 MHz)" % onset, onset
    if onset >= 84.0:
        return "partly", "partly shifted: onset at %g MHz (K1: 42 … 60 MHz; prediction 116 … 166 MHz)" % onset, onset
    return "not_shifted", "not shifted: onset at %g MHz, as K1 (42 … 60 MHz)" % onset, onset


def tci(v):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if len(v) < 2:
        return None
    se = v.std(ddof=1) / math.sqrt(len(v))
    h = stats.t.ppf(0.95, len(v) - 1) * se
    return v.mean(), se, v.mean() - h, v.mean() + h, len(v)


def main():
    prefix = sys.argv[sys.argv.index("--prefix") + 1] if "--prefix" in sys.argv else "K3"
    vals = {(f, k): run_values(R / name(prefix, k, f)) for f in K3_FREQ_MHZ for k in SEEDS}
    nvalid = {f: sum(vals[(f, k)][0] for k in SEEDS) for f in K3_FREQ_MHZ}
    rn9 = {f: float(np.mean([vals[(f, k)][3][0] for k in SEEDS if vals[(f, k)][0] and vals[(f, k)][3]
                             and vals[(f, k)][3][0] is not None])) if nvalid[f] else float("nan") for f in K3_FREQ_MHZ}
    beta = {(f, k): EF.beta(vals[(f, k)][2]) if vals[(f, k)][0] else float("nan") for f in K3_FREQ_MHZ for k in SEEDS}
    key, text, onset = reading(rn9, nvalid)
    L = ["# Frequency test at the corner K3 (%s; PROTOKOLL 7.106, rules fixed before the runs)" % prefix, "",
         "K3: K_eff 2.76 kJ/m³ (2.76 × K1), d/l_ex 212, K_σ 1.5 K_eff, φ 0.69, dx 3, 3 seeds. RN: non-repeating share of "
         "the dissipation over the last 3 cycles of a stage (mean over the valid seeds).", "",
         "| f (MHz) | valid seeds | RN at 9 mT | class | β 9 … 150 mT, mean ± SE (n) | β(f) − β(30 MHz), paired (90 % CI) |",
         "|---|---|---|---|---|---|"]
    for f in K3_FREQ_MHZ:
        b = tci([beta[(f, k)] for k in SEEDS])
        d = tci([beta[(f, k)] - beta[(30.0, k)] for k in SEEDS]) if f != 30.0 else None
        L.append("| %g | %d | %.3f | %s | %s | %s |" % (
            f, nvalid[f], rn9[f], klass(rn9[f]), "%.2f ± %.2f (%d)" % (b[0], b[1], b[4]) if b else "–",
            "%+.2f (%+.2f … %+.2f)" % (d[0], d[2], d[3]) if d else "–"))
    L += ["", "RN over B (mean over the valid seeds), amplitudes %s mT:" % " / ".join("%g" % b for b in AMPS), ""]
    for f in K3_FREQ_MHZ:
        per = [[vals[(f, k)][3][i] for k in SEEDS if vals[(f, k)][0] and len(vals[(f, k)][3]) > i
                and vals[(f, k)][3][i] is not None] for i in range(len(AMPS))]
        L.append("- %g MHz: %s" % (f, " / ".join("%.3f" % np.mean(p) if p else "–" for p in per)))
    bad = ["%s: %s" % (name(prefix, k, f), vals[(f, k)][1]) for f in K3_FREQ_MHZ for k in SEEDS if not vals[(f, k)][0]]
    L += ["", "Excluded runs: %s." % ("; ".join(bad) if bad else "none"), "",
          "Reference K1 (small particle, FD, 7.99): RN at 9 mT 0.037 at 42 MHz, 0.201 at 60 MHz.", "",
          "**Reading: %s.**" % text]
    out = {"prefix": prefix, "reading": key, "onset_MHz": onset, "rn9": {str(f): rn9[f] for f in K3_FREQ_MHZ},
           "valid": {str(f): nvalid[f] for f in K3_FREQ_MHZ},
           "beta": {"%g|%d" % (f, k): beta[(f, k)] for f in K3_FREQ_MHZ for k in SEEDS}}
    (pathlib.Path(__file__).resolve().parent / "results" / ("k3freq_%s.json" % prefix)).write_text(json.dumps(out, indent=1))
    print("\n".join(L))


if __name__ == "__main__":
    main()
