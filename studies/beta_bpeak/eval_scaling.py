#!/usr/bin/env python
"""
How the loss scales with A (request of Chris 2026-10-13, item 3; rules written before the result).

Exponent s(B) = ∂ ln P / ∂ ln(d/l_ex) at fixed K_eff (Q), stress, packing, per peak flux density B. At a fixed particle
size d, d/l_ex ∝ A^(−1/2), thus P ∝ A^(−s/2). Reading (fixed now; "equal" = the 90 % interval contains the value):
  s = −2   P ∝ A: the exchange energy density A/d² of textures on the particle scale (vortex, curling) sets the loss;
  s = −1   P ∝ √A: the energy of domain walls (√(A K) per area, area per volume ∝ 1/d);
  s = 0    P does not depend on A (anisotropy energy);
  s > 0    the loss rises when A falls, more than these energy scales give (a change of the state or of its mobility).
The A lever on β follows from s: β(A 20) − β(A 10) ≈ [s(150 mT) − s(9 mT)] · ln(300/212) / ln(150/9) (two ends).
Estimators:
  (a) design model (Gamma GLM of analyze_design, all α 0.1 design runs): mean over the design points of
      [ln μ(ln d/l_ex = ln 300) − ln μ(ln 212)] / ln(300/212) at each nominal B; 1000 bootstrap resamples (90 %).
  (b) the same with K_eff fixed at L_eff 30 / 19 / 12 l_ex (coded ln Q −1 / 0 / +1; stress and packing at the design
      points).
  (c) the pair K1 / K1A10 (eval_pair.py, n 3): ln[P(A 10)/P(A 20)] / ln(d/l_ex ratio), per seed, mean, 90 % CI (t).
Limits: coarse mesh dx 3 (all runs); the quadratic design model; the pair is one corner (L_eff 30).

    python eval_scaling.py > results/scaling_A.md      # also results/scaling_A.json (deck)
"""
import math
import pathlib

import numpy as np
from scipy import stats

import analyze
import analyze_design as AD
import eval_pair as EP

R = pathlib.Path(__file__).resolve().parent / "runs"
AMPS = list(AD.AMPS_MT)
DLN = math.log(300.0 / 212.0)
LEVELS = (-1.0, 0.0, 1.0)               # coded ln Q_eff: L_eff 30 / 19 / 12 l_ex (the K_eff levels of the deck, slide A)
N_BOOT = 1000


def L_of(q):
    lo, hi = AD.RANGES["lnQ"]
    return math.exp(-0.5 * (lo + (q + 1) / 2 * (hi - lo)))


def coded_lnQ(L):
    lo, hi = AD.RANGES["lnQ"]
    return 2.0 * (math.log(1.0 / L ** 2) - lo) / (hi - lo) - 1.0


def s_of_B(c, X, lnQ=None):
    """[7] exponent s at each nominal B from GLM coefficients c, mean over the design points X."""
    hi, lo = X.copy(), X.copy()
    j = AD.FACTORS.index("lnd")
    hi[:, j], lo[:, j] = 1.0, -1.0
    if lnQ is not None:
        q = AD.FACTORS.index("lnQ")
        hi[:, q] = lo[:, q] = lnQ
    return ((AD.terms(hi) - AD.terms(lo)) @ c.T).mean(0) / DLN


def main():
    analyze.STEADY_TAIL, analyze.STEADY_MIN = AD.BAL_EVENT, 2
    recs, _ = AD.load_runs(R, reruns=True)
    X = AD.coded(recs)
    T = AD.terms(X)
    W = np.exp(np.array([[r["lnw@%gmT" % m] for m in AMPS] for r in recs]))
    variants = [("all design points", None)] + [("L_eff %.0f l_ex" % L_of(q), q) for q in LEVELS]
    c0 = AD.glm_coefs(T, W)
    est = {v: s_of_B(c0, X, q) for v, q in variants}
    rng = np.random.default_rng(20261013)
    boots = {v: [] for v, _ in variants}
    for _ in range(N_BOOT):
        i = rng.integers(0, len(recs), len(recs))
        try:
            c = AD.glm_coefs(T[i], W[i])
        except Exception:                               # noqa: BLE001 - a non-converging resample is skipped
            continue
        for v, q in variants:
            boots[v].append(s_of_B(c, X, q))
    L = ["# Scaling of the loss with A: exponent s = ∂ ln P / ∂ ln(d/l_ex)", "",
         "At a fixed particle size P ∝ A^(−s/2). s = −2: P ∝ A (exchange energy of textures on the particle scale); "
         "s = −1: P ∝ √A (domain walls); s = 0: no A dependence. Design model: %d α 0.1 design runs, Gamma GLM, %d "
         "bootstrap resamples (90 %% interval). Pair: K1 / K1A10, n 3, 90 %% CI with t. Coarse mesh." %
         (len(recs), len(boots["all design points"])), "",
         "| B (mT) | " + " | ".join(v for v, _ in variants) + " | pair K1 / K1A10 (L_eff 30) |",
         "|---|" + "---|" * (len(variants) + 1)]
    rec = EP.records()
    import json
    dl = [(math.log(json.loads((EP.folder(10, k) / "config.json").read_text())["d_lex_used"]) -
           math.log(json.loads((EP.folder(20, k) / "config.json").read_text())["d_lex_used"])) for k in EP.SEEDS]
    pair = []
    for b in AMPS:
        v = np.array([(rec[10][k]["lnw@%gmT" % b] - rec[20][k]["lnw@%gmT" % b]) / dl[k] for k in EP.SEEDS])
        h = stats.t.ppf(0.95, len(v) - 1) * v.std(ddof=1) / math.sqrt(len(v))
        pair.append((v.mean(), v.mean() - h, v.mean() + h))
    for i, b in enumerate(AMPS):
        cells = []
        for v, _ in variants:
            bb = np.array(boots[v])[:, i]
            cells.append("%+.2f (%+.2f … %+.2f)" % (est[v][i], *np.percentile(bb, [5, 95])))
        cells.append("%+.2f (%+.2f … %+.2f)" % pair[i])
        L.append("| %g | %s |" % (b, " | ".join(cells)))
    L += ["", "Check of the relation to β (two ends 9 and 150 mT): [s(150) − s(9)] · ln(300/212) / ln(150/9) = A lever "
          "on β between the ends."]
    for v, _ in variants:
        L.append("- %s: %+.2f" % (v, -(est[v][-1] - est[v][0]) * DLN / math.log(150 / 9)))
    L.append("- pair: %+.2f" % (-(pair[-1][0] - pair[0][0]) * np.mean(dl) / math.log(150 / 9)))
    print("\n".join(L))
    import json as js
    out = {"amps": AMPS, "variants": {v: {"est": est[v].tolist(), "lo": np.percentile(np.array(boots[v]), 5, 0).tolist(),
                                          "hi": np.percentile(np.array(boots[v]), 95, 0).tolist()} for v, _ in variants},
           "pair": [list(x) for x in pair], "note": "s = d ln P / d ln(d/l_ex); A exponent a = -s/2"}
    (pathlib.Path(__file__).resolve().parent / "results" / "scaling_A.json").write_text(js.dumps(out, indent=1))


if __name__ == "__main__":
    main()
