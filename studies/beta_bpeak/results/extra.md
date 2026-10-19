/home/christoph/Documents/code/magnum.np-beta_bpeak/studies/beta_bpeak/eval_extra.py:69: OptimizeWarning: Covariance of the parameters could not be estimated
  (lc, ltau), _ = curve_fit(model, om, np.log(w), p0=(math.log(w.max() * 2), math.log(5e-9)), maxfev=20000)
**A. Ring-down fit** p(t) = p_inf + A exp(-t/τ) over the 2 periods of the hold (T = 33.3 ns)

| run | B | p(0)/p_cycle | τ [ns] | τ/T | p_inf/p_cycle | ΔM_par during the hold / M amplitude |
|---|---|---|---|---|---|---|
| PH_s921 | 9 mT | 0.61 | 14.4 | 0.43 | 0.036 | -0.001 |
| PH_s921 | 20 mT | 0.53 | 22.5 | 0.68 | -0.036 | -0.001 |
| PH_s921 | 151 mT | 0.52 | 12.5 | 0.37 | 0.005 | +0.001 |
| PH_s922 | 9 mT | 0.90 | 17.3 | 0.52 | 0.013 | -0.003 |
| PH_s922 | 20 mT | 0.88 | 8.5 | 0.25 | 0.013 | +0.006 |
| PH_s922 | 150 mT | 0.76 | 15.1 | 0.45 | 0.021 | +0.001 |

**B. Debye fit of the frequency test** w = c x/(1 + x²), x = ωτ (f/2, f, 2f = 15, 30, 60 MHz)

| seed | B | τ [ns] | f_relax = 1/(2πτ) [MHz] | residual rms of ln w |
|---|---|---|---|---|
| 921 | 9 mT | 0.00 | 6471583 | 0.434 |
| 922 | 9 mT | 1.52 | 105 | 0.755 |
| 923 | 9 mT | 0.00 | 542760 | 0.549 |
| 921 | 20 mT | 0.00 | 1088294 | 0.310 |
| 922 | 20 mT | 0.00 | 11861685 | 0.626 |
| 923 | 20 mT | 1.33 | 120 | 0.236 |

9 mT: τ = 0.51 ± 0.51 ns (mean ± SE, n = 3) → f_relax ≈ 314 MHz

20 mT: τ = 0.44 ± 0.44 ns (mean ± SE, n = 3) → f_relax ≈ 359 MHz

**C. Events per mesh** (kept cycles with balance outside 0.9 … 1.1; drops of the number of large-angle pairs by ≥ 20 % from one cycle to the next, all cycles)

| set | runs | kept cycles | cycles with balance outside 0.9 … 1.1 | n60 drops |
|---|---|---|---|---|
| mesh300 dx 1.5 | 12 | 36 | 8 (22 %) | 9 |
| mesh300 dx 3 | 12 | 36 | 8 (22 %) | 3 |
| baseline dx 3 | 8 | 280 | 30 (11 %) | 5 |
