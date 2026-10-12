# Same-state mesh test at K3 (A 20 pJ/m, d/l_ex 212)

Loss density P (W/cm³, geometric mean over the repeats), steady loops (≥ 2). n = repeats with a value.

| B (mT) | dx 1.5, own state | dx 3, fine state | dx 3, own state |
|---|---|---|---|
| 9 | 4.3 (n 6) | 2.6 (n 6) | 2.7 (n 6) |
| 50 | 213.1 (n 6) | 119.4 (n 6) | 148.4 (n 6) |
| 100 | 947.5 (n 6) | 732.5 (n 6) | 826.2 (n 6) |

Paired ratios per repeat (mean of ln ratio ± SE, as a factor; Welch-free paired t-test p):

| B (mT) | dx 3 fine state / dx 1.5 | dx 3 own state / dx 1.5 | dx 3 fine state / dx 3 own state |
|---|---|---|---|
| 9 | × 0.61 (ln -0.49 ± 0.16, n 6, p 0.027) | × 0.63 (ln -0.47 ± 0.09, n 6, p 0.003) | × 0.98 (ln -0.03 ± 0.09, n 6, p 0.783) |
| 50 | × 0.56 (ln -0.58 ± 0.14, n 6, p 0.009) | × 0.70 (ln -0.36 ± 0.13, n 6, p 0.043) | × 0.80 (ln -0.22 ± 0.16, n 6, p 0.239) |
| 100 | × 0.77 (ln -0.26 ± 0.08, n 6, p 0.026) | × 0.87 (ln -0.14 ± 0.09, n 6, p 0.189) | × 0.89 (ln -0.12 ± 0.10, n 6, p 0.298) |

β over 9 / 50 / 100 mT (power-law fit per repeat, mean ± SE):

- dx 1.5, own state: β 2.25 ± 0.04 (n 6)
- dx 3, fine state: β 2.32 ± 0.03 (n 6)
- dx 3, own state: β 2.37 ± 0.03 (n 6)

Decision rule (7.63): |mean ln ratio| ≤ 0.10 at all 3 amplitudes for dx 3 (fine state) / dx 1.5: not met ([np.False_, np.False_, np.False_])
