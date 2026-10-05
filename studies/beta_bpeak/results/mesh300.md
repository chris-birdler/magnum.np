**A. Loss per seed and mesh** (w = W/K_d per cycle; balance w_dis/w_loop)

| seed | B | w dx 1.5 | w dx 3 | Δ ln w (dx 3 − dx 1.5) | balance dx 1.5 / dx 3 |
|---|---|---|---|---|---|
| 921 | 9 mT | 1.483e-07 | 1.056e-07 | -0.339 | 1.08 / 1.70 |
| 921 | 50 mT | 5.968e-06 | 4.935e-06 | -0.190 | 1.00 / 1.07 |
| 921 | 100 mT | 1.88e-05 | 1.705e-05 | -0.098 | 1.02 / 2.35 |
| 922 | 9 mT | 1.893e-07 | 3.82e-07 | +0.702 | 1.56 / 0.99 |
| 922 | 50 mT | 5.013e-06 | 8.872e-06 | +0.571 | 1.02 / 1.11 |
| 922 | 100 mT | 1.849e-05 | 2.562e-05 | +0.326 | 1.03 / 1.02 |
| 923 | 9 mT | 2.079e-07 | 3.475e-07 | +0.514 | 1.04 / 0.97 |
| 923 | 50 mT | 4.691e-06 | 7.086e-06 | +0.412 | 1.00 / 1.00 |
| 923 | 100 mT | 1.752e-05 | 2.395e-05 | +0.313 | 0.99 / 1.00 |
| 924 | 9 mT | 2.454e-07 | 4.231e-07 | +0.545 | 1.04 / 0.98 |
| 924 | 50 mT | 5.42e-06 | 7.004e-06 | +0.256 | 1.42 / 1.03 |
| 924 | 100 mT | 1.912e-05 | 2.333e-05 | +0.199 | 0.97 / 1.01 |

| B | Δ ln w (paired, n = 4) | ratio w(dx 3)/w(dx 1.5) |
|---|---|---|
| 9 mT | +0.355 ± 0.235 | 1.43 |
| 50 mT | +0.262 ± 0.164 | 1.30 |
| 100 mT | +0.185 ± 0.099 | 1.20 |

**B. β of the segments, per mesh (mean ± SE over the seeds) and paired Δβ**

| segment | β dx 1.5 | β dx 3 | Δβ (dx 3 − dx 1.5), paired |
|---|---|---|---|
| 9→50 mT | 1.92 ± 0.08 | 1.87 ± 0.13 | -0.05 ± 0.06 |
| 50→100 mT | 1.81 ± 0.06 | 1.70 ± 0.06 | -0.11 ± 0.10 |
| 9→100 mT | 1.89 ± 0.04 | 1.82 ± 0.10 | -0.06 ± 0.06 |

**C. Rule of PLAN step 1b at 50 mT** (same mesh error at both ends of d/l_ex 212 … 300)

d/l_ex 300: +0.262 ± 0.164; d/l_ex 212 (7.10, clean): -0.06 ± 0.14; difference +0.322 ± 0.216; |difference| + 1.645 SE = 0.677 (rule: < 0.10)

**D. Protocol check, dx 3:** one job per amplitude from the virgin state (mesh300) vs amplitudes in sequence (baseline B300, same seeds)

| B | Δ ln w (mesh300 − baseline), paired n = 4 |
|---|---|
| 9 mT | -0.021 ± 0.016 |
| 50 mT | -0.068 ± 0.054 |
| 100 mT | -0.019 ± 0.011 |

d/l_ex 300 (1 µm), base point, α 0.1, f = 30 MHz, converged virgin state, cycles from 2 on (5 cycles per job). Errors: SE over the 4 seeds (paired differences).
