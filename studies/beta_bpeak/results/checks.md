**A. Frequency test at small B** (w per cycle; n from a fit of ln w over ln f per seed; mean ± SE, n = 3)

| B | w f/2 | w f | w 2f | exponent n (w ~ f^n) | balance f/2 / f / 2f |
|---|---|---|---|---|---|
| 9 mT | 7.5e-08 | 2.82e-07 | 3.89e-07 | 1.26 ± 0.19 | 1.01…1.05 / 0.98…1.40 / 1.00…1.15 |
| 20 mT | 4.28e-07 | 1.23e-06 | 1.78e-06 | 1.10 ± 0.15 | 1.01…1.01 / 1.01…2.07 / 0.99…1.04 |

| f | β(9 → 20 mT), mean ± SE |
|---|---|
| f/2 (15 MHz) | 2.26 ± 0.26 |
| f (30 MHz) | 1.94 ± 0.15 |
| 2f (60 MHz) | 1.98 ± 0.22 |

**B. Phase-resolved local response** (16 snapshots per cycle, last 3 cycles; shares of the measured dissipation p_dis)

| run | B | R1 (local fundamental) | R_np (non-periodic) | local phase lag mean / spread [°] |
|---|---|---|---|---|
| PH_s921 | 9 mT | 0.85 | 0.062 | 0.2 / 33.8 |
| PH_s921 | 20 mT | 0.81 | 0.009 | 1.8 / 32.0 |
| PH_s921 | 151 mT | 0.39 | 0.191 | -2.5 / 24.0 |
| PH_s922 | 9 mT | 0.95 | 0.000 | 0.5 / 58.7 |
| PH_s922 | 20 mT | 0.51 | 0.147 | -2.9 / 37.4 |
| PH_s922 | 150 mT | 0.44 | 0.049 | -2.9 / 20.4 |

The harmonics 2 … 7 (R_harm in phasemap.json) are not usable: motion faster than T/16 aliases into them and the weight (k ω)² amplifies it (R_harm = 3.6 … 44 > 1 at 9 … 20 mT).

**C. Ring-down** (H held at H max after the last cycle; p_dis relative to the mean dissipation of the cycles before; E_rd = dissipated energy in the 2 periods of the hold / loss per cycle)

| run | B | p(0) | p(0.01 T) | p(0.05 T) | p(0.1 T) | p(0.5 T) | p(2 T) | E_rd / w_cycle |
|---|---|---|---|---|---|---|---|---|
| PH_s921 | 9 mT | 0.609 | 0.572 | 0.443 | 0.358 | 0.172 | 0.027 | 0.266 |
| PH_s921 | 20 mT | 0.535 | 0.516 | 0.448 | 0.514 | 0.605 | 0.015 | 0.290 |
| PH_s921 | 151 mT | 0.519 | 0.417 | 0.464 | 0.503 | 0.110 | 0.004 | 0.196 |
| PH_s922 | 9 mT | 0.897 | 0.880 | 0.839 | 0.805 | 0.355 | 0.024 | 0.486 |
| PH_s922 | 20 mT | 0.878 | 0.824 | 1.060 | 0.960 | 0.160 | 0.007 | 0.412 |
| PH_s922 | 150 mT | 0.760 | 0.738 | 0.583 | 0.505 | 0.399 | 0.005 | 0.287 |

d/l_ex 300, dx 3, base point, α 0.1, converged virgin state. A: B300 gives f; all runs start from the same init.pt per seed and run 9 then 20 mT (7 cycles, kept from cycle 2).
