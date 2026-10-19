# Audit: time curves and field data (116 design runs; kept cycles of the main rule unless stated)

## (A) Integrity of the loop integrals

W_loop and W_dis of every cycle (all stages, all cycles) recomputed from samples_*.csv: largest relative difference to summary.json 4.8e-08 (W_loop) and 8.5e-08 (W_dis) over the 116 run folders used.

## (B) Loop shape: permeability, loss angle, harmonics (median over 116 runs; quartiles)

| B (mT) | μ_r | sin δ | δ (deg) | \|M2\|/\|M1\| | \|M3\|/\|M1\| | π μ0 H1 M1 sin δ / W_loop |
|---|---|---|---|---|---|---|
| 9 | 6.88 (5.87 … 8.15) | 0.00585 (0.0038 … 0.00879) | 0.33 | 0.00046 (0.000279 … 0.000741) | 0.00032 (0.000137 … 0.00075) | 0.9984 … 1.0012 |
| 20 | 6.88 (5.87 … 8.15) | 0.00667 (0.00507 … 0.00973) | 0.38 | 0.000668 (0.000475 … 0.000892) | 0.000437 (0.000257 … 0.000751) | 0.9992 … 1.0007 |
| 35 | 6.87 (5.87 … 8.18) | 0.00741 (0.00591 … 0.00991) | 0.42 | 0.000683 (0.000421 … 0.00106) | 0.000526 (0.00035 … 0.000759) | 0.9990 … 1.0008 |
| 50 | 6.87 (5.86 … 8.12) | 0.00802 (0.00632 … 0.00981) | 0.46 | 0.000671 (0.000471 … 0.000926) | 0.000547 (0.000381 … 0.000908) | 0.9993 … 1.0008 |
| 70 | 6.87 (5.87 … 8.13) | 0.00832 (0.00673 … 0.0102) | 0.48 | 0.000708 (0.000457 … 0.00103) | 0.000533 (0.000369 … 0.000775) | 0.9976 … 1.0007 |
| 100 | 6.86 (5.88 … 8.13) | 0.0082 (0.00684 … 0.0105) | 0.47 | 0.000766 (0.000512 … 0.0011) | 0.000479 (0.000349 … 0.000699) | 0.9992 … 1.0007 |
| 150 | 6.86 (5.88 … 8.14) | 0.00793 (0.00648 … 0.0099) | 0.45 | 0.000702 (0.000499 … 0.00103) | 0.000553 (0.000365 … 0.000751) | 0.9995 … 1.0005 |

Decomposition β = 2 + s_δ − s_μ (s = slope over ln B, 9 … 150 mT, per run): s_δ mean +0.116 (sd 0.250), s_μ mean -0.0008 (sd 0.0015); remainder β − (2 + s_δ − s_μ): mean -0.0000, max |.| 0.0026 (the peak B differs slightly from the fundamental). Correlation of β with s_δ: r = 1.000; with s_μ: r = 0.243.

Factor effects (OLS HC3, range effect = 2 × coefficient of the coded factor; last column: lnd×lnQ interaction):

| quantity | lnQ | r_p | φ | lnd | lnd×lnQ | R² |
|---|---|---|---|---|---|---|
| β (full) | +0.134 ± 0.037 | +0.135 ± 0.061 | -0.062 ± 0.046 | -0.555 ± 0.047 | +0.409 ± 0.053 | 0.66 |
| s_δ | +0.136 ± 0.037 | +0.135 ± 0.061 | -0.064 ± 0.046 | -0.555 ± 0.047 | +0.410 ± 0.053 | 0.66 |
| −s_μ | -0.002 ± 0.000 | -0.001 ± 0.001 | +0.001 ± 0.001 | +0.000 ± 0.001 | -0.001 ± 0.001 | 0.26 |
| ln μ_r @9 | -0.020 ± 0.008 | -0.025 ± 0.010 | +0.744 ± 0.013 | -0.027 ± 0.011 | -0.024 ± 0.013 | 0.98 |
| ln μ_r @150 | -0.015 ± 0.008 | -0.022 ± 0.010 | +0.742 ± 0.012 | -0.027 ± 0.010 | -0.022 ± 0.012 | 0.98 |
| ln sin δ @9 | -0.331 ± 0.087 | -0.223 ± 0.127 | +0.909 ± 0.110 | +1.094 ± 0.109 | -0.706 ± 0.130 | 0.65 |
| ln sin δ @50 | -0.107 ± 0.090 | +0.040 ± 0.121 | +0.674 ± 0.109 | +0.327 ± 0.099 | -0.087 ± 0.132 | 0.33 |
| ln sin δ @150 | +0.131 ± 0.048 | +0.177 ± 0.075 | +0.728 ± 0.060 | -0.255 ± 0.067 | +0.487 ± 0.078 | 0.68 |
| |M3|/|M1| @150 | +0.000 ± 0.000 | +0.000 ± 0.000 | +0.001 ± 0.000 | -0.000 ± 0.000 | +0.000 ± 0.000 | 0.42 |

## (C) Loss per particle by stress-axis angle (runs with r_p > 0.5; per particle volume)

Ratio Q = (w of the 2 particles at 37°) / (w of the 2 particles at 79°) per stage; mean of ln Q ± SE over runs; slope of ln Q against r_p (OLS). Also the fundamental susceptibility ratio.

| B (mT) | n runs | ln Q mean ± SE | Q | slope d lnQ / d r_p ± SE | χ1 ratio 37°/79° (median) | check Σ particles / box |
|---|---|---|---|---|---|---|
| 9 | 100 | +0.370 ± 0.043 | 1.45 | -0.046 ± 0.062 | 1.019 | 1.000 … 1.000 |
| 20 | 100 | +0.393 ± 0.047 | 1.48 | -0.114 ± 0.067 | 1.020 | 1.000 … 1.000 |
| 35 | 100 | +0.456 ± 0.042 | 1.58 | -0.097 ± 0.060 | 1.019 | 1.000 … 1.000 |
| 50 | 100 | +0.560 ± 0.038 | 1.75 | -0.016 ± 0.056 | 1.019 | 1.000 … 1.000 |
| 70 | 100 | +0.593 ± 0.035 | 1.81 | -0.040 ± 0.051 | 1.020 | 1.000 … 1.000 |
| 100 | 100 | +0.600 ± 0.033 | 1.82 | -0.040 ± 0.048 | 1.020 | 1.000 … 1.000 |
| 150 | 100 | +0.607 ± 0.035 | 1.84 | +0.002 ± 0.051 | 1.020 | 1.000 … 1.000 |

## (D) Local loss measures (acc.json: all cycles ≥ 2 incl. events; phasemap.json: 16 phases of the last 2 cycles)

| B (mT) | R1 acc median | R1 acc, stages without events | R1 phasemap median | R1 acc − phasemap (median, n) | core share | shell share | texture-top-10 % share | p-top-10 % share |
|---|---|---|---|---|---|---|---|---|
| 9 | 0.861 | 0.906 (n 70) | 0.921 | -0.011 (116) | 0.384 | 0.224 | 0.610 | 0.694 |
| 20 | 0.717 | 0.745 (n 98) | – | – | 0.367 | 0.258 | 0.625 | 0.718 |
| 35 | 0.557 | 0.588 (n 91) | – | – | 0.305 | 0.343 | 0.654 | 0.749 |
| 50 | 0.496 | 0.503 (n 102) | 0.509 | -0.001 (116) | 0.274 | 0.387 | 0.693 | 0.760 |
| 70 | 0.471 | 0.472 (n 112) | – | – | 0.240 | 0.424 | 0.690 | 0.759 |
| 100 | 0.439 | 0.444 (n 109) | – | – | 0.208 | 0.454 | 0.657 | 0.749 |
| 150 | 0.400 | 0.400 (n 114) | 0.399 | +0.001 (116) | 0.176 | 0.499 | 0.550 | 0.726 |

Factor effects on R1 (phasemap, last 2 cycles; OLS HC3, range effects):

| B | lnQ | r_p | φ | lnd | lnd×lnQ | R² |
|---|---|---|---|---|---|---|
| 9 mT | +0.001 ± 0.023 | +0.019 ± 0.032 | -0.026 ± 0.029 | -0.082 ± 0.027 | +0.075 ± 0.031 | 0.12 |
| 50 mT | +0.051 ± 0.031 | +0.061 ± 0.041 | +0.049 ± 0.037 | -0.030 ± 0.031 | +0.027 ± 0.039 | 0.09 |
| 150 mT | -0.055 ± 0.011 | -0.051 ± 0.014 | +0.016 ± 0.013 | +0.023 ± 0.014 | -0.051 ± 0.018 | 0.40 |

## (E) Switched volume per half cycle (domains.json, last cycle at 9 / 50 / 150 mT) against β

- V_sw/V: median 0.65 / 3.31 / 9.81 % at 9 / 50 / 150 mT; exponent median 0.98 (0.68 … 1.13).
- Correlation β_full vs V_sw exponent: r = 0.910. Adding the standardised exponent to the factor model: lnd effect -0.555 → -0.130, interaction +0.409 → +0.060, exponent coefficient +0.219 ± 0.014 (per sd), R² 0.66 → 0.90.
- Caveat: the snapshots of the reruns stay from the originals (13 runs; 3 of them at 9 mT in event cycles).

## (F) Events in the measured cycles (cycles 2 … 6) of the 116 original design runs

- release of stored energy (w_dis > w_loop or w_loop ≤ 0): 173; uptake: 2.

## (G) Damping test α 0.01 vs α 0.1 (C0 … C7, same start state)

| rule | 9 mT | 20 mT | 35 mT | 50 mT | 70 mT | 100 mT | 150 mT | Δβ (7-amplitude fit) |
|---|---|---|---|---|---|---|---|---|
| all ≥1 | × 0.35 (ln -1.06 ± 0.35, n 8) | × 0.75 (ln -0.29 ± 0.20, n 8) | × 0.76 (ln -0.27 ± 0.17, n 8) | × 0.84 (ln -0.17 ± 0.13, n 8) | × 0.68 (ln -0.38 ± 0.14, n 8) | × 0.78 (ln -0.25 ± 0.09, n 8) | × 0.60 (ln -0.51 ± 0.06, n 8) | +0.17 ± 0.11 (n 8; max +0.83) |
| steady ≥1 | × 0.94 (ln -0.07 ± 0.65, n 2) | × 0.80 (ln -0.23 ± 0.61, n 3) | × 0.78 (ln -0.24 ± 0.19, n 7) | × 0.82 (ln -0.20 ± 0.14, n 7) | × 0.73 (ln -0.31 ± 0.14, n 7) | × 0.78 (ln -0.24 ± 0.08, n 7) | × 0.63 (ln -0.46 ± 0.06, n 7) | -0.21 ± 0.15 (n 8; max +0.46) |
| steady ≥2 | n 0 | × 0.80 (ln -0.23 ± 0.61, n 3) | × 0.79 (ln -0.23 ± 0.28, n 5) | × 0.73 (ln -0.31 ± 0.18, n 5) | × 0.75 (ln -0.29 ± 0.16, n 6) | × 0.81 (ln -0.21 ± 0.08, n 6) | × 0.62 (ln -0.48 ± 0.09, n 4) | -0.08 ± 0.27 (n 8; max +0.99) |
