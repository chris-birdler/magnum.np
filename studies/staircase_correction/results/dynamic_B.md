# Dynamic arms at d/l_ex 96, φ 0.55, 9 mT (PROTOKOLL 3.5, D1 and D2)

| arm | meaning | P (W/cm³) per seed 941 / 942 / 943 / 944 / 945 / 946 / 947 / 948 / 949 / 950 / 951 / 952 | n | geometric mean (W/cm³) | SD of ln P | status / Bloch points of the prepared state |
|---|---|---|---|---|---|---|
| T3 | dx 3, voxel model, start from D15 | 1.03 / 0.794 / 0.605 / 0.776 / 0.606 / 11.2 / 0.784 / 0.6 / 1.22 / 0.981 / 0.793 / 1.02 | 12 | 1.01 | 0.79 | settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ? |
| T3Pc | dx 3, pairwise correction with centroid term, q 2, p 3, start from D15 (main arm) | 1.75 / 1.07 / 0.611 / 1.56 / 0.598 / 1.13 / 0.951 / 0.592 / 6.23 / 3.9 / 0.852 / 2.29 | 12 | 1.34 | 0.75 | settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ? |
| D15 | dx 1.5, voxel model of dx 1.5 (the source of the transfer) | 1.88 / 1.39 / 0.623 / 1.16 / 0.599 / 1.74 / 1.24 / 0.606 / 20.3 / 2.42 / 4.7 / 1.26 | 12 | 1.66 | 0.99 | settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0 |
| T15g | dx 1.5, voxel geometry of dx 3 locked, start from D15 | 1.16 / 0.857 / 0.615 / 0.874 / 0.62 / 1.06 / 0.867 / 0.61 / 1.58 / 1.08 / 0.856 / 1.13 | 12 | 0.906 | 0.29 | settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ? |
| X3 | dx 3, exact stray field of the fine staircase (reference of the correction, 3.18), start from D15 | 1.75 / 1.08 / 0.611 / 1.57 / 0.598 / 1.13 / 0.955 / 0.593 / 6.28 / 3.9 / 0.851 / 2.31 | 12 | 1.34 | 0.75 | settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ? |

Design B. Geometry effect of the fine staircase at dx 1.5, Δ_g = ln P(D15) − ln P(T15g): +0.604 ± 0.224 (SE, n 12)
Resolution of m at the dx 3 geometry, ln P(T15g) − ln P(T3): -0.113 ± 0.205 (SE, n 12)

| corrected arm | Δ_c = ln P(arm) − ln P(T3) | d = Δ_c − Δ_g (per seed) | 95 % interval of mean d (t) | reading |
|---|---|---|---|---|
| T3Pc | +0.275 ± 0.282 (SE, n 12) | -0.329 ± 0.289 (SE, n 12) | -0.97 … +0.31 | not determined |

## Proofs and mean effects (paired by seed; ln P; 95 % interval of the mean, t distribution)

| effect | meaning | per seed 941 / 942 / 943 / 944 / 945 / 946 / 947 / 948 / 949 / 950 / 951 / 952 | mean ± SE (n) | 95 % interval | factor on P (mean) | reading |
|---|---|---|---|---|---|---|
| N1 | resolution of m at fixed morphology (dx 3 geometry): ln P(T15g) − ln P(T3) | +0.119 / +0.077 / +0.015 / +0.119 / +0.022 / -2.357 / +0.101 / +0.017 / +0.260 / +0.092 / +0.076 / +0.100 | -0.113 ± 0.2049 (SE, n 12) | -0.564 … +0.338 | 0.8931 | not determined |
| N3 | correction against the exact fine stray field (same m): ln P(T3Pc) − ln P(X3) | -0.0001 / -0.0003 / -0.0001 / -0.0066 / +0.0005 / +0.0012 / -0.0042 / -0.0007 / -0.0071 / +0.0015 / +0.0009 / -0.0085 | -0.0020 ± 0.0010 (SE, n 12) | -0.0043 … +0.0003 | 0.9980 | within ±0.15 |
| Δ_g | geometry effect of the fine staircase at dx 1.5: ln P(D15) − ln P(T15g) | +0.484 / +0.484 / +0.013 / +0.281 / -0.034 / +0.495 / +0.359 / -0.007 / +2.552 / +0.810 / +1.703 / +0.113 | +0.604 ± 0.2239 (SE, n 12) | +0.112 … +1.097 | 1.8302 | information |
| Δ_c | geometry effect of the correction at dx 3: ln P(T3Pc) − ln P(T3) | +0.528 / +0.303 / +0.009 / +0.696 / -0.013 / -2.296 / +0.194 / -0.012 / +1.630 / +1.381 / +0.071 / +0.810 | +0.275 ± 0.2818 (SE, n 12) | -0.345 … +0.895 | 1.3167 | information |
| m | resolution of m on the fine staircase: ln P(D15) − ln P(X3) | +0.074 / +0.257 / +0.019 / -0.302 / +0.002 / +0.436 / +0.263 / +0.022 / +1.175 / -0.477 / +1.709 / -0.605 | +0.214 ± 0.1900 (SE, n 12) | -0.204 … +0.633 | 1.2390 | information |
| total | corrected dx 3 against dx 1.5 with its own staircase: ln P(T3Pc) − ln P(D15) | -0.074 / -0.258 / -0.019 / +0.295 / -0.001 / -0.435 / -0.267 / -0.023 / -1.182 / +0.479 / -1.708 / +0.597 | -0.216 ± 0.1899 (SE, n 12) | -0.634 … +0.202 | 0.8055 | information |
