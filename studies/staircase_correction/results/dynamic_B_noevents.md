# Dynamic arms at d/l_ex 96, φ 0.55, 9 mT (PROTOKOLL 3.5, D1 and D2)

**Reading on the seed pairs without an event (decision of Chris, PROTOKOLL 3.28; rule fixed after the n-12 data):** the runs with an event (a change of the prepared state; results/events_B.json) are left out: D15_s949, D15_s951, T3_s946, T3_s949, T3_s951, T3Pc_s944, T3Pc_s947, T3Pc_s949, T3Pc_s950, X3_s944, X3_s947, X3_s949, X3_s950. The event numbers of each arm are a separate result.

| arm | meaning | P (W/cm³) per seed 941 / 942 / 943 / 944 / 945 / 946 / 947 / 948 / 949 / 950 / 951 / 952 | n | geometric mean (W/cm³) | SD of ln P | status / Bloch points of the prepared state |
|---|---|---|---|---|---|---|
| T3 | dx 3, voxel model, start from D15 | 1.03 / 0.794 / 0.605 / 0.776 / 0.606 / – / 0.784 / 0.6 / – / 0.981 / – / 1.02 | 9 | 0.782 | 0.22 | settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ? |
| T3Pc | dx 3, pairwise correction with centroid term, q 2, p 3, start from D15 (main arm) | 1.75 / 1.07 / 0.611 / – / 0.598 / 1.13 / – / 0.592 / – / – / 0.852 / 2.29 | 8 | 0.986 | 0.51 | settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ? |
| D15 | dx 1.5, voxel model of dx 1.5 (the source of the transfer) | 1.88 / 1.39 / 0.623 / 1.16 / 0.599 / 1.74 / 1.24 / 0.606 / – / 2.42 / – / 1.26 | 10 | 1.16 | 0.50 | settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0 |
| T15g | dx 1.5, voxel geometry of dx 3 locked, start from D15 | 1.16 / 0.857 / 0.615 / 0.874 / 0.62 / 1.06 / 0.867 / 0.61 / 1.58 / 1.08 / 0.856 / 1.13 | 12 | 0.906 | 0.29 | settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ? |
| X3 | dx 3, exact stray field of the fine staircase (reference of the correction, 3.18), start from D15 | 1.75 / 1.08 / 0.611 / – / 0.598 / 1.13 / – / 0.593 / – / – / 0.851 / 2.31 | 8 | 0.987 | 0.51 | settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ?; settled ? |

Design B. Geometry effect of the fine staircase at dx 1.5, Δ_g = ln P(D15) − ln P(T15g): +0.300 ± 0.088 (SE, n 10)
Resolution of m at the dx 3 geometry, ln P(T15g) − ln P(T3): +0.074 ± 0.014 (SE, n 9)

| corrected arm | Δ_c = ln P(arm) − ln P(T3) | d = Δ_c − Δ_g (per seed) | 95 % interval of mean d (t) | reading |
|---|---|---|---|---|
| T3Pc | +0.271 ± 0.140 (SE, n 6) | +0.096 ± 0.125 (SE, n 6) | -0.22 … +0.42 | not determined (event-free pairs) |

## Proofs and mean effects (paired by seed; ln P; 95 % interval of the mean, t distribution)

| effect | meaning | per seed 941 / 942 / 943 / 944 / 945 / 946 / 947 / 948 / 949 / 950 / 951 / 952 | mean ± SE (n) | 95 % interval | factor on P (mean) | reading |
|---|---|---|---|---|---|---|
| N1 | resolution of m at fixed morphology (dx 3 geometry): ln P(T15g) − ln P(T3) | +0.119 / +0.077 / +0.015 / +0.119 / +0.022 / – / +0.101 / +0.017 / – / +0.092 / – / +0.100 | +0.074 ± 0.0145 (SE, n 9) | +0.040 … +0.107 | 1.0765 | within ±0.15 (event-free pairs) |
| N3 | correction against the exact fine stray field (same m): ln P(T3Pc) − ln P(X3) | -0.0001 / -0.0003 / -0.0001 / – / +0.0005 / +0.0012 / – / -0.0007 / – / – / +0.0009 / -0.0085 | -0.0009 ± 0.0011 (SE, n 8) | -0.0035 … +0.0017 | 0.9991 | within ±0.15 (event-free pairs) |
| Δ_g | geometry effect of the fine staircase at dx 1.5: ln P(D15) − ln P(T15g) | +0.484 / +0.484 / +0.013 / +0.281 / -0.034 / +0.495 / +0.359 / -0.007 / – / +0.810 / – / +0.113 | +0.300 ± 0.0878 (SE, n 10) | +0.101 … +0.498 | 1.3497 | information |
| Δ_c | geometry effect of the correction at dx 3: ln P(T3Pc) − ln P(T3) | +0.528 / +0.303 / +0.009 / – / -0.013 / – / – / -0.012 / – / – / – / +0.810 | +0.271 ± 0.1399 (SE, n 6) | -0.089 … +0.631 | 1.3114 | information |
| m | resolution of m on the fine staircase: ln P(D15) − ln P(X3) | +0.074 / +0.257 / +0.019 / – / +0.002 / +0.436 / – / +0.022 / – / – / – / -0.605 | +0.029 ± 0.1218 (SE, n 7) | -0.269 … +0.327 | 1.0296 | information |
| total | corrected dx 3 against dx 1.5 with its own staircase: ln P(T3Pc) − ln P(D15) | -0.074 / -0.258 / -0.019 / – / -0.001 / -0.435 / – / -0.023 / – / – / – / +0.597 | -0.030 ± 0.1206 (SE, n 7) | -0.326 … +0.265 | 0.9701 | information |
