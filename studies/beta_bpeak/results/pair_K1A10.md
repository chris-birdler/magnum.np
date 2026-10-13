# Pair K1 (A 20 pJ/m, d/l_ex 212) against K1A10 (A 10 pJ/m, d/l_ex 300)

Corner K1: L_eff 30 l_ex (K_eff 1.0 kJ/m³ at A 20 pJ/m, d = 1 µm), K_σ 1.5 K_eff, φ 0.62, α 0.1, f 30 MHz, seeds 2011 … 2013 (K1_0 … K1_2). Steady rule ≥ 2 calm cycles. Paired by seed; 90 % CI with t(n − 1); fewer than 3 seeds: no statement. The cell-by-cell axes of the pair differ (the Herzer cells scale with the box).

Runs used: K1_0_mc2, K1_1, K1_2, K1A10_0, K1A10_1, K1A10_2

## Loss density per B

| B (mT) | P at A 20 (W/cm³) | P at A 10 (W/cm³) | P(A 10) / P(A 20), paired (90 % CI) | settled A 20 / A 10 |
|---|---|---|---|---|
| 9 | 2.24 | 6.96 | × 3.11 (2.32 … 4.17), n 3 | 3 / 3 |
| 20 | 17.5 | 28.2 | × 1.61 (1.03 … 2.51), n 3 | 3 / 3 |
| 35 | 64.5 | 96 | × 1.49 (0.68 … 3.26), n 3 | 3 / 3 |
| 50 | 134 | 168 | × 1.26 (0.68 … 2.33), n 3 | 3 / 3 |
| 70 | 358 | 290 | × 0.81 (0.67 … 0.99), n 3 | 3 / 3 |
| 100 | 933 | 537 | × 0.58 (0.49 … 0.68), n 3 | 3 / 3 |
| 150 | 2.2e+03 | 1.13e+03 | × 0.51 (0.50 … 0.52), n 3 | 3 / 3 |

## β fits

| fit (mT) | β at A 20, per seed (mean / SD, n) | β at A 10, per seed | Δβ = A 10 − A 20, paired (seed offsets) |
|---|---|---|---|
| 9 … 150 | 2.45 / 0.12 (n 3) | 1.81 / 0.05 (n 3) | -0.64 ± 0.05 (90 % CI -0.74 … -0.55; 21 values, 3 seeds) |
| 20 … 150 (without 9 mT) | 2.44 / 0.20 (n 3) | 1.79 / 0.08 (n 3) | -0.65 ± 0.09 (90 % CI -0.80 … -0.50; 18 values, 3 seeds) |
| 9 / 20 / 35 | 2.48 / 0.08 (n 3) | 1.92 / 0.18 (n 3) | -0.56 ± 0.16 (90 % CI -0.89 … -0.23; 9 values, 3 seeds) |
| 20 / 35 / 50 | 2.23 / 0.46 (n 3) | 1.97 / 0.20 (n 3) | -0.26 ± 0.23 (90 % CI -0.73 … +0.21; 9 values, 3 seeds) |
| 35 / 50 / 70 | 2.47 / 0.40 (n 3) | 1.60 / 0.23 (n 3) | -0.87 ± 0.22 (90 % CI -1.32 … -0.42; 9 values, 3 seeds) |
| 50 / 70 / 100 | 2.80 / 0.41 (n 3) | 1.68 / 0.20 (n 3) | -1.12 ± 0.22 (90 % CI -1.57 … -0.67; 9 values, 3 seeds) |
| 70 / 100 / 150 | 2.38 / 0.18 (n 3) | 1.78 / 0.05 (n 3) | -0.60 ± 0.12 (90 % CI -0.84 … -0.35; 9 values, 3 seeds) |

## Measured against the prediction of the design model

Design model: Gamma GLM (quadratic in the 4 coded factors) on 116 α 0.1 design runs, evaluated at the coded factors of the pair runs (coded ln d/l_ex: A 20 -0.89, -0.89, -0.89, A 10 1.07, 1.07, 1.07); 1000 bootstrap resamples (90 % interval).

| output | measured (paired, 90 % CI) | design model (90 % interval) |
|---|---|---|
| Δβ 9 … 150 mT | -0.64 (-0.74 … -0.55) | -0.92 (-1.04 … -0.82) |
| Δβ 20 … 150 (without 9 mT) mT | -0.65 (-0.80 … -0.50) | -1.00 (-1.14 … -0.86) |
| Δβ 9 / 20 / 35 mT | -0.56 (-0.89 … -0.23) | -0.74 (-0.94 … -0.53) |
| Δβ 20 / 35 / 50 mT | -0.26 (-0.73 … +0.21) | -0.96 (-1.27 … -0.63) |
| Δβ 35 / 50 / 70 mT | -0.87 (-1.32 … -0.42) | -1.38 (-1.65 … -1.07) |
| Δβ 50 / 70 / 100 mT | -1.12 (-1.57 … -0.67) | -1.14 (-1.44 … -0.93) |
| Δβ 70 / 100 / 150 mT | -0.60 (-0.84 … -0.35) | -0.58 (-0.90 … -0.36) |
| P(A 10) / P(A 20) at 9 mT | × 3.11 (2.32 … 4.17) | × 5.61 (4.33 … 7.36) |
| P(A 10) / P(A 20) at 20 mT | × 1.61 (1.03 … 2.51) | × 3.27 (2.64 … 4.19) |
| P(A 10) / P(A 20) at 35 mT | × 1.49 (0.68 … 3.26) | × 2.05 (1.68 … 2.60) |
| P(A 10) / P(A 20) at 50 mT | × 1.26 (0.68 … 2.33) | × 1.34 (1.10 … 1.75) |
| P(A 10) / P(A 20) at 70 mT | × 0.81 (0.67 … 0.99) | × 0.79 (0.66 … 1.00) |
| P(A 10) / P(A 20) at 100 mT | × 0.58 (0.49 … 0.68) | × 0.61 (0.52 … 0.72) |
| P(A 10) / P(A 20) at 150 mT | × 0.51 (0.50 … 0.52) | × 0.50 (0.43 … 0.58) |

## Loss channels and location (mean over the seeds with a value, n in brackets)

| B (mT) | R1 A 20 | R1 A 10 | top-10 % share A 20 | top-10 % share A 10 | jump share A 20 | jump share A 10 |
|---|---|---|---|---|---|---|
| 9 | 0.97 (3) | 0.89 (3) | 0.53 (3) | 0.70 (3) | 0.01 (3) | 0.03 (3) |
| 20 | 0.80 (3) | 0.75 (3) | 0.59 (3) | 0.72 (3) | 0.06 (3) | 0.07 (3) |
| 35 | 0.72 (2) | 0.52 (3) | 0.62 (2) | 0.76 (3) | 0.07 (2) | 0.11 (2) |
| 50 | 0.56 (3) | 0.52 (3) | 0.68 (3) | 0.75 (3) | 0.12 (3) | 0.09 (2) |
| 70 | 0.50 (3) | 0.52 (3) | 0.70 (3) | 0.75 (3) | 0.13 (3) | 0.09 (3) |
| 100 | 0.42 (2) | 0.51 (3) | 0.72 (2) | 0.74 (3) | 0.12 (2) | 0.08 (3) |
| 150 | 0.37 (3) | 0.48 (3) | 0.71 (3) | 0.72 (3) | 0.10 (3) | 0.06 (3) |

R1: share of the cell loss at f (calm accumulator windows). Top-10 % share: share of the loss in the 10 % of the particle cells with the highest loss. Jump share: RF / 2 (energy above the smooth periodic curve, last 3 calm cycles).
