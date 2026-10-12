# Audit: β from power-law fits per run (116 design runs, steady loops ≥ 2, reruns)

## Full range 9 … 150 mT (7 amplitudes, realised B)

- β_full: mean 2.116, sd 0.250, SE 0.023, median 2.135, min 1.59 (S027), max 2.60 (S012); above 2: 78 of 116.
- Difference to the study output β_all (nominal-interpolated ln w, equal weights): max |Δ| 0.0092 (realised vs nominal B). To β_pow of analyze.features (weights ∝ number of kept cycles): mean -0.0023, max |Δ| 0.078.
- Fit quality: rms residual of ln P median 0.134 (quartiles 0.087 … 0.175), max |residual| median 0.227; SE of β_full from the residuals median 0.067.
- χ²/dof of the single power law against the cycle error of ln P: median 838 (quartiles 125 … 24963): the cycle error (≈ 1 %) is far smaller than the misfit; a single power law is not exact for any run.

## Curvature dβ/d ln B (quadratic fit in ln-ln, per run)

- mean -0.077 ± 0.023 (SE over runs, n 116); median -0.077; runs with |t| > 2.78 (95 %, dof 4): 24 of 116 (18 negative, 6 positive).

## Local β in windows (per-run power-law fits; mean ± SE over the 116 runs, sd, fraction > 2)

| window | mean ± SE | sd | median | fraction > 2 |
|---|---|---|---|---|
| 9-35 mT | 2.150 ± 0.027 | 0.288 | 2.143 | 0.72 |
| 20-50 mT | 2.151 ± 0.041 | 0.446 | 2.128 | 0.62 |
| 35-70 mT | 2.163 ± 0.046 | 0.493 | 2.147 | 0.62 |
| 50-100 mT | 2.067 ± 0.041 | 0.440 | 2.020 | 0.52 |
| 70-150 mT | 1.980 ± 0.034 | 0.366 | 1.984 | 0.49 |
| 9-50(4) mT | 2.155 ± 0.026 | 0.278 | 2.176 | 0.72 |
| 20-70(4) mT | 2.150 ± 0.037 | 0.398 | 2.126 | 0.61 |
| 35-100(4) mT | 2.113 ± 0.037 | 0.395 | 2.063 | 0.61 |
| 50-150(4) mT | 2.023 ± 0.032 | 0.346 | 1.963 | 0.46 |
| segment 9 → 20 mT | 2.163 ± 0.036 | 0.392 | 2.117 | 0.68 |
| segment 20 → 35 mT | 2.126 ± 0.045 | 0.483 | 2.108 | 0.63 |
| segment 35 → 50 mT | 2.196 ± 0.063 | 0.679 | 2.076 | 0.55 |
| segment 50 → 70 mT | 2.128 ± 0.060 | 0.643 | 2.104 | 0.59 |
| segment 70 → 100 mT | 2.011 ± 0.045 | 0.487 | 1.928 | 0.47 |
| segment 100 → 150 mT | 1.954 ± 0.042 | 0.449 | 1.887 | 0.38 |

Paired differences of local β within the runs (mean ± SE, t-test p, n 116):

- 20-50 − 9-35: +0.001 ± 0.039 (p 0.97)
- 35-70 − 9-35: +0.013 ± 0.050 (p 0.8)
- 70-150 − 35-70: -0.182 ± 0.054 (p 0.00091)
- 70-150 − 9-35: -0.169 ± 0.041 (p 7.5e-05)
- 50-100 − 20-50: -0.084 ± 0.056 (p 0.14)

## Broken power law (continuous, knot at the best inner amplitude)

- knot location (count of runs): 20 mT: 33, 35 mT: 22, 50 mT: 17, 70 mT: 31, 100 mT: 13
- β below / above the knot: 2.19 ± 0.04 / 2.00 ± 0.04 (mean ± SE); F (1, 3) > 10.1 (95 %): 35 of 116 runs.

## Dependence of the per-run β on the factors (OLS, HC3; coded factors ±1 = range ends; range effect = 2 × slope)

| output | lnQ (K_eff 1.0→6.2 kJ/m³) | r_p (0→3) | φ (0.55→0.69) | lnd (d/l_ex 212→300) | lnd×lnQ (interaction, coded) | R² |
|---|---|---|---|---|---|---|
| beta_full | +0.134 ± 0.037 | +0.135 ± 0.061 | -0.062 ± 0.046 | -0.555 ± 0.047 | +0.409 ± 0.053 | 0.66 |
| b9-35 | +0.130 ± 0.071 | +0.126 ± 0.101 | -0.068 ± 0.081 | -0.261 ± 0.088 | +0.449 ± 0.094 | 0.26 |
| b20-50 | +0.077 ± 0.107 | +0.079 ± 0.143 | -0.170 ± 0.137 | -0.756 ± 0.120 | +0.196 ± 0.173 | 0.28 |
| b35-70 | +0.061 ± 0.111 | +0.281 ± 0.136 | -0.115 ± 0.129 | -0.979 ± 0.125 | +0.421 ± 0.167 | 0.43 |
| b50-100 | +0.080 ± 0.111 | +0.068 ± 0.140 | +0.157 ± 0.127 | -0.694 ± 0.133 | +0.597 ± 0.193 | 0.37 |
| b70-150 | +0.323 ± 0.103 | +0.061 ± 0.130 | -0.009 ± 0.127 | -0.386 ± 0.119 | +0.390 ± 0.174 | 0.26 |
| curv | +0.049 ± 0.074 | -0.023 ± 0.098 | +0.059 ± 0.083 | -0.195 ± 0.088 | +0.001 ± 0.103 | 0.07 |

(The interaction column is 2 × the coefficient of x_lnQ·x_lnd: the change of the lnd range effect per unit of coded lnQ.)

## Mean local β against B by groups (mean over the runs of the group)

| group | n | 9-35 | 20-50 | 35-70 | 50-100 | 70-150 | full |
|---|---|---|---|---|---|---|---|
| L_eff 12 l_ex | 27 | 2.23 ± 0.06 | 2.18 ± 0.09 | 2.19 ± 0.10 | 2.13 ± 0.10 | 2.18 ± 0.09 | 2.20 ± 0.03 |
| L_eff 18 l_ex | 38 | 2.15 ± 0.04 | 2.15 ± 0.06 | 2.15 ± 0.07 | 2.02 ± 0.06 | 1.96 ± 0.05 | 2.10 ± 0.03 |
| L_eff 24 l_ex | 24 | 2.07 ± 0.07 | 2.05 ± 0.10 | 2.08 ± 0.10 | 2.01 ± 0.08 | 1.87 ± 0.06 | 2.03 ± 0.06 |
| L_eff 30 l_ex | 27 | 2.14 ± 0.05 | 2.21 ± 0.09 | 2.24 ± 0.10 | 2.12 ± 0.09 | 1.90 ± 0.07 | 2.13 ± 0.06 |
| A 9.7 … 12.5 pJ/m | 39 | 2.05 ± 0.04 | 1.92 ± 0.06 | 1.86 ± 0.05 | 1.83 ± 0.06 | 1.88 ± 0.06 | 1.94 ± 0.04 |
| A 12.5 … 15.8 | 38 | 2.16 ± 0.04 | 2.06 ± 0.07 | 2.06 ± 0.06 | 2.03 ± 0.06 | 1.93 ± 0.06 | 2.07 ± 0.03 |
| A 15.8 … 20.4 | 39 | 2.25 ± 0.05 | 2.47 ± 0.06 | 2.57 ± 0.07 | 2.34 ± 0.06 | 2.13 ± 0.06 | 2.34 ± 0.02 |
| r_p < 1 | 32 | 2.10 ± 0.06 | 2.14 ± 0.08 | 2.08 ± 0.09 | 1.96 ± 0.08 | 1.91 ± 0.07 | 2.06 ± 0.05 |
| r_p 1 … 2 | 50 | 2.15 ± 0.03 | 2.14 ± 0.06 | 2.21 ± 0.07 | 2.21 ± 0.06 | 2.04 ± 0.05 | 2.15 ± 0.03 |
| r_p ≥ 2 | 34 | 2.19 ± 0.05 | 2.18 ± 0.08 | 2.17 ± 0.08 | 1.95 ± 0.07 | 1.97 ± 0.07 | 2.12 ± 0.04 |
| L_eff 30 and A > 15.8 | 12 | 2.36 ± 0.05 | 2.53 ± 0.10 | 2.69 ± 0.10 | 2.55 ± 0.09 | 2.18 ± 0.07 | 2.44 ± 0.03 |
| L_eff 30 and A ≤ 12.5 | 7 | 1.84 ± 0.08 | 1.69 ± 0.12 | 1.62 ± 0.03 | 1.70 ± 0.06 | 1.73 ± 0.06 | 1.74 ± 0.05 |
| L_eff 12 and A > 15.8 | 9 | 2.15 ± 0.14 | 2.58 ± 0.12 | 2.51 ± 0.19 | 2.09 ± 0.17 | 2.08 ± 0.16 | 2.27 ± 0.05 |
| L_eff 12 and A ≤ 12.5 | 11 | 2.24 ± 0.04 | 2.02 ± 0.13 | 2.16 ± 0.09 | 2.17 ± 0.16 | 2.16 ± 0.14 | 2.17 ± 0.04 |
