# N0b (reference stray-field solver, PLAN 4.11): mesh convergence at 9 mT with a transferred state (PLAN 4.10, rules fixed before the runs)

Base microstructure (L_eff 12 l_ex, r_p 0, φ 0.65, α 0.1) at d/l_ex 96, 9 mT, seeds 941, 942, 943, 944. F: dx 0.75 l_ex. T: start from the relaxed F state (block average), then relaxed again. I: independent virgin relaxation from the same seed. Value: ln P (P in W/cm³) of the calm tail (≥ 2 calm cycles). Errors: SE over the seeds, 90 % CI with t(n − 1); SD of the paired difference with its 90 % CI (χ²).

## Runs

| arm | settled / runs | P (W/cm³), geometric mean | SD of ln P over the seeds | relaxation iterations (mean, range) | b_peak / b_target (range) | missing |
|---|---|---|---|---|---|---|
| F dx 0.75 | 4 / 4 | 4.5 | 0.663 | 11250 (8000 … 14500) | 0.9948 … 1.0080 | none |
| I dx 1.5 | 4 / 4 | 1.38 | 0.539 | 10000 (8500 … 13000) | 0.9933 … 1.0015 | none |
| I dx 3 | 4 / 4 | 1.13 | 0.247 | 4875 (4000 … 6000) | 1.0001 … 1.0004 | none |
| T dx 1.5 | 4 / 4 | 1.37 | 0.162 | 2750 (2500 … 3000) | 0.9989 … 0.9995 | none |
| T dx 3 | 4 / 4 | 1.79 | 1.126 | 2625 (2000 … 3500) | 0.9993 … 1.0052 | none |

Transfer: 8 runs, coarse particle cells without a fine sub-cell of the same particle: 0 (max), zero vectors: 0 (max).

Cycle noise (pooled SD of ln w_loop over the calm cycles): 0.009 (dof 40).

## Paired differences of ln P (same seed)

| difference | mean ± SE (90 % CI) | SD of the difference (90 % CI) | n |
|---|---|---|---|
| dx 1.5 − dx 0.75, transferred | -1.186 ± 0.382 (90 % CI -2.085 … -0.287, not inside ± 0.10) | 0.764 (90 % CI 0.473 … 2.231) | 4 |
| dx 3 − dx 1.5, transferred | +0.266 ± 0.500 (90 % CI -0.911 … +1.444, not inside ± 0.10) | 1.001 (90 % CI 0.620 … 2.922) | 4 |
| dx 3 − dx 0.75, transferred | -0.920 ± 0.650 (90 % CI -2.449 … +0.609, not inside ± 0.10) | 1.300 (90 % CI 0.805 … 3.795) | 4 |
| dx 1.5 − dx 0.75, independent (same seed) | -1.184 ± 0.551 (90 % CI -2.481 … +0.113, not inside ± 0.10) | 1.102 (90 % CI 0.683 … 3.219) | 4 |
| dx 3 − dx 1.5, independent (same seed) | -0.193 ± 0.168 (90 % CI -0.589 … +0.203, not inside ± 0.10) | 0.337 (90 % CI 0.209 … 0.983) | 4 |
| dx 1.5: transferred − independent | -0.002 ± 0.191 (90 % CI -0.451 … +0.447, not inside ± 0.10) | 0.382 (90 % CI 0.236 … 1.114) | 4 |
| dx 3: transferred − independent | +0.457 ± 0.468 (90 % CI -0.644 … +1.558, not inside ± 0.10) | 0.936 (90 % CI 0.580 … 2.733) | 4 |

## Unpaired differences of the arm means (Welch)

| difference | mean ± SE (90 % CI) | n |
|---|---|---|
| I dx 1.5 − F dx 0.75 | -1.184 ± 0.427 (-2.021 … -0.348) | 4 / 4 |
| I dx 3 − I dx 1.5 | -0.193 ± 0.296 (-0.816 … +0.430) | 4 / 4 |
| I dx 3 − F dx 0.75 | -1.377 ± 0.354 (-2.142 … -0.612) | 4 / 4 |

## Reading (pre-registered)

**stop: SD_T 0.764 ≥ 0.10, the transferred state does not remove the scatter → stop and re-scope**

Information only (9 mT, one amplitude, no β): Δ(1.5 − 0.75) of the T arm inside ± 0.10 means dx 1.5 agrees with dx 0.75 at 9 mT within 10 %; it is not the 4.10 criterion for β.

## Stage 1 reading (PLAN 4.11, pre-registered)

| difference | mean ± SE (90 % CI) | SD of the difference (90 % CI) | n |
|---|---|---|---|
| ln P(I3) − ln P(I15), N0b | -0.193 ± 0.168 (90 % CI -0.589 … +0.203, not inside ± 0.10) | 0.337 (90 % CI 0.209 … 0.983) | 4 |
| Δ_N − Δ_FD (B0b, same seed) | +0.599 ± 0.192 (90 % CI +0.148 … +1.050, not inside ± 0.10) | 0.383 (90 % CI 0.237 … 1.119) | 4 |

**much smaller: |mean Δ_N| 0.193 ≤ 0.30**
