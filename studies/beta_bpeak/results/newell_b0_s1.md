# N0b (reference stray-field solver, PLAN 4.11): mesh convergence at 9 mT with a transferred state (PLAN 4.10, rules fixed before the runs)

Base microstructure (L_eff 12 l_ex, r_p 0, φ 0.65, α 0.1) at d/l_ex 96, 9 mT, seeds 941, 942, 943, 944. F: dx 0.75 l_ex. T: start from the relaxed F state (block average), then relaxed again. I: independent virgin relaxation from the same seed. Value: ln P (P in W/cm³) of the calm tail (≥ 2 calm cycles). Errors: SE over the seeds, 90 % CI with t(n − 1); SD of the paired difference with its 90 % CI (χ²).

## Runs

| arm | settled / runs | P (W/cm³), geometric mean | SD of ln P over the seeds | relaxation iterations (mean, range) | b_peak / b_target (range) | missing |
|---|---|---|---|---|---|---|
| F dx 0.75 | 0 / 4 | – | – | – | – | s941: not run; s942: not run; s943: not run; s944: not run |
| I dx 1.5 | 4 / 4 | 1.34 | 0.534 | 10000 (8500 … 13000) | 0.9410 … 1.0015 | none |
| I dx 3 | 4 / 4 | 1.11 | 0.289 | 4875 (4000 … 6000) | 0.9495 … 1.0004 | none |
| T dx 1.5 | 0 / 4 | – | – | – | – | s941: not run; s942: not run; s943: not run; s944: not run |
| T dx 3 | 0 / 4 | – | – | – | – | s941: not run; s942: not run; s943: not run; s944: not run |

Cycle noise (pooled SD of ln w_loop over the calm cycles): 0.015 (dof 16).

## Paired differences of ln P (same seed)

| difference | mean ± SE (90 % CI) | SD of the difference (90 % CI) | n |
|---|---|---|---|
| dx 1.5 − dx 0.75, transferred | n 0: no statement |
| dx 3 − dx 1.5, transferred | n 0: no statement |
| dx 3 − dx 0.75, transferred | n 0: no statement |
| dx 1.5 − dx 0.75, independent (same seed) | n 0: no statement |
| dx 3 − dx 1.5, independent (same seed) | -0.194 ± 0.133 (90 % CI -0.506 … +0.119, not inside ± 0.10) | 0.266 (90 % CI 0.165 … 0.776) | 4 |
| dx 1.5: transferred − independent | n 0: no statement |
| dx 3: transferred − independent | n 0: no statement |

## Unpaired differences of the arm means (Welch)

| difference | mean ± SE (90 % CI) | n |
|---|---|---|
| I dx 1.5 − F dx 0.75 | too few values | – |
| I dx 3 − I dx 1.5 | -0.194 ± 0.304 (-0.817 … +0.430) | 4 / 4 |
| I dx 3 − F dx 0.75 | too few values | – |

## Reading (pre-registered)

**not determined: 0 seeds with a value in F and T15 (need 3)**

Information only (9 mT, one amplitude, no β): Δ(1.5 − 0.75) of the T arm inside ± 0.10 means dx 1.5 agrees with dx 0.75 at 9 mT within 10 %; it is not the 4.10 criterion for β.

## Stage 1 reading (PLAN 4.11, pre-registered)

| difference | mean ± SE (90 % CI) | SD of the difference (90 % CI) | n |
|---|---|---|---|
| ln P(I3) − ln P(I15), N0b | -0.194 ± 0.133 (90 % CI -0.506 … +0.119, not inside ± 0.10) | 0.266 (90 % CI 0.165 … 0.776) | 4 |
| Δ_N − Δ_FD (B0b, same seed) | +0.598 ± 0.202 (90 % CI +0.122 … +1.074, not inside ± 0.10) | 0.404 (90 % CI 0.251 … 1.181) | 4 |

**much smaller: |mean Δ_N| 0.194 ≤ 0.30**
