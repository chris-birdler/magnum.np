# N1g (final stray-field setting, PLAN 4.12 V3 / V3g): mesh convergence at 9 mT with a transferred state (PLAN 4.10, rules fixed before the runs)

Base microstructure (L_eff 12 l_ex, r_p 0, φ 0.65, α 0.1) at d/l_ex 96, 9 mT, seeds 941, 942, 943, 944. F: dx 0.75 l_ex. T: start from the relaxed F state (block average), then relaxed again. I: independent virgin relaxation from the same seed. Value: ln P (P in W/cm³) of the calm tail (≥ 2 calm cycles). Errors: SE over the seeds, 90 % CI with t(n − 1); SD of the paired difference with its 90 % CI (χ²).

## Runs

| arm | settled / runs | P (W/cm³), geometric mean | SD of ln P over the seeds | relaxation iterations (mean, range) | b_peak / b_target (range) | missing |
|---|---|---|---|---|---|---|
| F dx 0.75 | 4 / 4 | 1.42 | 0.173 | 4125 (4000 … 4500) | 0.9997 … 1.0007 | none |
| I dx 1.5 | 4 / 4 | 1.3 | 0.227 | 4500 (3500 … 6000) | 1.0000 … 1.0013 | none |
| I dx 3 | 4 / 4 | 1.11 | 0.183 | 5375 (4000 … 6500) | 0.9974 … 1.0003 | none |
| T dx 1.5 | 4 / 4 | 1.41 | 0.157 | 1500 (1500 … 1500) | 1.0000 … 1.0008 | none |
| T dx 3 | 4 / 4 | 1.19 | 0.220 | 2250 (1500 … 4000) | 0.9933 … 1.0000 | none |

Transfer: 8 runs, coarse particle cells without a fine sub-cell of the same particle: 0 (max), zero vectors: 0 (max).

Cycle noise (pooled SD of ln w_loop over the calm cycles): 0.000 (dof 40).

## Paired differences of ln P (same seed)

| difference | mean ± SE (90 % CI) | SD of the difference (90 % CI) | n |
|---|---|---|---|
| dx 1.5 − dx 0.75, transferred | -0.008 ± 0.011 (90 % CI -0.033 … +0.018, inside ± 0.10) | 0.022 (90 % CI 0.013 … 0.063) | 4 |
| dx 3 − dx 1.5, transferred | -0.168 ± 0.059 (90 % CI -0.307 … -0.029, not inside ± 0.10) | 0.118 (90 % CI 0.073 … 0.345) | 4 |
| dx 3 − dx 0.75, transferred | -0.175 ± 0.065 (90 % CI -0.329 … -0.022, not inside ± 0.10) | 0.131 (90 % CI 0.081 … 0.381) | 4 |
| dx 1.5 − dx 0.75, independent (same seed) | -0.092 ± 0.086 (90 % CI -0.295 … +0.111, not inside ± 0.10) | 0.173 (90 % CI 0.107 … 0.504) | 4 |
| dx 3 − dx 1.5, independent (same seed) | -0.155 ± 0.050 (90 % CI -0.273 … -0.038, not inside ± 0.10) | 0.100 (90 % CI 0.062 … 0.292) | 4 |
| dx 1.5: transferred − independent | +0.085 ± 0.078 (90 % CI -0.099 … +0.268, not inside ± 0.10) | 0.156 (90 % CI 0.097 … 0.456) | 4 |
| dx 3: transferred − independent | +0.072 ± 0.024 (90 % CI +0.016 … +0.128, not inside ± 0.10) | 0.047 (90 % CI 0.029 … 0.138) | 4 |

## Unpaired differences of the arm means (Welch)

| difference | mean ± SE (90 % CI) | n |
|---|---|---|
| I dx 1.5 − F dx 0.75 | -0.092 ± 0.142 (-0.372 … +0.188) | 4 / 4 |
| I dx 3 − I dx 1.5 | -0.155 ± 0.146 (-0.441 … +0.130) | 4 / 4 |
| I dx 3 − F dx 0.75 | -0.247 ± 0.126 (-0.492 … -0.003) | 4 / 4 |

## Reading (pre-registered)

**go: SD_T 0.022 < 0.10 and SD_I 0.173 ≥ 0.15, the transferred state removes the state scatter → variant B (separate go)**

Information only (9 mT, one amplitude, no β): Δ(1.5 − 0.75) of the T arm inside ± 0.10 means dx 1.5 agrees with dx 0.75 at 9 mT within 10 %; it is not the 4.10 criterion for β.

## Convergence reading (PLAN 4.11 / 4.12, pre-registered)

- T arm (decisive): Δ(1.5 − 0.75) -0.008 ± 0.011 (90 % CI -0.033 … +0.018, inside ± 0.10) | 0.022 (90 % CI 0.013 … 0.063) | 4; Δ(3 − 1.5) -0.168 ± 0.059 (90 % CI -0.307 … -0.029, not inside ± 0.10) | 0.118 (90 % CI 0.073 … 0.345) | 4. **Reading: converged at dx 1.5.**
- I arm (information): Δ(1.5 − 0.75) -0.092 ± 0.086 (90 % CI -0.295 … +0.111, not inside ± 0.10) | 0.173 (90 % CI 0.107 … 0.504) | 4; Δ(3 − 1.5) -0.155 ± 0.050 (90 % CI -0.273 … -0.038, not inside ± 0.10) | 0.100 (90 % CI 0.062 … 0.292) | 4. **Reading: not determined.**

## Information: μ_r = b_peak / h and the energy balance

| arm | μ_r (geometric mean) | SD of ln μ_r |
|---|---|---|
| F | 8.285 | 0.0093 |
| I15 | 8.265 | 0.0078 |
| I3 | 8.170 | 0.0060 |
| T15 | 8.269 | 0.0097 |
| T3 | 8.107 | 0.0124 |
- ln μ_r(T3) − ln μ_r(T15): -0.0198 ± 0.0031 (SE, n 4)
- ln μ_r(T15) − ln μ_r(F): -0.0020 ± 0.0002 (SE, n 4)
- ln μ_r(I3) − ln μ_r(I15): -0.0115 ± 0.0035 (SE, n 4)
- ln μ_r(I15) − ln μ_r(F): -0.0025 ± 0.0013 (SE, n 4)
- ε of the F tails (s941 … s944): +0.025 / +0.031 / +0.024 / +0.045.
