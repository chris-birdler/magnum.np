# Audit: conditioning test, determinism, slow drifts

## (1) Conditioning test recomputed (centre point, C0 … C5)

| run | design P (W/cm³, calm tail of cycles 2 … 6) | late P (cycles 28, 29) | d = ln late/design | events (cycle) | trend of ln P over the calm cycles ≥ 5 (per 100 cycles, ± SE) | closure first / last calm cycle |
|---|---|---|---|---|---|---|
| cond_C0_b9 | 8.06 | 8.48 | +0.051 | - | +0.057 ± 0.013 | 1.1e-04 / 1.9e-07 |
| cond_C1_b9 | 6.75 | 6.55 | -0.030 | 2,3 | +0.026 ± 0.016 | 2.4e-04 / 3.5e-07 |
| cond_C2_b9 | 3.65 | 3.66 | +0.004 | - | -0.001 ± 0.001 | 2.8e-06 / 1.1e-08 |
| cond_C3_b9 | 4.47 | 4.5 | +0.008 | - | +0.001 ± 0.000 | 2.0e-06 / 3.1e-07 |
| cond_C4_b9 | 3.83 | 3.86 | +0.008 | - | +0.001 ± 0.001 | 2.0e-06 / 1.7e-07 |
| cond_C5_b9 | 6.54 | 6.32 | -0.035 | 2 | +0.023 ± 0.018 | 4.2e-04 / 1.2e-06 |
| cond_C0_b20 | 38.1 | 37.9 | -0.006 | - | -0.033 ± 0.005 | 2.7e-04 / 4.8e-07 |
| cond_C1_b20 | 40.6 | 39.8 | -0.019 | - | -0.015 ± 0.058 | 3.1e-04 / 3.3e-04 |
| cond_C2_b20 | 32.6 | 32.9 | +0.008 | - | +0.008 ± 0.098 | 2.9e-03 / 2.9e-03 |
| cond_C3_b20 | 24.8 | 24.8 | -0.002 | - | -0.000 ± 0.000 | 1.8e-05 / 3.9e-07 |
| cond_C4_b20 | 19.6 | 19.4 | -0.008 | - | -0.000 ± 0.000 | 1.7e-07 / 2.4e-07 |
| cond_C5_b20 | 23.1 | 23.6 | +0.018 | 3,4 | +0.034 ± 0.017 | 1.3e-05 / 2.2e-07 |
| cond_C0_b35 | 117 | 115 | -0.012 | - | +0.084 ± 0.105 | 9.9e-04 / 1.9e-03 |
| cond_C1_b35 | 85.1 | 85 | -0.001 | - | +0.027 ± 0.033 | 1.2e-03 / 4.8e-04 |
| cond_C2_b35 | 68.5 | 66.2 | -0.034 | 4 | -0.061 ± 0.034 | 1.6e-04 / 1.7e-07 |
| cond_C3_b35 | 113 | 112 | -0.007 | - | -0.022 ± 0.046 | 6.6e-04 / 9.0e-05 |
| cond_C4_b35 | 94.3 | 97.5 | +0.034 | - | +0.055 ± 0.121 | 2.2e-03 / 2.1e-03 |
| cond_C5_b35 | 70.6 | 70.4 | -0.003 | - | +0.000 ± 0.000 | 6.8e-05 / 8.7e-08 |

| B (mT) | mean d ± SE (n 6) | 90 % t-CI | mean trend per 100 cycles ± SE | 95 % t-CI of the trend |
|---|---|---|---|---|
| 9 | +0.0010 ± 0.0126 | -0.024 … +0.027 | +0.018 ± 0.009 | -0.006 … +0.042 |
| 20 | -0.0015 ± 0.0054 | -0.012 … +0.009 | -0.001 ± 0.009 | -0.025 … +0.022 |
| 35 | -0.0038 ± 0.0089 | -0.022 … +0.014 | +0.014 ± 0.021 | -0.041 … +0.069 |

## (2) Drift inside the calm tails of the 116 design runs (main rule, reruns; stages with ≥ 4 kept cycles)

| B (mT) | n stages | mean slope of ln w per cycle ± SE | stages with negative slope | median closure of the kept cycles |
|---|---|---|---|---|
| 9 | 100 | -0.0032 ± 0.0022 | 55 | 7.3e-05 |
| 20 | 107 | -0.0009 ± 0.0023 | 58 | 1.5e-04 |
| 35 | 97 | -0.0012 ± 0.0013 | 56 | 2.7e-04 |
| 50 | 107 | +0.0018 ± 0.0028 | 57 | 5.2e-04 |
| 70 | 112 | +0.0014 ± 0.0016 | 47 | 4.2e-04 |
| 100 | 112 | -0.0008 ± 0.0014 | 56 | 6.7e-04 |
| 150 | 116 | -0.0000 ± 0.0011 | 58 | 6.5e-04 |

## (3) Determinism

- cond runs: the earlier stages and the first 7 cycles of the tested stage are bit-identical to C0 … C5 (all 12 checks), unless listed above.
