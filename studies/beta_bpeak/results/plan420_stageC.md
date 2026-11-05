# PLAN 4.20 stage C: variant V2 in the 4.16 geometry

| job | done | relaxed conv. | Bloch relaxed | Bloch samples (cycles) | complete | w | P (W/cm^3) | run time (h) |
|---|---|---|---|---|---|---|---|---|
| DV2_s921_b9 | True | True | False | 0 [] | True | 2.3382e-07 | 6.28 | 2.03 |
| DV2_s921_b50 | True | True | False | 0 [] | True | 4.3508e-06 | 116.9 | 2.44 |
| DV2_s921_b100 | True | True | False | 0 [] | True | 2.1597e-05 | 580 | 2.26 |
| DV2_s922_b9 | True | True | True | 150 [0] | True | 2.7881e-07 | 7.488 | 1.87 |
| DV2_s922_b50 | True | True | True | 67 [0] | True | 5.7911e-06 | 155.5 | 1.90 |
| DV2_s922_b100 | True | True | True | 55 [0] | True | 1.9861e-05 | 533.4 | 2.06 |
| DV2_s923_b9 | True | True | False | 0 [] | True | 2.3554e-07 | 6.326 | 1.31 |
| DV2_s923_b50 | True | True | False | 0 [] | True | 4.575e-06 | 122.9 | 1.35 |
| DV2_s923_b100 | True | True | False | 0 [] | True | 1.759e-05 | 472.4 | 1.48 |
| DV2_s924_b9 | True | True | False | 0 [] | True | 3.6088e-07 | 9.692 | 1.53 |
| DV2_s924_b50 | True | True | False | 0 [] | True | 5.3049e-06 | 142.5 | 1.56 |
| DV2_s924_b100 | True | True | False | 0 [] | True | 2.4124e-05 | 647.9 | 1.74 |

- C-a (q 8): counted 9 jobs, ran their cycles 9, fail none: holds
- C-a (q 16): D16 has Bloch samples in seeds [921, 923, 924]; counted 3, ran 3, fail none: holds

**Verdict: no Bloch points in these runs**

## C-b: Delta ln w and Delta beta (paired by seed, mean +- SE, 95 % t-interval)

- w Dv/D_9mT: factor 0.934 [0.890, 0.979], n 4
- w Dv/D_50mT: factor 1.029 [0.969, 1.094], n 4
- w Dv/D_100mT: factor 0.968 [0.901, 1.039], n 4
- w Dv/C_9mT: factor 1.276 [0.772, 2.109], n 4
- w Dv/C_50mT: factor 1.135 [0.753, 1.712], n 4
- w Dv/C_100mT: factor 1.394 [1.129, 1.722], n 4
- beta Dv-D_9_50: +0.057 +- 0.010 [+0.024, +0.090], n 4
- beta Dv-D_50_100: -0.089 +- 0.038 [-0.211, +0.033], n 4
- beta Dv-D_9_100: +0.015 +- 0.010 [-0.018, +0.048], n 4
- beta Dv-C_9_50: -0.068 +- 0.032 [-0.170, +0.033], n 4
- beta Dv-C_50_100: +0.296 +- 0.164 [-0.224, +0.817], n 4
- beta Dv-C_9_100: +0.037 +- 0.055 [-0.140, +0.213], n 4

- R0: w 2.2130246e-05 against 2.2130246e-05 of 4.16, relative nan: passes
- HFP_s922: 25500 iterations, E/K_d 0.00016218: equal to 4.16
