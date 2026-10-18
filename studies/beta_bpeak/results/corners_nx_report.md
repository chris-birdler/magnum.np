# Corner test: the 12 corner runs with the production setup of PLAN 4.13 (eval_corners_nx.py)

β = slope of ln w over ln B, 9 … 150 mT (main-fit definition). Old: DemagFieldPBC, atol 1e-5 (the _mc2 rerun where it exists); new: PLAN 4.13. Δβ = new − old, paired per run.

| run (old partner) | β old | β new | Δβ |
|---|---|---|---|
| K1_0 (K1_0_mc2) | 2.328 | 2.111 | -0.216 |
| K2_0 (K2_0) | 2.093 | 2.284 | +0.192 |
| K3_0 (K3_0) | 2.378 | 2.313 | -0.065 |
| K4_0 (K4_0_mc2) | 1.775 | 2.020 | +0.245 |
| K1_1 (K1_1) | 2.553 | 2.127 | -0.426 |
| K2_1 (K2_1) | 2.330 | 2.387 | +0.056 |
| K3_1 (K3_1) | 2.464 | 1.759 | -0.705 |
| K4_1 (K4_1_mc2) | 1.993 | 2.224 | +0.232 |
| K1_2 (K1_2) | 2.482 | 2.404 | -0.077 |
| K2_2 (K2_2) | 2.209 | 2.103 | -0.107 |
| K3_2 (K3_2) | 2.294 | 2.078 | -0.216 |
| K4_2 (K4_2) | 1.813 | 1.962 | +0.148 |

## Per corner (mean ± SE, n 3)

| corner | β old | β new | Δβ |
|---|---|---|---|
| K1 | 2.454 ± 0.067 | 2.214 ± 0.095 | -0.240 ± 0.101 |
| K2 | 2.211 ± 0.069 | 2.258 ± 0.083 | +0.047 ± 0.086 |
| K3 | 2.379 ± 0.049 | 2.050 ± 0.161 | -0.328 ± 0.193 |
| K4 | 1.860 ± 0.067 | 2.069 ± 0.080 | +0.208 ± 0.030 |

Corner differences β(Ki) − β(K1), old / new: K2 -0.243 / +0.044; K3 -0.076 / -0.164; K4 -0.594 / -0.146.

Mean Δ ln w per amplitude (new − old, SE, n): 9 mT -0.170 ± 0.104 (12); 20 mT -0.277 ± 0.152 (12); 35 mT -0.158 ± 0.210 (12); 50 mT -0.163 ± 0.195 (12); 70 mT -0.407 ± 0.176 (12); 100 mT -0.376 ± 0.147 (12); 150 mT -0.370 ± 0.114 (12).

Mean Δβ over the runs: -0.078 ± 0.082 (SE, n 12; 90 % CI -0.226 … +0.070). **Reading: not determined.**

## Information (PLAN 4.13)

- Flags (S9) of the new stages: flag ε 15, flag w 9 of 84 stages. Mean Δβ without the flagged stages: -0.041 ± 0.127 (SE, n 7).

Energy-terms check (S5) on the calm-tail cycles of the new runs: d = ε + ΔE_int / w_loop.

| amplitude | cycles | median \|ε\| | median \|d\| | max \|d\| |
|---|---|---|---|---|
| 9 mT | 59 | 0.0007 | 0.0031 | 0.0166 |
| 20 mT | 60 | 0.0016 | 0.0009 | 0.0036 |
| 35 mT | 58 | 0.0050 | 0.0009 | 0.0064 |
| 50 mT | 58 | 0.0028 | 0.0006 | 0.0104 |
| 70 mT | 57 | 0.0026 | 0.0013 | 0.0080 |
| 100 mT | 58 | 0.0067 | 0.0030 | 0.0130 |
| 150 mT | 58 | 0.0098 | 0.0018 | 0.0178 |

Relaxation iterations of the new runs: K1_0 8000; K2_0 32500; K3_0 22500; K4_0 19500; K1_1 12000; K2_1 33000; K3_1 14500; K4_1 13000; K1_2 14000; K2_2 43000; K3_2 25500; K4_2 48000.
