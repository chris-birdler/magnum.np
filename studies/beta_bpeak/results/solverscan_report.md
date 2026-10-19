# Solver scan: DemagFieldPBC against Newell at all amplitudes from the same states (eval_solverscan.py)

Δ_k = ln w_FD − ln w_Newell per stage (same state at the start of the stage).

| corner run | 9 mT | 20 mT | 35 mT | 50 mT | 70 mT | 100 mT | 150 mT | β Newell | β FD | Δβ_solver | β old − β new (corner test) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| K1_0_nx | +0.092 | -0.196 | +0.444 | +0.740 | +0.997 | +0.654 | +0.548 | 2.112 | 2.407 | +0.296 | +0.216 |
| K2_0_nx | +0.514 | +0.646 | -0.200 | -0.315 | -0.440 | -0.201 | +0.199 | 2.285 | 2.027 | -0.258 | -0.192 |
| K3_0_nx | +0.390 | +0.329 | +0.052 | -0.031 | +0.143 | +0.503 | +0.513 | 2.313 | 2.349 | +0.036 | +0.065 |
| K4_0_nx | +0.150 | -0.130 | +0.099 | +0.018 | -0.187 | -0.257 | -0.221 | 2.020 | 1.891 | -0.129 | -0.245 |

Mean Δ_k over the corners (± SE, n 4): 9 mT +0.287 ± 0.100; 20 mT +0.162 ± 0.199; 35 mT +0.099 ± 0.132; 50 mT +0.103 ± 0.225; 70 mT +0.128 ± 0.313; 100 mT +0.175 ± 0.235; 150 mT +0.260 ± 0.178.

Corner differences β(Ki) − β(K1), Newell / FD: K2 +0.173 / -0.380; K3 +0.201 / -0.058; K4 -0.092 / -0.516.

Mean Δβ_solver: -0.014 ± 0.119 (SE, n 4; 90 % CI -0.295 … +0.267). **Reading: not determined.**
