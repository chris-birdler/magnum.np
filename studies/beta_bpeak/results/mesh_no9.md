# Mesh check without the 9 mT point (post hoc, request of Chris 2026-10-12)

β from 50 and 100 mT only (and for K3 also from 20 … 100 mT). dx in l_ex. ± SE, 90 % CI with t. Fewer than 3 seeds: no statement. The full table with 9 mT: results/mesh212.md.

## β (50 → 100 mT) per mesh and mesh error

| point | β fine dx 1.5 | β coarse dx 3 | Δβ = coarse − fine (paired) |
|---|---|---|---|
| d/l_ex 212 (0.71 µm) | +1.66 ± 0.12 (CI +1.38 … +1.95, n 4) | +2.40 ± 0.15 (CI +2.06 … +2.74, n 4) | +0.73 ± 0.20 (CI +0.27 … +1.19, n 4) |
| d/l_ex 300 (1.0 µm) | +1.81 ± 0.06 (CI +1.68 … +1.94, n 4) | +1.72 ± 0.07 (CI +1.56 … +1.88, n 4) | -0.10 ± 0.11 (CI -0.36 … +0.17, n 4) |
| K3: d/l_ex 212, φ 0.69, L_eff 18, r_p 1.5 | +2.13 ± 0.17 (CI +1.79 … +2.48, n 6) | +2.49 ± 0.15 (CI +2.20 … +2.79, n 6) | +0.36 ± 0.11 (CI +0.14 … +0.58, n 6) |
| K3, 20 … 100 mT | +2.38 ± 0.08 (CI +2.23 … +2.54, n 6) | +2.43 ± 0.05 (CI +2.33 … +2.53, n 6) | +0.05 ± 0.07 (CI -0.09 … +0.19, n 6) |
| K3, 50 → 100 mT, eval_k3fill runs | +2.13 ± 0.17 (CI +1.79 … +2.48, n 6) | +2.49 ± 0.15 (CI +2.20 … +2.79, n 6) | +0.36 ± 0.11 (CI +0.14 … +0.58, n 6) |

## Size lever 212 → 300 (base microstructure; value at 1 µm minus at 0.71 µm, Welch)

| output | lever fine dx 1.5 | lever coarse dx 3 | mesh share (coarse − fine) |
|---|---|---|---|
| β (50 → 100 mT) | +0.15 ± 0.13 (CI -0.13 … +0.43) | -0.68 ± 0.16 (CI -1.02 … -0.34) | -0.83 ± 0.23 (CI -1.29 … -0.37) |
