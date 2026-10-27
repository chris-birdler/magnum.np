# Dynamic arms at d/l_ex 96, φ 0.55, 9 mT (PROTOKOLL 3.5, D1 and D2)

| arm | meaning | P (W/cm³) per seed 941 / 942 / 943 / 944 | geometric mean | SD of ln P | status / Bloch points of the prepared state |
|---|---|---|---|---|---|
| D3 | dx 3, voxel model (production) | 1.02 / 0.799 / 0.697 / 0.775 | 0.815 | 0.16 | settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0 |
| D3f | dx 3, volume fraction alone, q 2 | 0.882 / 12.1 / 0.64 / 0.717 | 1.49 | 1.40 | settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0 |
| D3Pc | dx 3, pairwise correction with centroid term, q 2, p 3 (main arm, D4/D8) | 1.27 / 0.594 / 3.03 / 0.585 | 1.08 | 0.78 | settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0 |
| D3Pc4 | dx 3, pairwise correction with centroid term, q 4, p 3 | – / – / – / – | – | – | not run ?; not run ?; not run ?; not run ? |
| D3Lc | dx 3, cell-local correction with centroid term, q 2, p 3 (not usable by R3; information) | 0.857 / 0.581 / – / – | 0.706 | 0.27 | settled +0 -0; settled +0 -0; started, no summary (refused or crashed) ?; not run ? |
| D3P | dx 3, pairwise correction without centroid term, q 2, p 3 (information) | 1.43 / – / – / – | 1.43 | – | settled +0 -0; not run ?; not run ?; not run ? |
| D15 | dx 1.5, voxel model of dx 1.5 | 1.88 / 1.39 / 0.623 / 1.16 | 1.17 | 0.47 | settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0 |
| D15g | dx 1.5, voxel geometry of dx 3 locked | 1.08 / 0.871 / 1.1 / 0.875 | 0.976 | 0.13 | settled +0 -0; settled +0 -0; settled +0 -0; settled +0 -0 |

Design A. Geometry effect of the fine staircase at dx 1.5, Δ_g = ln P(D15) − ln P(D15g): +0.183 ± 0.256 (SE, n 4)
Resolution of m at the dx 3 geometry, ln P(D15g) − ln P(D3): +0.180 ± 0.092 (SE, n 4)

| corrected arm | Δ_c = ln P(arm) − ln P(D3) | d = Δ_c − Δ_g (per seed) | 95 % interval of mean d (t) | reading |
|---|---|---|---|---|
| D3Pc | +0.278 ± 0.415 (SE, n 4) | +0.095 ± 0.654 (SE, n 4) | -1.99 … +2.17 | not determined |
| D3Pc4 | n 0 | n 0 | – | not determined (n 0) (information) |
| D3Lc | -0.247 ± 0.072 (SE, n 2) | -0.757 ± 0.030 (SE, n 2) | -1.14 … -0.37 | not determined (n 2) (information) |
| D3P | n 1 | n 1 | – | not determined (n 1) (information) |
| D3f | +0.602 ± 0.705 (SE, n 4) | +0.418 ± 0.658 (SE, n 4) | -1.68 … +2.51 | not determined (information) |
