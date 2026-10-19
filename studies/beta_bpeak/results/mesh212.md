# Mesh error at d/l_ex 212 and 300 (base microstructure), dx 3 against dx 1.5

Base microstructure: L_eff 12 l_ex, r_p 0, φ 0.65, α 0.1, A 10 pJ/m (l_ex 3.34 nm; the µm values hold for this A). Steady loops (≥ 2 calm cycles), reruns replace their original. dx in units of l_ex. Errors: SE over the seeds, 90 % CI with t(n − 1); fewer than 3 seeds: no statement.

## Loss density P (W/cm³, geometric mean over the seeds)

| group | mesh | 9 mT | 50 mT | 100 mT |
|---|---|---|---|---|
| d/l_ex 212 (0.71 µm) | dx 1.5 | 6.6 (n 4) | 227 (n 4) | 727 (n 4) |
| d/l_ex 212 (0.71 µm) | dx 3 | 2.83 (n 4) | 206 (n 4) | 1.12e+03 (n 4) |
| d/l_ex 300 (1.0 µm) | dx 1.5 | 5.23 (n 4) | 134 (n 4) | 496 (n 4) |
| d/l_ex 300 (1.0 µm) | dx 3 | 7.73 (n 4) | 188 (n 4) | 619 (n 4) |
| K3: d/l_ex 212, φ 0.69, L_eff 18, r_p 1.5 | dx 1.5 | 4.28 (n 6) | 213 (n 6) | 948 (n 6) |
| K3: d/l_ex 212, φ 0.69, L_eff 18, r_p 1.5 | dx 3 | 2.68 (n 6) | 148 (n 6) | 826 (n 6) |

## Mesh error (paired per seed)

| group | P(dx 3)/P(dx 1.5) 9 mT | 50 mT | 100 mT | Δβ (9 / 50 / 100 mT) | Δβ (9 → 50 mT) |
|---|---|---|---|---|---|
| d/l_ex 212 (0.71 µm) | × 0.43 (ln -0.85 ± 0.19, CI -1.29 … -0.40, n 4) | × 0.91 (ln -0.10 ± 0.12, CI -0.39 … +0.20, n 4) | × 1.54 (ln +0.43 ± 0.16, CI +0.06 … +0.80, n 4) | +0.51 ± 0.13 (CI +0.20 … +0.81, n 4) | +0.43 ± 0.18 (CI +0.01 … +0.86, n 4) |
| d/l_ex 300 (1.0 µm) | × 1.48 (ln +0.39 ± 0.20, CI -0.08 … +0.86, n 4) | × 1.40 (ln +0.34 ± 0.16, CI -0.03 … +0.71, n 4) | × 1.25 (ln +0.22 ± 0.06, CI +0.07 … +0.37, n 4) | -0.06 ± 0.05 (CI -0.19 … +0.07, n 4) | -0.02 ± 0.03 (CI -0.09 … +0.05, n 4) |
| K3: d/l_ex 212, φ 0.69, L_eff 18, r_p 1.5 | × 0.63 (ln -0.47 ± 0.09, CI -0.64 … -0.29, n 6) | × 0.70 (ln -0.36 ± 0.13, CI -0.63 … -0.09, n 6) | × 0.87 (ln -0.14 ± 0.09, CI -0.32 … +0.04, n 6) | +0.12 ± 0.03 (CI +0.06 … +0.18, n 6) | +0.05 ± 0.05 (CI -0.04 … +0.14, n 6) |

Stages without a steady value (≥ 2 calm cycles): M212: none; M300: none; MK3: none

## Mesh share of the size lever (d/l_ex 212 → 300, base microstructure)

Size lever on one mesh: mean ln P (or β) at 300 minus at 212, over the seeds with a value on both meshes (different seeds at 212 and 300, Welch). Mesh share = lever at dx 3 − lever at dx 1.5 = Δ(300) − Δ(212). Sign change of the mesh error: opposite signs and both 90 % CIs exclude 0.

| output | lever at dx 1.5 | lever at dx 3 | mesh share Δ(300) − Δ(212) | sign change |
|---|---|---|---|---|
| lnP 9 | -0.23 ± 0.22 (CI -0.69 … +0.23) | +1.00 ± 0.30 (CI +0.34 … +1.67) | +1.24 ± 0.27 (CI +0.70 … +1.77) | no |
| lnP 50 | -0.53 ± 0.09 (CI -0.71 … -0.35) | -0.09 ± 0.17 (CI -0.42 … +0.23) | +0.44 ± 0.20 (CI +0.05 … +0.83) | no |
| lnP 100 | -0.38 ± 0.12 (CI -0.65 … -0.11) | -0.59 ± 0.11 (CI -0.82 … -0.37) | -0.21 ± 0.17 (CI -0.57 … +0.15) | no |
| beta | -0.09 ± 0.12 (CI -0.35 … +0.17) | -0.65 ± 0.13 (CI -0.91 … -0.40) | -0.56 ± 0.14 (CI -0.86 … -0.26) | no |
| beta 9-50 | -0.18 ± 0.16 (CI -0.49 … +0.13) | -0.63 ± 0.16 (CI -0.95 … -0.32) | -0.46 ± 0.18 (CI -0.87 … -0.04) | no |

## Old pairs (relaxation not converged; hints only, PROTOKOLL 7.64)

| group | P(dx a)/P(dx 1.5) 9 mT | 50 mT | 100 mT | Δβ (9 / 50 / 100 mT) | Δβ (9 → 50 mT) |
|---|---|---|---|---|---|
| d/l_ex 96, dx a = 3 | n 2: no statement | × 0.33 (ln -1.12 ± 0.67, CI -3.06 … +0.82, n 3) | – | – | n 2: no statement |
| d/l_ex 96, dx a = 2 | n 2: no statement | × 1.30 (ln +0.26 ± 0.33, CI -0.70 … +1.22, n 3) | – | – | n 2: no statement |
| d/l_ex 96, dx a = 1 | – | × 1.64 (ln +0.50 ± 0.36, CI -0.55 … +1.55, n 3) | – | – | – |
| d/l_ex 212, dx a = 3 | × 0.70 (ln -0.36 ± 0.34, CI -1.15 … +0.43, n 4) | × 0.94 (ln -0.06 ± 0.09, CI -0.28 … +0.16, n 4) | – | – | +0.15 ± 0.23 (CI -0.38 … +0.68, n 4) |
