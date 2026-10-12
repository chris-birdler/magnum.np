# Mesh error at d/l_ex 212 and 300 (base microstructure), dx 3 against dx 1.5

Base microstructure: L_eff 12 l_ex, r_p 0, φ 0.65, α 0.1, A 10 pJ/m (l_ex 3.34 nm; the µm values hold for this A). Steady loops (≥ 2 calm cycles), reruns replace their original. dx in units of l_ex. Errors: SE over the seeds, 90 % CI with t(n − 1); fewer than 3 seeds: no statement.

## Loss density P (W/cm³, geometric mean over the seeds)

| group | mesh | 9 mT | 50 mT | 100 mT |
|---|---|---|---|---|
| d/l_ex 212 (0.71 µm) | dx 1.5 | – | – | – |
| d/l_ex 212 (0.71 µm) | dx 3 | – | – | – |
| d/l_ex 300 (1.0 µm) | dx 1.5 | 5.23 (n 4) | 134 (n 4) | 496 (n 4) |
| d/l_ex 300 (1.0 µm) | dx 3 | 7.73 (n 4) | 188 (n 4) | 619 (n 4) |
| K3: d/l_ex 212, φ 0.69, L_eff 18, r_p 1.5 | dx 1.5 | 4.28 (n 6) | 213 (n 6) | 948 (n 6) |
| K3: d/l_ex 212, φ 0.69, L_eff 18, r_p 1.5 | dx 3 | 2.68 (n 6) | 148 (n 6) | 826 (n 6) |

## Mesh error (paired per seed)

| group | P(dx 3)/P(dx 1.5) 9 mT | 50 mT | 100 mT | Δβ (9 / 50 / 100 mT) | Δβ (9 → 50 mT) |
|---|---|---|---|---|---|
| d/l_ex 212 (0.71 µm) | – | – | – | – | – |
| d/l_ex 300 (1.0 µm) | × 1.48 (ln +0.39 ± 0.20, CI -0.08 … +0.86, n 4) | × 1.40 (ln +0.34 ± 0.16, CI -0.03 … +0.71, n 4) | × 1.25 (ln +0.22 ± 0.06, CI +0.07 … +0.37, n 4) | -0.06 ± 0.05 (CI -0.19 … +0.07, n 4) | -0.02 ± 0.03 (CI -0.09 … +0.05, n 4) |
| K3: d/l_ex 212, φ 0.69, L_eff 18, r_p 1.5 | × 0.63 (ln -0.47 ± 0.09, CI -0.64 … -0.29, n 6) | × 0.70 (ln -0.36 ± 0.13, CI -0.63 … -0.09, n 6) | × 0.87 (ln -0.14 ± 0.09, CI -0.32 … +0.04, n 6) | +0.12 ± 0.03 (CI +0.06 … +0.18, n 6) | +0.05 ± 0.05 (CI -0.04 … +0.14, n 6) |

Stages without a steady value (≥ 2 calm cycles): M212: no runs yet; M300: none; MK3: none

## Mesh share of the size lever (d/l_ex 212 → 300, base microstructure)

Size lever on one mesh: mean ln P (or β) at 300 minus at 212, over the seeds with a value on both meshes (different seeds at 212 and 300, Welch). Mesh share = lever at dx 3 − lever at dx 1.5 = Δ(300) − Δ(212). Sign change of the mesh error: opposite signs and both 90 % CIs exclude 0.

| output | lever at dx 1.5 | lever at dx 3 | mesh share Δ(300) − Δ(212) | sign change |
|---|---|---|---|---|
| lnP 9 | – | – | – | not determined |
| lnP 50 | – | – | – | not determined |
| lnP 100 | – | – | – | not determined |
| beta | – | – | – | not determined |
| beta 9-50 | – | – | – | not determined |

## Old pairs (relaxation not converged; hints only, PROTOKOLL 7.64)

| group | P(dx a)/P(dx 1.5) 9 mT | 50 mT | 100 mT | Δβ (9 / 50 / 100 mT) | Δβ (9 → 50 mT) |
|---|---|---|---|---|---|
| d/l_ex 96, dx a = 3 | n 2: no statement | × 0.33 (ln -1.12 ± 0.67, CI -3.06 … +0.82, n 3) | – | – | n 2: no statement |
| d/l_ex 96, dx a = 2 | n 2: no statement | × 1.30 (ln +0.26 ± 0.33, CI -0.70 … +1.22, n 3) | – | – | n 2: no statement |
| d/l_ex 96, dx a = 1 | – | × 1.64 (ln +0.50 ± 0.36, CI -0.55 … +1.55, n 3) | – | – | – |
| d/l_ex 212, dx a = 3 | × 0.70 (ln -0.36 ± 0.34, CI -1.15 … +0.43, n 4) | × 0.94 (ln -0.06 ± 0.09, CI -0.28 … +0.16, n 4) | – | – | +0.15 ± 0.23 (CI -0.38 … +0.68, n 4) |
