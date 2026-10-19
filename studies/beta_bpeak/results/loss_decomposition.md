# Loss decomposition P = a1·x + a2·x² + a3·x³ (x = B / 50 mT), 116 design runs

Steady loops (≥ 2 calm cycles) with reruns, coarse mesh, P in W/cm³ per core volume at 30 MHz. Non-negative fit of the relative residual. Term names are hypotheses: x¹ switching, x² linear viscous, x³ Rayleigh.

## Fit quality

- rms of the ln P residual: 3 terms median 0.137 (quartiles 0.105 … 0.204); single power law median 0.134 (0.086 … 0.173).
- Terms present (a_j > 0) / determined (a_j > 2 SE): x¹ switching 54 / 17, x² viscous 115 / 91, x³ Rayleigh 80 / 38.

## Share of each term in P (median over runs; quartiles)

| group | term | 9 mT | 50 mT | 150 mT |
|---|---|---|---|---|
| all | x¹ switching | 0.00 (0.00 … 0.33) | 0.00 (0.00 … 0.06) | 0.00 (0.00 … 0.02) |
| all | x² viscous | 0.87 (0.67 … 0.96) | 0.80 (0.62 … 0.91) | 0.61 (0.36 … 0.94) |
| all | x³ Rayleigh | 0.03 (0.00 … 0.10) | 0.17 (0.00 … 0.37) | 0.39 (0.00 … 0.64) |
| low A (10 … 12.5 pJ/m) | x¹ switching | 0.13 (0.00 … 0.50) | 0.03 (0.00 … 0.15) | 0.01 (0.00 … 0.06) |
| low A (10 … 12.5 pJ/m) | x² viscous | 0.87 (0.50 … 0.96) | 0.86 (0.75 … 0.93) | 0.93 (0.68 … 0.97) |
| low A (10 … 12.5 pJ/m) | x³ Rayleigh | 0.00 (0.00 … 0.02) | 0.00 (0.00 … 0.11) | 0.00 (0.00 … 0.27) |
| high A (15.8 … 20 pJ/m) | x¹ switching | 0.00 (0.00 … 0.10) | 0.00 (0.00 … 0.01) | 0.00 (0.00 … 0.00) |
| high A (15.8 … 20 pJ/m) | x² viscous | 0.85 (0.75 … 0.90) | 0.55 (0.47 … 0.66) | 0.30 (0.23 … 0.39) |
| high A (15.8 … 20 pJ/m) | x³ Rayleigh | 0.11 (0.08 … 0.14) | 0.43 (0.34 … 0.52) | 0.70 (0.60 … 0.76) |

## Factor effects on the shares (OLS, HC3; range effect = 2 × coefficient; A 10 → 20 pJ/m is − lnd)

| share | B | K_eff (lnQ) | r_p | φ | A (− lnd) | lnd × lnQ | R² |
|---|---|---|---|---|---|---|---|
| x¹ switching | 9 mT | -0.14 ± 0.06 | -0.14 ± 0.08 | +0.09 ± 0.07 | -0.26 ± 0.07 | -0.32 ± 0.07 | 0.31 |
| x¹ switching | 50 mT | -0.06 ± 0.02 | -0.06 ± 0.03 | +0.04 ± 0.02 | -0.11 ± 0.02 | -0.11 ± 0.02 | 0.41 |
| x¹ switching | 150 mT | -0.02 ± 0.01 | -0.03 ± 0.01 | +0.01 ± 0.01 | -0.05 ± 0.01 | -0.04 ± 0.01 | 0.47 |
| x² viscous | 9 mT | +0.14 ± 0.06 | +0.13 ± 0.08 | -0.08 ± 0.07 | +0.13 ± 0.07 | +0.22 ± 0.08 | 0.18 |
| x² viscous | 50 mT | -0.01 ± 0.05 | +0.02 ± 0.07 | -0.02 ± 0.05 | -0.34 ± 0.06 | -0.17 ± 0.07 | 0.34 |
| x² viscous | 150 mT | -0.13 ± 0.06 | -0.08 ± 0.08 | +0.00 ± 0.06 | -0.61 ± 0.07 | -0.39 ± 0.09 | 0.54 |
| x³ Rayleigh | 9 mT | +0.00 ± 0.01 | +0.01 ± 0.01 | -0.00 ± 0.01 | +0.13 ± 0.01 | +0.10 ± 0.02 | 0.55 |
| x³ Rayleigh | 50 mT | +0.07 ± 0.04 | +0.05 ± 0.06 | -0.02 ± 0.04 | +0.45 ± 0.05 | +0.28 ± 0.06 | 0.56 |
| x³ Rayleigh | 150 mT | +0.15 ± 0.06 | +0.11 ± 0.08 | -0.01 ± 0.06 | +0.66 ± 0.07 | +0.43 ± 0.09 | 0.58 |

## Damping check (centre C0 … C7, steady loops ≥ 2 at both α)

| seed | amplitudes | a1 α0.1 → α0.01 | a2 α0.1 → α0.01 | a3 α0.1 → α0.01 |
|---|---|---|---|---|
| C0 | 4 | 25.3 → 102.4 | 184.6 → 81.6 | 0.0 → 0.0 |
| C1 | 5 | 69.8 → 0.0 | 65.9 → 65.3 | 22.4 → 18.7 |
| C3 | 4 | 0.0 → 32.5 | 145.8 → 206.5 | 58.9 → 0.0 |
| C4 | 4 | 15.2 → 43.2 | 174.4 → 82.2 | 0.0 → 0.0 |

Sum over the seeds, ratio α 0.01 / α 0.1: a1 × 1.61, a2 × 0.76, a3 × 0.23 (viscous expectation for a2: × 0.1).
