# Design model against the direct coarse-mesh tests of the A effect

A effect = value at A 20 pJ/m minus value at A 10 pJ/m (d = 1 µm), coarse mesh. Model: Gamma GLM, 116 design runs, 1000 bootstrap resamples, 90 % interval. Measured: 90 % interval (t; Welch for the base point).

| point | output | measured | design model at the same factors |
|---|---|---|---|
| K1 pair | beta 9/50/100 | +0.67 (+0.61 … +0.73, n 3, paired) | +0.91 (+0.78 … +1.03) |
| K1 pair | beta 9-150 | +0.64 (+0.53 … +0.76, n 3, paired) | +0.92 (+0.82 … +1.04) |
| K1 pair | lnP 9 | -1.13 (-1.43 … -0.84, n 3, paired) | -1.72 (-2.00 … -1.47) |
| K1 pair | lnP 50 | -0.23 (-0.85 … +0.39, n 3, paired) | -0.29 (-0.56 … -0.10) |
| K1 pair | lnP 100 | +0.55 (+0.39 … +0.72, n 3, paired) | +0.50 (+0.33 … +0.65) |
| base point | beta 9/50/100 | +0.65 (+0.40 … +0.91, n 4 + 4, Welch) | +0.37 (+0.14 … +0.58) |
| base point | beta 9-150 | not measured | +0.41 (+0.21 … +0.59) |
| base point | lnP 9 | -1.00 (-1.67 … -0.34, n 4 + 4, Welch) | -0.68 (-1.15 … -0.05) |
| base point | lnP 50 | +0.09 (-0.23 … +0.42, n 4 + 4, Welch) | -0.13 (-0.54 … +0.63) |
| base point | lnP 100 | +0.59 (+0.37 … +0.82, n 4 + 4, Welch) | +0.23 (-0.07 … +0.63) |

Coded factors (lnQ, r_p, φ, ln d/l_ex): K1 pair A 20 [-1.0, 0.0, 0.0, -0.89], A 10 [-1.0, 0.0, 0.01, 1.07]; base point A 20 [1.0, -1.0, 0.44, -1.02], A 10 [1.0, -1.0, 0.43, 1.01]
