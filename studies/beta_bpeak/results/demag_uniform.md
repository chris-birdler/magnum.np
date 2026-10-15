# Demag solver tests with uniform m and against the upstream code (check_demag_uniform.py)

Energy per box volume in units of K_d, double precision, local CPU. Relative error = E / E_exact − 1. Newell: one image shell, compensated outer surface charges (check_mesh_energy.newell_pbc).

## 1. B0b geometry, exact φ_vox (1 − φ_vox) / 3

| dx (l_ex) | φ_vox | exact | m along | DemagFieldPBC | rel. error | Newell | rel. error |
|---|---|---|---|---|---|---|---|
| 3 | 0.64757 | 0.076074 | (100) | 0.076074 | -2.00e-08 | 0.076074 | -2.0e-08 |
| 3 | 0.64757 | 0.076074 | (110) | 0.079685 | +4.75e-02 | 0.076069 | -6.9e-05 |
| 3 | 0.64757 | 0.076074 | (111) | 0.083296 | +9.49e-02 | 0.076064 | -1.4e-04 |
| 1.5 | 0.64945 | 0.075888 | (100) | 0.075888 | -2.00e-08 | 0.075888 | -2.0e-08 |
| 1.5 | 0.64945 | 0.075888 | (110) | 0.077685 | +2.37e-02 | 0.075887 | -1.7e-05 |
| 1.5 | 0.64945 | 0.075888 | (111) | 0.079482 | +4.74e-02 | 0.075886 | -3.5e-05 |
| 0.75 | 0.64992 | 0.075841 | (100) | 0.075841 | -2.00e-08 | – | – |
| 0.75 | 0.64992 | 0.075841 | (110) | 0.076742 | +1.19e-02 | – | – |
| 0.75 | 0.64992 | 0.075841 | (111) | 0.077643 | +2.38e-02 | – | – |

## 2. Aligned cube, edge a/2, exact 7/192 = 0.036458

| cells per cube edge | m along | DemagFieldPBC | rel. error | Newell | rel. error |
|---|---|---|---|---|---|
| 4 | (100) | 0.036458 | -2.22e-16 | 0.036458 | +1.3e-14 |
| 4 | (110) | 0.038957 | +6.85e-02 | 0.036458 | +1.3e-14 |
| 4 | (111) | 0.041456 | +1.37e-01 | 0.036458 | +1.3e-14 |
| 8 | (100) | 0.036458 | +2.22e-16 | 0.036458 | +2.4e-12 |
| 8 | (110) | 0.037295 | +2.29e-02 | 0.036458 | +2.4e-12 |
| 8 | (111) | 0.038131 | +4.59e-02 | 0.036458 | +2.4e-12 |
| 16 | (100) | 0.036458 | +0.00e+00 | – | – |
| 16 | (110) | 0.036720 | +7.19e-03 | – | – |
| 16 | (111) | 0.036983 | +1.44e-02 | – | – |
| 32 | (100) | 0.036458 | -3.33e-16 | – | – |
| 32 | (110) | 0.036537 | +2.16e-03 | – | – |
| 32 | (111) | 0.036616 | +4.33e-03 | – | – |

## 3. Our DemagFieldPBC against the upstream code (git fde0c62)

Random m on the B0b geometry, dx 3: max |h_ours − h_upstream| / max |h_upstream| = 5.8e-16.
