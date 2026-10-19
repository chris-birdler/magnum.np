# Corrugation of the vortex energy by the surface (check_mesh_energy.py corrugation, d/l_ex 48, φ 0.45)

E(δ): total energy per box volume in K_d with the vortex axes shifted by δ = 0 … 2 dx (16 steps). Corrugation = peak to peak of E(δ) minus its quadratic fit. Reference: anisotropy energy of the same field.

| surface | dx (l_ex) | E(0) | E(2 dx) − E(0) | corrugation | corrugation / abs. anisotropy energy |
|---|---|---|---|---|---|
| voxel | 3 | 0.013229 | +4.04e-03 | 4.31e-05 | 0.042 |
| voxel | 1.5 | 0.009798 | +1.07e-03 | 3.41e-05 | 0.033 |
| voxel | 0.75 | 0.008005 | +2.75e-04 | 2.44e-06 | 0.002 |
| fraction | 3 | 0.008999 | +3.64e-03 | 4.70e-05 | 0.046 |
| fraction | 1.5 | 0.007435 | +1.01e-03 | 3.42e-05 | 0.033 |
| fraction | 0.75 | 0.006814 | +2.72e-04 | 2.32e-06 | 0.002 |
