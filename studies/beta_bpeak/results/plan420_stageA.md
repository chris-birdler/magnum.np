# PLAN 4.20 stage A: relaxed states of the tied variants (staircase geometry 48^3, q 8)

Pairs above 90 deg by type (ff full-full, cf cut-full, cc cut-cut; thin: with a thin cell; tied: with a tied cell). Bloch points (+, -) of the relaxed state. Overlap: mean m . m_B8 per particle.

| arm | seed | conv. | iterations | tied | free cut tilt max (deg) | free cut > 45 deg | pairs > 90 deg (all / ff / cf / cc / thin / tied) | largest pair angle (deg) | Bloch | overlap with B8 |
|---|---|---|---|---|---|---|---|---|---|---|
| B8 | 941 | True | 3500 | 0 | 72.9 | 9 | 5 / 0 / 1 / 4 / 5 / 0 | 133.7 | [0, 0] | - |
| B8 | 942 | True | 9000 | 0 | 75.4 | 11 | 4 / 0 / 1 / 3 / 4 / 0 | 124.2 | [0, 0] | - |
| B8 | 943 | True | 6000 | 0 | 76.8 | 6 | 2 / 0 / 0 / 2 / 2 / 0 | 102.1 | [0, 0] | - |
| B8 | 944 | True | 10000 | 0 | 82.1 | 12 | 6 / 0 / 2 / 4 / 6 / 0 | 118.7 | [0, 0] | - |
| C8t | 941 | True | 6000 | 0 | 25.6 | 0 | 0 / 0 / 0 / 0 / 0 / 0 | 74.5 | [0, 0] | - |
| C8t | 942 | True | 12500 | 0 | 23.2 | 0 | 0 / 0 / 0 / 0 / 0 / 0 | 75.0 | [0, 0] | - |
| C8t | 943 | True | 10500 | 0 | 22.6 | 0 | 0 / 0 / 0 / 0 / 0 / 0 | 69.7 | [0, 0] | - |
| C8t | 944 | True | 7000 | 0 | 22.1 | 0 | 0 / 0 / 0 / 0 / 0 / 0 | 71.2 | [0, 0] | - |
| V1 | 941 | True | 13000 | 7104 | 40.9 | 0 | 1 / 0 / 0 / 1 / 0 / 0 | 92.5 | [0, 0] | 0.9997 0.9998 0.9995 0.9996 |
| V1 | 942 | True | 4000 | 7104 | 38.4 | 0 | 1 / 0 / 0 / 1 / 0 / 0 | 90.4 | [0, 0] | 0.9997 0.9997 0.9994 0.9995 |
| V1 | 943 | True | 15000 | 7104 | 40.1 | 0 | 0 / 0 / 0 / 0 / 0 / 0 | 86.7 | [0, 0] | 0.9996 0.9997 0.9997 0.9997 |
| V1 | 944 | True | 13000 | 7104 | 39.4 | 0 | 0 / 0 / 0 / 0 / 0 / 0 | 89.3 | [0, 0] | 0.9996 0.9994 0.9997 0.9993 |
| V2 | 941 | True | 11500 | 8256 | 42.2 | 0 | 1 / 0 / 0 / 1 / 0 / 0 | 90.0 | [0, 0] | 0.9998 0.9999 0.9996 0.9997 |
| V2 | 942 | True | 13000 | 8256 | 38.5 | 0 | 0 / 0 / 0 / 0 / 0 / 0 | 89.5 | [0, 0] | 0.9998 0.9997 0.9995 0.9997 |
| V2 | 943 | True | 6500 | 8256 | 40.5 | 0 | 0 / 0 / 0 / 0 / 0 / 0 | 85.8 | [0, 0] | 0.9997 0.9998 0.9998 0.9997 |
| V2 | 944 | True | 27000 | 8256 | 39.3 | 0 | 0 / 0 / 0 / 0 / 0 / 0 | 89.3 | [0, 0] | 0.9996 0.9995 0.9998 0.9995 |

## Gate A (a variant passes with no pair above 90 deg and no Bloch point in at least 3 of 4 seeds)

- V1 (thin_face): seeds without a pair above 90 deg and a Bloch point: [943, 944]; fails
- V2 (ghost_mirror): seeds without a pair above 90 deg and a Bloch point: [942, 943, 944]; passes
