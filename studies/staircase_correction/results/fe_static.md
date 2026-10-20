# Static FE comparison: coarse m in the smooth spheres (PROTOKOLL 3.24)

FE: magnum.pi 3.13.3, P2 potential, periodic box mesh of n^3 cubes with 6 tetrahedra each, mean field zero. FD: magnum.np, Newell with full PBC. Values in units of S = E_vox − E_X3,q2 of the state (column S, in K_d V_box). E_s: energy of the smooth spheres from the difference method.

Checks of the difference method: ε_q2, ε_q4 on the 5-nm and 2.5-nm staircases (exact FD values); ε_u on the smooth spheres with uniform m (exact target φ (1 − φ) / 3). Continuous m: the raw FE value without a correction; its FE error is not known; the raw FE error of the smooth spheres with uniform m is given for scale.

| state | mesh | S (K_d V_box) | ε_q2 | ε_q4 | ε_u | X3 q2 − E_s | X3 q4 − E_s | T3Pc q2 − E_s | T3Pc q4 − E_s | continuous m (raw) − E_s | raw FE error, smooth, uniform m |
|---|---|---|---|---|---|---|---|---|---|---|---|
| T3Pc_s941 | box:48 | 2.14e-03 | -0.403 | -0.115 | – | +0.666 | +0.146 | +0.674 | +0.158 | – | – |
| T3Pc_s942 | box:48 | 2.14e-03 | -0.406 | -0.117 | – | +0.671 | +0.148 | +0.679 | +0.159 | – | – |
| T3Pc_s943 | box:48 | 2.14e-03 | – | – | – | – | – | – | – | – | – |

Cost per FE run (cells, wall time in s, peak RAM in GB):

- T3Pc_s941 box:48: vox/state 663552 cells, 60 s, 4.4 GB; full/state 663552 cells, 60 s, 4.4 GB; q2/state 663552 cells, 122 s, 4.4 GB; q4/state 663552 cells, 122 s, 4.4 GB; smooth/state 663552 cells, 1220 s, 4.4 GB
- T3Pc_s942 box:48: vox/state 663552 cells, 59 s, 4.4 GB; full/state 663552 cells, 59 s, 4.4 GB; q2/state 663552 cells, 122 s, 4.4 GB; q4/state 663552 cells, 124 s, 4.4 GB; smooth/state 663552 cells, 1213 s, 4.4 GB
- T3Pc_s943 box:48: vox/state 663552 cells, 59 s, 4.4 GB; full/state 663552 cells, 60 s, 4.4 GB
