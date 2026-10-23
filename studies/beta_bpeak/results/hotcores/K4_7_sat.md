# Hot-spot follow-up: structure per snapshot, vortex core or Bloch line, Thiele scale (eval_hotcores.py)

Main regions of the loss outside f as in eval_hotspots.py. line / wall / 3D: shares of the snapshots of each type; turn: angle of the line direction of a snapshot to the mean one (median / max). centre: |m.t| of the core cell (information); ring: (max - min) / 2 of m.t at 2 ... 3 cells from the core (medians); class from the ring (threshold 0.30) and the type shares. Thiele: loss of the window from the accumulators over the loss of the rigid motion along the core path, P (all) and q (P minus the loss at f); snap / acc: loss at f and at 2f ... 7f from the 16 snapshot phases over the same from the accumulators; v harm: share of k >= 2 in the mean square core velocity; step: largest core step between 2 phases (cells); amp: core path, l_ex.

| run | stage | share | particle | line | wall | 3D | turn, deg | centre | ring | class | amp, l_ex | step | jumps | v harm | P acc / rigid | q acc / rigid | f snap / acc | 2f ... 7f snap / acc |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| K4_7_sat | amp00 | 0.560 | 0 | 1.00 | 0.00 | 0.00 | 0 / 0 | 0.96 | 0.15 | vortex core | 6.0 | 1.1 | 0.00 | 0.77 | 0.17 | 0.01 | 1.02 | 0.98 |
| K4_7_sat | amp03 | 0.551 | 0 | 1.00 | 0.00 | 0.00 | 2 / 4 | 0.97 | 0.13 | vortex core | 14.4 | 2.1 | 0.00 | 0.67 | 0.67 | 0.53 | 0.97 | 1.14 |
| K4_7_sat | amp06 | 0.163 | 0 | 0.81 | 0.00 | 0.19 | 26 / 39 | 0.92 | 0.20 | vortex core | 33.9 | 4.8 | 0.56 | 0.44 | 1.52 | 1.72 | 1.04 | 1.12 |
| K4_7_sat | amp06 | 0.123 | 3 | 1.00 | 0.00 | 0.00 | 2 / 5 | 0.11 | 0.11 | vortex core | 0.2 | 0.0 | 0.17 | 0.40 | 94.92 | 80.01 | 0.13 | 0.04 |
