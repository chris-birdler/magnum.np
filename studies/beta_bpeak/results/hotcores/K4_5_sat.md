# Hot-spot follow-up: structure per snapshot, vortex core or Bloch line, Thiele scale (eval_hotcores.py)

Main regions of the loss outside f as in eval_hotspots.py. line / wall / 3D: shares of the snapshots of each type; turn: angle of the line direction of a snapshot to the mean one (median / max). centre: |m.t| of the core cell (information); ring: (max - min) / 2 of m.t at 2 ... 3 cells from the core (medians); class from the ring (threshold 0.30) and the type shares. Thiele: loss of the window from the accumulators over the loss of the rigid motion along the core path, P (all) and q (P minus the loss at f); snap / acc: loss at f and at 2f ... 7f from the 16 snapshot phases over the same from the accumulators; v harm: share of k >= 2 in the mean square core velocity; step: largest core step between 2 phases (cells); amp: core path, l_ex.

| run | stage | share | particle | line | wall | 3D | turn, deg | centre | ring | class | amp, l_ex | step | jumps | v harm | P acc / rigid | q acc / rigid | f snap / acc | 2f ... 7f snap / acc |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| K4_5_sat | amp00 | 0.661 | 0 | 1.00 | 0.00 | 0.00 | 1 / 1 | 0.98 | 0.11 | vortex core | 8.4 | 1.3 | 0.00 | 0.73 | 0.21 | 0.02 | 1.00 | 1.02 |
| K4_5_sat | amp03 | 0.166 | 0 | 1.00 | 0.00 | 0.00 | 8 / 11 | 0.98 | 0.12 | vortex core | 39.1 | 3.6 | 0.02 | 0.28 | 0.85 | 0.97 | 1.00 | 1.06 |
| K4_5_sat | amp03 | 0.089 | 3 | 1.00 | 0.00 | 0.00 | 8 / 12 | 0.92 | 0.25 | vortex core | 18.9 | 2.3 | 0.03 | 0.55 | 1.24 | 1.27 | 1.01 | 1.12 |
| K4_5_sat | amp06 | 0.107 | 3 | 1.00 | 0.00 | 0.00 | 20 / 28 | 0.86 | 0.24 | vortex core | 18.2 | 5.8 | 0.95 | 0.91 | 2.42 | 1.80 | 0.98 | 1.01 |
| K4_5_sat | amp06 | 0.104 | 1 | 1.00 | 0.00 | 0.00 | 13 / 35 | 0.45 | 0.62 | Bloch line | 25.9 | 4.1 | 0.52 | 0.70 | 3.00 | 3.31 | 1.00 | 1.09 |
