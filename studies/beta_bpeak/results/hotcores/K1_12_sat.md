# Hot-spot follow-up: structure per snapshot, vortex core or Bloch line, Thiele scale (eval_hotcores.py)

Main regions of the loss outside f as in eval_hotspots.py. line / wall / 3D: shares of the snapshots of each type; turn: angle of the line direction of a snapshot to the mean one (median / max). centre: |m.t| of the core cell (information); ring: (max - min) / 2 of m.t at 2 ... 3 cells from the core (medians); class from the ring (threshold 0.30) and the type shares. Thiele: loss of the window from the accumulators over the loss of the rigid motion along the core path, P (all) and q (P minus the loss at f); snap / acc: loss at f and at 2f ... 7f from the 16 snapshot phases over the same from the accumulators; v harm: share of k >= 2 in the mean square core velocity; step: largest core step between 2 phases (cells); amp: core path, l_ex.

| run | stage | share | particle | line | wall | 3D | turn, deg | centre | ring | class | amp, l_ex | step | jumps | v harm | P acc / rigid | q acc / rigid | f snap / acc | 2f ... 7f snap / acc |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| K1_12_sat | amp00 | 0.181 | 3 | 1.00 | 0.00 | 0.00 | 0 / 0 | 0.93 | 0.16 | vortex core | 0.0 | 0.0 | 0.00 | 0.18 | 0.05 | 0.00 | 1.00 | 0.99 |
| K1_12_sat | amp03 | 0.416 | 0 | 1.00 | 0.00 | 0.00 | 1 / 1 | 0.96 | 0.15 | vortex core | 11.8 | 1.3 | 0.00 | 0.59 | 0.26 | 0.02 | 1.00 | 1.01 |
| K1_12_sat | amp06 | 0.429 | 0 | 0.94 | 0.00 | 0.06 | 17 / 29 | 0.92 | 0.17 | vortex core | 59.4 | 5.5 | 0.70 | 0.18 | 1.12 | 1.31 | 1.08 | 1.03 |
