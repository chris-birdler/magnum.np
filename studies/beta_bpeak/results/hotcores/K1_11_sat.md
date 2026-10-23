# Hot-spot follow-up: structure per snapshot, vortex core or Bloch line, Thiele scale (eval_hotcores.py)

Main regions of the loss outside f as in eval_hotspots.py. line / wall / 3D: shares of the snapshots of each type; turn: angle of the line direction of a snapshot to the mean one (median / max). centre: |m.t| of the core cell (information); ring: (max - min) / 2 of m.t at 2 ... 3 cells from the core (medians); class from the ring (threshold 0.30) and the type shares. Thiele: loss of the window from the accumulators over the loss of the rigid motion along the core path, P (all) and q (P minus the loss at f); snap / acc: loss at f and at 2f ... 7f from the 16 snapshot phases over the same from the accumulators; v harm: share of k >= 2 in the mean square core velocity; step: largest core step between 2 phases (cells); amp: core path, l_ex.

| run | stage | share | particle | line | wall | 3D | turn, deg | centre | ring | class | amp, l_ex | step | jumps | v harm | P acc / rigid | q acc / rigid | f snap / acc | 2f ... 7f snap / acc |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| K1_11_sat | amp00 | 0.079 | 0 | 1.00 | 0.00 | 0.00 | 0 / 0 | 0.95 | 0.11 | vortex core | 0.0 | 0.0 | 0.00 | 0.02 | 0.01 | 0.00 | 1.00 | 0.98 |
| K1_11_sat | amp03 | 0.292 | 3 | 1.00 | 0.00 | 0.00 | 1 / 1 | 0.91 | 0.25 | vortex core | 2.9 | 0.9 | 0.00 | 0.85 | 0.07 | 0.03 | 1.00 | 1.02 |
| K1_11_sat | amp06 | 0.165 | 0 | 1.00 | 0.00 | 0.00 | 1 / 2 | 0.60 | 0.17 | vortex core | 6.4 | 2.0 | 0.17 | 0.86 | 0.60 | 0.21 | 1.00 | 1.16 |
| K1_11_sat | amp06 | 0.099 | 3 | 1.00 | 0.00 | 0.00 | 2 / 3 | 0.93 | 0.15 | vortex core | 17.0 | 2.8 | 0.20 | 0.69 | 0.19 | 0.05 | 1.00 | 1.02 |
