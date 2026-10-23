# Hot-spot follow-up: structure per snapshot, vortex core or Bloch line, Thiele scale (eval_hotcores.py)

Main regions of the loss outside f as in eval_hotspots.py. line / wall / 3D: shares of the snapshots of each type; turn: angle of the line direction of a snapshot to the mean one (median / max). centre: |m.t| of the core cell (information); ring: (max - min) / 2 of m.t at 2 ... 3 cells from the core (medians); class from the ring (threshold 0.30) and the type shares. Thiele: loss of the window from the accumulators over the loss of the rigid motion along the core path, P (all) and q (P minus the loss at f); snap / acc: loss at f and at 2f ... 7f from the 16 snapshot phases over the same from the accumulators; v harm: share of k >= 2 in the mean square core velocity; step: largest core step between 2 phases (cells); amp: core path, l_ex.

| run | stage | share | particle | line | wall | 3D | turn, deg | centre | ring | class | amp, l_ex | step | jumps | v harm | P acc / rigid | q acc / rigid | f snap / acc | 2f ... 7f snap / acc |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| K4_9_sat | amp00 | 0.612 | 3 | 1.00 | 0.00 | 0.00 | 0 / 0 | 0.99 | 0.09 | vortex core | 6.2 | 1.0 | 0.00 | 0.76 | 0.21 | 0.01 | 1.01 | 1.00 |
| K4_9_sat | amp03 | 0.340 | 0 | 0.81 | 0.00 | 0.19 | 6 / 11 | 0.96 | 0.17 | vortex core | 32.2 | 3.0 | 0.00 | 0.37 | 1.08 | 1.10 | 1.02 | 1.10 |
| K4_9_sat | amp06 | 0.288 | 0 | 0.88 | 0.00 | 0.12 | 16 / 56 | 0.91 | 0.22 | vortex core | 23.5 | 4.8 | 0.97 | 0.87 | 2.34 | 2.19 | 0.93 | 1.16 |
