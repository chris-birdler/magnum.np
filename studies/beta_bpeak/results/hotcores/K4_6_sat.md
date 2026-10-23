# Hot-spot follow-up: structure per snapshot, vortex core or Bloch line, Thiele scale (eval_hotcores.py)

Main regions of the loss outside f as in eval_hotspots.py. line / wall / 3D: shares of the snapshots of each type; turn: angle of the line direction of a snapshot to the mean one (median / max). centre: |m.t| of the core cell (information); ring: (max - min) / 2 of m.t at 2 ... 3 cells from the core (medians); class from the ring (threshold 0.30) and the type shares. Thiele: loss of the window from the accumulators over the loss of the rigid motion along the core path, P (all) and q (P minus the loss at f); snap / acc: loss at f and at 2f ... 7f from the 16 snapshot phases over the same from the accumulators; v harm: share of k >= 2 in the mean square core velocity; step: largest core step between 2 phases (cells); amp: core path, l_ex.

| run | stage | share | particle | line | wall | 3D | turn, deg | centre | ring | class | amp, l_ex | step | jumps | v harm | P acc / rigid | q acc / rigid | f snap / acc | 2f ... 7f snap / acc |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| K4_6_sat | amp00 | 0.566 | 3 | 0.00 | 0.00 | 1.00 | 1 / 1 | 0.25 | 0.57 | unclear | 3.0 | 0.8 | 0.16 | 0.78 | 0.05 | 0.01 | 1.00 | 1.13 |
| K4_6_sat | amp03 | 0.369 | 0 | 0.00 | 0.00 | 1.00 | 12 / 73 | 0.24 | 0.50 | unclear | 13.8 | 2.8 | 0.47 | 0.76 | 0.86 | 0.40 | 0.99 | 1.15 |
| K4_6_sat | amp06 | 0.263 | 0 | 0.88 | 0.00 | 0.12 | 20 / 68 | 0.65 | 0.33 | Bloch line | 24.0 | 3.5 | 0.61 | 0.74 | 3.82 | 4.56 | 1.00 | 1.04 |
