# Hot-spot follow-up: structure per snapshot, vortex core or Bloch line, Thiele scale (eval_hotcores.py)

Main regions of the loss outside f as in eval_hotspots.py. line / wall / 3D: shares of the snapshots of each type; turn: angle of the line direction of a snapshot to the mean one (median / max). centre: |m.t| of the core cell (information); ring: (max - min) / 2 of m.t at 2 ... 3 cells from the core (medians); class from the ring (threshold 0.30) and the type shares. Thiele: loss of the window from the accumulators over the loss of the rigid motion along the core path, P (all) and q (P minus the loss at f); snap / acc: loss at f and at 2f ... 7f from the 16 snapshot phases over the same from the accumulators; v harm: share of k >= 2 in the mean square core velocity; step: largest core step between 2 phases (cells); amp: core path, l_ex.

| run | stage | share | particle | line | wall | 3D | turn, deg | centre | ring | class | amp, l_ex | step | jumps | v harm | P acc / rigid | q acc / rigid | f snap / acc | 2f ... 7f snap / acc |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| K4_8_sat | amp00 | 0.322 | 0 | 1.00 | 0.00 | 0.00 | 0 / 1 | 0.99 | 0.13 | vortex core | 6.6 | 1.0 | 0.00 | 0.79 | 0.19 | 0.03 | 1.03 | 0.95 |
| K4_8_sat | amp03 | 0.457 | 0 | 1.00 | 0.00 | 0.00 | 4 / 5 | 0.99 | 0.11 | vortex core | 22.0 | 2.8 | 0.00 | 0.45 | 1.05 | 1.07 | 0.97 | 1.18 |
| K4_8_sat | amp06 | 0.214 | 3 | 0.81 | 0.00 | 0.19 | 20 / 81 | 0.64 | 0.35 | Bloch line | 30.0 | 4.6 | 0.72 | 0.69 | 3.99 | 5.10 | 1.00 | 1.07 |
| K4_8_sat | amp06 | 0.183 | 0 | 0.62 | 0.00 | 0.38 | 59 / 88 | 0.48 | 0.41 | Bloch line | 7.5 | 2.5 | 0.86 | 0.89 | 5.04 | 5.12 | 0.58 | 0.90 |
