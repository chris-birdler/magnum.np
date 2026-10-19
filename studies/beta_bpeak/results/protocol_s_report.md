# Protocol S: start saturated, AC from saturation, AC from 300 mT (eval_protocol_s.py)

Per pass: h0 (H/M_s), half cycles, Bloch points max during the pass / after its relaxation, relaxation iterations, time (h). Pass 0: from saturation; pass 1: from 300 mT; pass 2: the further pass from 300 mT.

| run | pass 0 | pass 1 | pass 2 | BP after the last relaxation | remanence M/M_s | BP max at 9 mT | ln w 9 mT − corner run | ln w 9 mT − AC-300 test |
|---|---|---|---|---|---|---|---|---|
| K2_0_nx | 0.38, 34, 628 / 1, 2500, 2.01 | 0.03093, 23, 12 / 1, 2500, 0.90 | 0.03093, 23, 12 / 1, 5500, 0.91 | 1 | 0.0018 | 1 | -0.111 | -0.034 |
| K2_1_nx | 0.38, 34, 677 / 0, 2500, 1.95 | 0.03103, 23, 11 / 0, 6500, 0.89 | – | 0 | 0.0013 | 0 | +0.903 | +0.836 |
| K4_1_nx | 0.45, 35, 1362 / 0, 3500, 1.88 | 0.04009, 24, 1 / 0, 2000, 0.94 | – | 0 | 0.0003 | 0 | +0.191 | -0.919 |
| K3_1_nx | 0.31, 33, 200 / 0, 1500, 0.71 | 0.01886, 20, 1 / 0, 3000, 0.34 | – | 0 | 0.0021 | 0 | -0.872 | -0.352 |

**Reading: not clean.** K2_0_nx: 1 Bloch points after the last relaxation; K2_0_nx: Bloch points at 9 mT
