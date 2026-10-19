# Bloch-point cut: tube to the particle surface, then relaxation (eval_smooth.py)

Per step: Bloch points before the cut / after the cut / after the relaxation (relaxation iterations; tube lengths in cells).

| run | loaded | step 0 | step 1 | step 2 | after the last relaxation | 9 mT: max (share of samples) | E/K_d loaded → end | ln w 9 mT − corner run |
|---|---|---|---|---|---|---|---|---|
| K1_2_nx | 1 | 1 / 1 / 1 (5000 it; 21) | 1 / 1 / 0 (8000 it; 19) | – | 0 | 0 (0 %) | 0.00174752 → 0.00173203 | +0.074 |
| K2_0_nx | 9 | 9 / 11 / 4 (24000 it; 11, 38, 34, 37, 27, 32, 29, 24, 36) | 4 / 8 / 4 (7500 it; 36, 37, 34, 37) | 4 / 8 / 6 (6000 it; 36, 38, 34, 37) | 6 | 6 (100 %) | -0.00743144 → -0.00744436 | -0.205 |
| K2_1_nx | 5 | 5 / 6 / 5 (13500 it; 2, 36, 24, 14, 25) | 5 / 6 / 2 (17500 it; 1, 37, 26, 13, 26) | 2 / 2 / 1 (13500 it; 37, 25) | 1 | 1 (100 %) | -0.00743489 → -0.00748 | +0.091 |
| K2_2_nx | – | – | – | – | not done | – | – | – |
| K3_1_nx | – | – | – | – | not done | – | – | – |
| K3_2_nx | – | – | – | – | not done | – | – | – |
| K4_0_nx | – | – | – | – | not done | – | – | – |
| K4_1_nx | – | – | – | – | not done | – | – | – |
| K4_2_nx | – | – | – | – | not done | – | – | – |

Clean: 1 of 3 runs done (K1_2_nx).

**Reading: none: test stopped by Chris (2026-10-19) after 3 of 9 runs.**

## Test after one AC pass (option C with one attempt, then the cut)

| run | after the AC pass (option C run, attempt 0) | step 0 | step 1 | step 2 | after the last relaxation | 9 mT: max (share of samples) | E/K_d after the AC pass → end | ln w 9 mT − option C run |
|---|---|---|---|---|---|---|---|---|
| K3_1_ad | 2 (2) | 2 / 2 / 2 (2500 it; 8, 8) | 2 / 3 / 2 (3000 it; 8, 8) | 2 / 1 / 2 (9000 it; 8, 8) | 2 | 2 (100 %) | 0.000715846 → 0.000716045 | +0.375 |
| K4_0_ad | 1 (1) | 1 / 3 / 1 (2500 it; 47) | 1 / 1 / 1 (4500 it; 46) | 1 / 1 / 1 (10500 it; 44) | 1 | 1 (100 %) | -0.000352049 → -0.000351839 | -0.373 |
| K2_2_ad | 1 (1) | 1 / 1 / 1 (2000 it; 30) | 1 / 1 / 1 (3500 it; 30) | 1 / 2 / 1 (15000 it; 29) | 1 | 1 (100 %) | -0.00749104 → -0.00749103 | +0.085 |

Tested (Bloch point after the AC pass): 3 of 3 runs done (K3_1_ad, K4_0_ad, K2_2_ad); clean: 0 (none).

**Reading: does not work.**
