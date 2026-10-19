# Cost of the cycle tolerance and of the final stray-field setting on dx 3 (eval_tolbench.py)

| run | wall per cycle (s) | steps per cycle | mean step (ps) | ε (cycles 1, 2) | P (W/cm³) | GPU |
|---|---|---|---|---|---|---|
| TB_K1_0_a5 | 104 (+0 %) | 4608 (+0 %) | 7.234 | -0.053 / +0.035 | 2.705 | Tesla V100-SXM2-16GB |
| TB_K1_0_a6 | 120 (+15 %) | 4672 (+1 %) | 7.135 | -0.092 / +0.006 | 2.705 | Tesla V100-SXM2-16GB |
| TB_K1_0_tight | 111 (+7 %) | 5272 (+14 %) | 6.322 | -0.092 / +0.006 | 2.705 | Tesla V100-SXM2-16GB |
| K1_0, old run (DemagFieldPBC, 9 mT stage) | – | – | – | – | – | Tesla V100-SXM2-16GB |

K1_0: wall time per step 15.85 ms (DemagFieldPBC, old run) against 22.57 ms (final setting, a5): +42 %.

| TB_K2_0_a5 | 273 (+0 %) | 4839 (+0 %) | 6.888 | -0.059 / +0.003 | 4.496 | Tesla V100-SXM2-16GB |
| TB_K2_0_a6 | 274 (+0 %) | 4835 (-0 %) | 6.894 | -0.061 / +0.001 | 4.496 | Tesla V100-SXM2-16GB |
| TB_K2_0_tight | 238 (-13 %) | 5453 (+13 %) | 6.113 | -0.061 / +0.001 | 4.496 | Tesla V100-SXM2-16GB |
| K2_0, old run (DemagFieldPBC, 9 mT stage) | – | – | – | – | – | Tesla V100-SXM2-16GB |

K2_0: wall time per step 39.73 ms (DemagFieldPBC, old run) against 56.32 ms (final setting, a5): +42 %.

