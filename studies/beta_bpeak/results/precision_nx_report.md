# Precision test of PLAN 4.13: fp64 against fp32 in the drive cycles (eval_precision.py)

9 mT stage from the same init.pt; Δ_p = ln w(p64) − ln w(p32) on the calm tails.

| run | ln w p64 | ln w p32 | Δ_p | ε p64 | ε p32 | wall per cycle p64 / p32 (s) | ln w corner run − p32 |
|---|---|---|---|---|---|---|---|
| K1_0_nx | -16.1526 | -16.1526 | -0.0000 | +0.0063 | +0.0065 | 241 / 155 | +0.0000 |
| K2_0_nx | -15.6166 | -15.6167 | +0.0000 | +0.0002 | +0.0003 | 1030 / 333 | +0.0000 |
| K3_0_nx | -16.5499 | -16.5500 | +0.0001 | +0.0039 | +0.0041 | 239 / 150 | +0.0000 |
| K4_0_nx | -15.5691 | -15.5690 | -0.0000 | +0.0109 | +0.0110 | 1202 / 332 | +0.0000 |

Mean Δ_p: +0.0000 ± 0.0000 (SE, n 4; 90 % CI -0.0000 … +0.0001). **Reading: fp32 sufficient.**
