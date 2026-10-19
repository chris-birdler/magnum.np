# Solver test at 150 mT: DemagFieldPBC against the final setting from the same state (eval_solvertest.py)

| run | ln w fd | ln w newell | Δ_s (fd − newell) | ε fd | ε newell | ln w corner run − newell restart |
|---|---|---|---|---|---|---|
| K1_0_nx | -9.5251 | -10.0732 | +0.5482 | +0.0048 | +0.0065 | -0.0000 |
| K1_1_nx | -9.2814 | -9.8444 | +0.5629 | +0.0142 | +0.0189 | +0.0000 |
| K1_2_nx | -9.4171 | -10.0257 | +0.6086 | +0.0097 | +0.0142 | -0.0000 |

Mean Δ_s: +0.573 ± 0.018 (SE, n 3; 90 % CI +0.520 … +0.626). **Reading: solver effect.** For comparison: ln w old − new of the corner runs K1 at 150 mT +0.57 (n 3).
