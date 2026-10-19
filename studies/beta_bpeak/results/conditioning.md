# Conditioning test (PLAN 4.5 D): drift of the loss density over 30 cycles

Centre point, seeds C0 … C5. Drift d = ln P(cycles 28 … 29) − ln P(design estimator: calm tail of the cycles ≤ 6). Decision: 90 % CI of the mean inside ±0.05 → design value stands; 95 % CI excludes 0 and |mean| ≥ 0.05 → systematic error; else not determined (t, n − 1).

## Drift per seed and amplitude (events = cycles 2 … 29 with |w_dis/w_loop − 1| > 0.20)

| seed | 9 mT: d (events) | 20 mT: d (events) | 35 mT: d (events) |
|---|---|---|---|
| C0 | +0.051 (0) | -0.006 (0) | -0.012 (0) |
| C1 | -0.030 (2: 2, 3) | -0.019 (0) | -0.001 (0) |
| C2 | +0.004 (0) | +0.008 (0) | -0.034 (1: 4) |
| C3 | +0.008 (0) | -0.002 (0) | -0.007 (0) |
| C4 | +0.008 (0) | -0.008 (0) | +0.034 (0) |
| C5 | -0.035 (1: 2) | +0.018 (2: 3, 4) | -0.003 (0) |

## Decision per amplitude

| B (mT) | n | mean d ± SE | 90 % CI | 95 % CI | factor P(late)/P(design) | decision |
|---|---|---|---|---|---|---|
| 9 | 6 | +0.001 ± 0.013 | -0.024 … +0.027 | -0.031 … +0.034 | × 1.001 | design value stands |
| 20 | 6 | -0.001 ± 0.005 | -0.012 … +0.009 | -0.015 … +0.012 | × 0.999 | design value stands |
| 35 | 6 | -0.004 ± 0.009 | -0.022 … +0.014 | -0.027 … +0.019 | × 0.996 | design value stands |

## β from the design estimators and from the late values (per seed, mean ± SE)

| window | n | β design | β late | late − design |
|---|---|---|---|---|
| 9 → 20 mT | 6 | 2.127 ± 0.159 | 2.124 ± 0.156 | -0.003 ± 0.018 |
| 20 → 35 mT | 6 | 2.021 ± 0.262 | 2.016 ± 0.273 | -0.004 ± 0.021 |

Determinism: the earlier stages of cond_Ck_b20 / _b35 repeat the design stages of Ck; the gate check_reruns.py compared them cycle by cycle (exit code 0, runs/mc2b_logs/queue_mc2b.log).
