# Loop closure and settling of the main run (PROTOKOLL 7.101, eval_closure_state.py)

## 1. Protocol (116 runs)

| protocol | cycles_per_amp | max_cycles_per_amp | min_calm | runs |
|---|---|---|---|---|
| virgin_asc | 7 | 7 | 0 | 103 |
| virgin_asc | 7 | 15 | 2 | 13 |

Cycles per stage in the runs with min_calm > 0: 7 … 12.

## 2. Steady tails (812 stages)

| B (mT) | stages with tail ≥ 2 | smallest tail | all measured cycles calm |
|---|---|---|---|
| 9 | 116 / 116 | 2 | 70 |
| 20 | 116 / 116 | 2 | 98 |
| 35 | 116 / 116 | 2 | 91 |
| 50 | 116 / 116 | 2 | 102 |
| 70 | 116 / 116 | 2 | 112 |
| 100 | 116 / 116 | 2 | 109 |
| 150 | 116 / 116 | 4 | 114 |

## 3. Closure (mean of the steady-tail cycles per stage; stages above a limit)

| B (mT) | > 2 % | > 0.5 % | > 0.2 % | > 0.1 % | median | median of closure / half width | cycles > 25 % of the half width | > 50 % |
|---|---|---|---|---|---|---|---|---|
| 9 | 0 | 1 | 9 | 23 | 0.02 % | 0.014 | 8.6 % | 3.1 % |
| 20 | 0 | 1 | 5 | 20 | 0.03 % | 0.027 | 11.9 % | 2.6 % |
| 35 | 0 | 0 | 13 | 35 | 0.05 % | 0.045 | 15.0 % | 2.1 % |
| 50 | 0 | 0 | 7 | 35 | 0.06 % | 0.077 | 11.4 % | 2.9 % |
| 70 | 0 | 1 | 12 | 35 | 0.05 % | 0.060 | 13.9 % | 3.0 % |
| 100 | 0 | 1 | 12 | 40 | 0.08 % | 0.094 | 15.5 % | 2.8 % |
| 150 | 0 | 0 | 12 | 47 | 0.08 % | 0.100 | 18.0 % | 2.9 % |
| sum | 0 | 4 | 70 | 235 | | | | |

## 4. Monotonic change or recurrent jumps (criteria fixed before the data were looked at)

| B (mT) | Δe / w_loop mean ± SE | same-sign share (pairs) | binomial SE | stages with an event | events release / uptake | reading |
|---|---|---|---|---|---|---|
| 9 | -0.74 ± 0.25 % | 38.9 % (398) | 2.5 % | 46 | 62 / 2 | recurrent jumps |
| 20 | -0.68 ± 0.17 % | 34.8 % (431) | 2.4 % | 18 | 29 / 0 | recurrent jumps |
| 35 | -0.18 ± 0.17 % | 36.4 % (409) | 2.5 % | 25 | 48 / 0 | recurrent jumps |
| 50 | -0.27 ± 0.09 % | 24.6 % (439) | 2.4 % | 14 | 20 / 0 | recurrent jumps |
| 70 | -0.31 ± 0.11 % | 29.6 % (453) | 2.3 % | 4 | 6 / 0 | recurrent jumps |
| 100 | -0.23 ± 0.13 % | 29.1 % (450) | 2.4 % | 7 | 9 / 0 | recurrent jumps |
| 150 | +0.04 ± 0.06 % | 24.7 % (462) | 2.3 % | 2 | 2 / 0 | recurrent jumps |

## 5. Same-sign share against 1/3 (fixed state) and 1/2 (wandering state); SE: bootstrap over whole stages, 2000 draws (added after the data were looked at)

| group | share | SE bootstrap | SE binomial at the share | pairs | distance to 1/2 (SE) | distance to 1/3 (SE) |
|---|---|---|---|---|---|---|
| 9 mT | 38.9 % | 2.7 % | 2.4 % | 398 | -4.1 | +2.1 |
| 20 mT | 34.8 % | 2.5 % | 2.3 % | 431 | -6.0 | +0.6 |
| 35 mT | 36.4 % | 2.8 % | 2.4 % | 409 | -4.9 | +1.1 |
| 50 mT | 24.6 % | 2.1 % | 2.1 % | 439 | -12.0 | -4.1 |
| 70 mT | 29.6 % | 2.3 % | 2.1 % | 453 | -8.9 | -1.6 |
| 100 mT | 29.1 % | 2.4 % | 2.1 % | 450 | -8.7 | -1.8 |
| 150 mT | 24.7 % | 2.3 % | 2.0 % | 462 | -11.2 | -3.8 |
| all B | 30.9 % | 0.9 % | 0.8 % | 3042 | -20.3 | -2.6 |
| C′ 15 MHz large | 27.5 % | 6.9 % | 5.0 % | 80 | -3.3 | -0.8 |
| C′ 15 MHz small | 23.8 % | 5.1 % | 4.8 % | 80 | -5.2 | -1.9 |
| C′ 30 MHz large | 28.0 % | 5.5 % | 5.2 % | 75 | -4.0 | -1.0 |
| C′ 30 MHz small | 33.3 % | 3.8 % | 5.3 % | 78 | -4.4 | +0.0 |
| C′ 60 MHz large | 33.8 % | 4.1 % | 5.3 % | 80 | -3.9 | +0.1 |
| C′ 60 MHz small | 26.4 % | 4.1 % | 5.2 % | 72 | -5.7 | -1.7 |

## 6. C′ per frequency (corner K1, 3 seeds × 7 amplitudes)

| f | particle | stages with tail ≥ 2 | largest number of cycles | at the limit without tail | events release / uptake | median closure | median cycle SD of w |
|---|---|---|---|---|---|---|---|
| 15 MHz | large | 21 / 21 | 7 | 0 | 4 / 0 | 0.006 % | 1.2 % |
| 15 MHz | small | 21 / 21 | 7 | 0 | 3 / 0 | 0.004 % | 1.6 % |
| 30 MHz | large | 21 / 21 | 9 | 0 | 9 / 0 | 0.048 % | 2.8 % |
| 30 MHz | small | 21 / 21 | 8 | 0 | 4 / 0 | 0.053 % | 2.3 % |
| 60 MHz | large | 21 / 21 | 7 | 0 | 4 / 0 | 0.068 % | 2.4 % |
| 60 MHz | small | 21 / 21 | 12 | 0 | 18 / 0 | 0.101 % | 4.4 % |

## 7. K1_k at 9 mT (main run, atol 1e-5): closure and w_dis / w_loop per cycle

- K1_0: closure 0.00862 / 0.00482 / 0.00035 / 0.00026 / 0.00003 / 0.00000 / 0.00000; w_dis / w_loop 6.134 / 1.473 / 0.963 / 1.006 / 1.018 / 1.015 / 1.015; c (eval_atol30) 0.055 %
- K1_1: closure 0.00140 / 0.00110 / 0.00018 / 0.00003 / 0.00000 / 0.00000 / 0.00000; w_dis / w_loop 0.836 / 0.923 / 1.041 / 1.034 / 1.030 / 1.031 / 1.030; c (eval_atol30) 0.057 %
- K1_2: closure 0.00201 / 0.00122 / 0.00000 / 0.00004 / 0.00000 / 0.00000 / 0.00000; w_dis / w_loop 0.756 / 0.966 / 1.021 / 1.016 / 1.017 / 1.016 / 1.017; c (eval_atol30) 0.054 %

## 8. B0 (atol 1e-5) against B0b (atol 1e-6): common complete cycles

- s941, cycles [0, 1, 2]: closure 1e-6 / 1e-5 0.9999 / 1.0002 / 1.0003; W_loop 1e-6 / 1e-5 1.000006 / 1.000038 / 0.999982; median ratio 1.0002
- s942, cycles [0, 1]: closure 1e-6 / 1e-5 1.0014 / 0.9999; W_loop 1e-6 / 1e-5 0.999988 / 0.999999; median ratio 1.0012
- s943, cycles [0, 1, 2]: closure 1e-6 / 1e-5 1.0013 / 1.0000 / 1.0000; W_loop 1e-6 / 1e-5 0.999952 / 0.999997 / 0.999999; median ratio 1.0013
- s944, cycles [0, 1, 2, 3]: closure 1e-6 / 1e-5 1.0001 / 1.0001 / 0.9998 / 1.0011; W_loop 1e-6 / 1e-5 1.000022 / 0.999995 / 1.000002 / 1.000002; median ratio 1.0001

**Reading: physical in 4 of 4 seeds, numerical in 0.**
