# Long tolerance test at 60 MHz, 9 mT, small particle (reading fixed in PROTOKOLL 7.100)

| seed | atol | class | c_settle | P last 10 cycles (W/cm³, mean ± SD) | C′ value (W/cm³) | start: E/K_d and w_loop of cycle 0, relative to L5 |
|---|---|---|---|---|---|---|
| 0 | 1e-6 | stays irregular | – | 35.7 ± 1.8 | 35.3 | +0.0e+00 / -8.1e-07 |
| 0 | 1e-5 | stays irregular | – | 35.7 ± 1.8 | 35.3 | reference |
| 1 | 1e-6 | stays irregular | – | 28.5 ± 1.1 | 28.9 | +0.0e+00 / -2.5e-07 |
| 1 | 1e-5 | stays irregular | – | 28.5 ± 1.1 | 28.9 | reference |
| 2 | 1e-6 | marginal | – | 38.1 ± 4.4 | 36.4 | +0.0e+00 / -8.5e-06 |
| 2 | 1e-5 | marginal | – | 38.1 ± 4.4 | 36.4 | reference |

**Reading (motion): stays irregular**

**Reading (loss value): the C′ value holds**

## Exploratory (not pre-registered; added 2026-10-15 after the readings)

| seed | largest abs. w_loop(L6) / w_loop(L5) − 1 in cycles 0–9 / 10–19 / 20–29 / 30–39 | largest abs. closure ratio − 1 in the same blocks | running mean of P within ±5 % from cycle (L6 / L5) |
|---|---|---|---|
| 0 | 6.0e-06 / 5.6e-06 / 4.6e-05 / 1.2e-04 | 6.9e-04 / 2.9e-03 / 5.6e-02 / 8.3e-03 | 2 / 2 |
| 1 | 5.4e-06 / 4.2e-06 / 4.6e-06 / 7.3e-06 | 1.2e-03 / 4.3e-04 / 1.2e-03 / 6.3e-04 | 10 / 10 |
| 2 | 8.5e-06 / 4.5e-06 / 7.4e-06 / 5.8e-06 | 9.2e-04 / 8.0e-03 / 2.6e-02 / 9.0e-03 | 6 / 6 |

| seed | run | cycles kept | R1 | R2 | R3 | R4 … R7 | rest | loss share, texture top 10 % | loss share, r ≥ 0.8 |
|---|---|---|---|---|---|---|---|---|---|
| 0 | L6 | 36 | 0.915 | 0.052 | 0.027 | 0.062 | -0.056 | 0.558 | 0.314 |
| 0 | L5 | 37 | 0.915 | 0.053 | 0.028 | 0.064 | -0.059 | 0.557 | 0.313 |
| 0 | C′ (window) | 4 | 0.882 | 0.083 | 0.048 | 0.114 | -0.127 | 0.558 | 0.320 |
| 1 | L6 | 37 | 0.920 | 0.033 | 0.019 | 0.040 | -0.013 | 0.603 | 0.280 |
| 1 | L5 | 37 | 0.919 | 0.033 | 0.019 | 0.040 | -0.011 | 0.602 | 0.280 |
| 1 | C′ (window) | 4 | 0.854 | 0.067 | 0.047 | 0.135 | -0.104 | 0.610 | 0.296 |
| 2 | L6 | 34 | 0.923 | 0.048 | 0.012 | 0.038 | -0.020 | 0.520 | 0.285 |
| 2 | L5 | 34 | 0.922 | 0.048 | 0.012 | 0.038 | -0.020 | 0.520 | 0.285 |
| 2 | C′ (window) | 4 | 0.860 | 0.069 | 0.029 | 0.114 | -0.072 | 0.531 | 0.314 |

rest = 1 − Σ R_k; a negative rest is the per-cell over-completeness of the harmonic sum (7.90).
