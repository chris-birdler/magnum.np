# Sobol design v2: analysis (PLAN 4.0)

Runs in the fits: 116 (flagged 0); bootstrap 1000 resamples (0 failed). Range effect = output at factor max minus at min, other factors averaged over the design points. φ: the effect on ln E[w] contains the dilution 0.227 (value without it in brackets).

**Declared deviation (PROTOKOLL 7.54): steady-tail rule.** Only the measured cycles after the last cycle with |w_dis / w_loop − 1| > 0.20 are used. Dropped: 227 of 4060 measured cycles; amplitudes with fewer than 2 kept cycles (missing): 0 in 0 runs.

## Primary outputs (8 tests, Holm 5 %)

| output | factor | main fit (95 % CI) | 90 % CI | class | fit A | fit B |
|---|---|---|---|---|---|---|
| beta_all | lnQ | +0.150 (+0.067 … +0.222) | +0.080 … +0.209 | lever | +0.131 ± 0.039 | +0.126 ± 0.040 |
| beta_all | r_p | +0.145 (+0.016 … +0.257) | +0.043 … +0.237 | not determined | +0.131 ± 0.058 | +0.131 ± 0.053 |
| beta_all | phi | -0.064 (-0.163 … +0.025) | -0.150 … +0.008 | not determined | -0.059 ± 0.049 | -0.054 ± 0.050 |
| beta_all | lnd | -0.572 (-0.668 … -0.475) | -0.656 … -0.493 | lever | -0.582 ± 0.045 | -0.594 ± 0.048 |
| lnw@50mT | lnQ | -0.065 (-0.254 … +0.117) | -0.224 … +0.087 | not determined | -0.103 ± 0.102 | -0.119 ± 0.074 |
| lnw@50mT | r_p | +0.036 (-0.202 … +0.263) | -0.158 … +0.233 | not determined | +0.027 ± 0.126 | +0.063 ± 0.097 |
| lnw@50mT | phi | +0.064 [-0.163] (-0.144 … +0.286) | -0.112 … +0.236 | not determined | +0.029 ± 0.123 | -0.007 ± 0.093 |
| lnw@50mT | lnd | +0.308 (+0.091 … +0.508) | +0.140 … +0.471 | lever | +0.323 ± 0.111 | +0.311 ± 0.088 |

## Secondary outputs (hints only)

| output | factor | main fit (95 % CI) | fit A | fit B | without flagged runs | with w_dis |
|---|---|---|---|---|---|---|
| beta@9-35mT | lnQ | +0.136 (-0.012 … +0.280) | +0.128 ± 0.075 | +0.130 ± 0.061 | +0.137 | +0.132 |
| beta@9-35mT | r_p | +0.146 (-0.067 … +0.363) | +0.117 ± 0.114 | +0.045 ± 0.080 | +0.137 | +0.127 |
| beta@9-35mT | phi | -0.104 (-0.275 … +0.072) | -0.075 ± 0.088 | -0.079 ± 0.077 | -0.109 | -0.108 |
| beta@9-35mT | lnd | -0.308 (-0.503 … -0.119) | -0.280 ± 0.094 | -0.297 ± 0.073 | -0.314 | -0.290 |
| beta@35-70mT | lnQ | +0.115 (-0.143 … +0.330) | +0.065 ± 0.115 | +0.103 ± 0.091 | +0.089 | +0.085 |
| beta@35-70mT | r_p | +0.291 (+0.017 … +0.547) | +0.273 ± 0.129 | +0.291 ± 0.120 | +0.261 | +0.252 |
| beta@35-70mT | phi | -0.070 (-0.371 … +0.178) | -0.113 ± 0.134 | -0.105 ± 0.115 | -0.078 | -0.076 |
| beta@35-70mT | lnd | -0.998 (-1.251 … -0.720) | -1.033 ± 0.123 | -1.038 ± 0.109 | -1.031 | -1.040 |
| beta@70-150mT | lnQ | +0.305 (+0.103 … +0.516) | +0.308 ± 0.121 | +0.299 ± 0.084 | +0.267 | +0.261 |
| beta@70-150mT | r_p | +0.020 (-0.236 … +0.254) | +0.057 ± 0.144 | +0.099 ± 0.110 | +0.027 | +0.030 |
| beta@70-150mT | phi | -0.006 (-0.266 … +0.222) | +0.000 ± 0.154 | +0.047 ± 0.106 | -0.021 | -0.024 |
| beta@70-150mT | lnd | -0.351 (-0.592 … -0.149) | -0.400 ± 0.134 | -0.362 ± 0.101 | -0.383 | -0.381 |
| lnw@10mT | lnQ | -0.309 (-0.476 … -0.144) | -0.302 ± 0.089 | -0.274 ± 0.084 | -0.296 | -0.283 |
| lnw@10mT | r_p | -0.228 (-0.445 … +0.012) | -0.190 ± 0.126 | -0.206 ± 0.110 | -0.195 | -0.178 |
| lnw@10mT | phi | +0.294 (+0.099 … +0.521) | +0.257 ± 0.118 | +0.266 ± 0.106 | +0.290 | +0.290 |
| lnw@10mT | lnd | +1.107 (+0.920 … +1.307) | +1.122 ± 0.106 | +1.171 ± 0.101 | +1.142 | +1.111 |
| lnw@100mT | lnQ | -0.017 (-0.159 … +0.091) | -0.046 ± 0.069 | -0.052 ± 0.059 | -0.029 | -0.028 |
| lnw@100mT | r_p | +0.101 (-0.071 … +0.247) | +0.090 ± 0.087 | +0.085 ± 0.077 | +0.096 | +0.103 |
| lnw@100mT | phi | +0.160 (+0.019 … +0.328) | +0.145 ± 0.088 | +0.181 ± 0.074 | +0.139 | +0.144 |
| lnw@100mT | lnd | -0.190 (-0.339 … -0.016) | -0.180 ± 0.085 | -0.150 ± 0.071 | -0.188 | -0.187 |

## Scatter model (fit A residuals, Holm over 4 slopes)

- beta_all: slopes lnQ +0.16, r_p +0.02, phi +0.35, lnd +0.24; p lnQ 0.587, r_p 0.948, phi 0.343, lnd 0.482; significant: none
- lnw@50mT: slopes lnQ +0.82, r_p +0.02, phi +0.62, lnd -0.17; p lnQ 0.008, r_p 0.962, phi 0.106, lnd 0.636; significant: lnQ

## Mesh-check trigger

- hits: beta@70-150mT × lnQ, lnw@10mT × lnQ, beta@35-70mT × r_p, lnw@10mT × phi, lnw@100mT × phi, beta@9-35mT × lnd, beta@35-70mT × lnd, beta@70-150mT × lnd, beta_all × lnd, lnw@10mT × lnd, lnw@50mT × lnd, lnw@100mT × lnd
- corners to check: K2, K3

## D3 check (fit A residuals, d/L_eff < 10 vs ≥ 10)

- beta_all: n(<10) = 35, Δ = +0.002, p = 0.951
- lnw@50mT: n(<10) = 35, Δ = +0.004, p = 0.947

## α 0.01 centre test (pairs: 8)

- ln(w α0.01 / w α0.1) per amplitude: 9 mT +nan ± nan, 20 mT +nan ± nan, 35 mT +nan ± nan, 50 mT +nan ± nan, 70 mT +nan ± nan, 100 mT +nan ± nan, 150 mT +nan ± nan
- lnR@9mT: R = nan (95 % CI nan … nan, n = 8): mixed
- lnR@50mT: R = nan (95 % CI nan … nan, n = 8): mixed
