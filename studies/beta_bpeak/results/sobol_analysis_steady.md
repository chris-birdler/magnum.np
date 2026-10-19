# Sobol design v2: analysis (PLAN 4.0)

Runs in the fits: 116 (flagged 0); bootstrap 1000 resamples (0 failed). Range effect = output at factor max minus at min, other factors averaged over the design points. φ: the effect on ln E[w] contains the dilution 0.227 (value without it in brackets).

**Declared deviation (PROTOKOLL 7.54): steady-tail rule.** Only the measured cycles after the last cycle with |w_dis / w_loop − 1| > 0.20 are used. Dropped: 223 of 4060 measured cycles; amplitudes without a kept cycle (missing): 4 in 4 runs.

## Primary outputs (8 tests, Holm 5 %)

| output | factor | main fit (95 % CI) | 90 % CI | class | fit A | fit B |
|---|---|---|---|---|---|---|
| beta_all | lnQ | +0.151 (+0.068 … +0.224) | +0.082 … +0.211 | lever | +0.132 ± 0.038 | +0.127 ± 0.039 |
| beta_all | r_p | +0.145 (+0.016 … +0.256) | +0.043 … +0.238 | not determined | +0.130 ± 0.058 | +0.133 ± 0.051 |
| beta_all | phi | -0.065 (-0.165 … +0.023) | -0.152 … +0.007 | not determined | -0.062 ± 0.050 | -0.056 ± 0.049 |
| beta_all | lnd | -0.572 (-0.669 … -0.476) | -0.656 … -0.492 | lever | -0.581 ± 0.045 | -0.592 ± 0.047 |
| lnw@50mT | lnQ | -0.066 (-0.254 … +0.117) | -0.225 … +0.088 | not determined | -0.104 ± 0.102 | -0.120 ± 0.073 |
| lnw@50mT | r_p | +0.037 (-0.200 … +0.264) | -0.156 … +0.235 | not determined | +0.027 ± 0.126 | +0.063 ± 0.096 |
| lnw@50mT | phi | +0.062 [-0.165] (-0.144 … +0.284) | -0.114 … +0.234 | not determined | +0.030 ± 0.123 | -0.006 ± 0.092 |
| lnw@50mT | lnd | +0.308 (+0.093 … +0.506) | +0.140 … +0.470 | lever | +0.324 ± 0.111 | +0.312 ± 0.088 |

## Secondary outputs (hints only)

| output | factor | main fit (95 % CI) | fit A | fit B | without flagged runs | with w_dis |
|---|---|---|---|---|---|---|
| beta@9-35mT | lnQ | +0.125 (-0.027 … +0.273) | +0.120 ± 0.076 | +0.122 ± 0.062 | +0.128 | +0.124 |
| beta@9-35mT | r_p | +0.137 (-0.086 … +0.359) | +0.108 ± 0.114 | +0.032 ± 0.081 | +0.123 | +0.111 |
| beta@9-35mT | phi | -0.108 (-0.278 … +0.074) | -0.071 ± 0.089 | -0.074 ± 0.078 | -0.110 | -0.108 |
| beta@9-35mT | lnd | -0.298 (-0.498 … -0.108) | -0.274 ± 0.094 | -0.287 ± 0.074 | -0.304 | -0.281 |
| beta@35-70mT | lnQ | +0.124 (-0.130 … +0.345) | +0.067 ± 0.115 | +0.105 ± 0.091 | +0.095 | +0.096 |
| beta@35-70mT | r_p | +0.304 (+0.027 … +0.560) | +0.282 ± 0.128 | +0.306 ± 0.120 | +0.281 | +0.278 |
| beta@35-70mT | phi | -0.069 (-0.376 … +0.191) | -0.132 ± 0.134 | -0.126 ± 0.115 | -0.087 | -0.082 |
| beta@35-70mT | lnd | -1.015 (-1.269 … -0.747) | -1.044 ± 0.123 | -1.052 ± 0.109 | -1.050 | -1.066 |
| beta@70-150mT | lnQ | +0.314 (+0.111 … +0.531) | +0.318 ± 0.122 | +0.315 ± 0.081 | +0.277 | +0.271 |
| beta@70-150mT | r_p | +0.018 (-0.239 … +0.257) | +0.061 ± 0.145 | +0.112 ± 0.106 | +0.031 | +0.034 |
| beta@70-150mT | phi | +0.008 (-0.257 … +0.235) | +0.013 ± 0.155 | +0.075 ± 0.102 | -0.008 | -0.013 |
| beta@70-150mT | lnd | -0.349 (-0.591 … -0.144) | -0.400 ± 0.136 | -0.351 ± 0.097 | -0.383 | -0.380 |
| lnw@10mT | lnQ | -0.309 (-0.477 … -0.148) | -0.311 ± 0.091 | -0.281 ± 0.084 | -0.296 | -0.288 |
| lnw@10mT | r_p | -0.226 (-0.445 … +0.015) | -0.191 ± 0.128 | -0.208 ± 0.110 | -0.194 | -0.178 |
| lnw@10mT | phi | +0.293 (+0.099 … +0.522) | +0.255 ± 0.121 | +0.260 ± 0.107 | +0.291 | +0.289 |
| lnw@10mT | lnd | +1.103 (+0.912 … +1.303) | +1.128 ± 0.109 | +1.191 ± 0.101 | +1.140 | +1.114 |
| lnw@100mT | lnQ | -0.016 (-0.158 … +0.091) | -0.045 ± 0.069 | -0.050 ± 0.058 | -0.028 | -0.029 |
| lnw@100mT | r_p | +0.102 (-0.072 … +0.248) | +0.091 ± 0.087 | +0.086 ± 0.076 | +0.096 | +0.103 |
| lnw@100mT | phi | +0.158 (+0.019 … +0.326) | +0.146 ± 0.088 | +0.182 ± 0.073 | +0.140 | +0.144 |
| lnw@100mT | lnd | -0.191 (-0.339 … -0.016) | -0.179 ± 0.085 | -0.149 ± 0.070 | -0.188 | -0.188 |

## Scatter model (fit A residuals, Holm over 4 slopes)

- beta_all: slopes lnQ +0.17, r_p +0.06, phi +0.28, lnd +0.28; p lnQ 0.547, r_p 0.877, phi 0.443, lnd 0.407; significant: none
- lnw@50mT: slopes lnQ +0.82, r_p -0.01, phi +0.62, lnd -0.14; p lnQ 0.008, r_p 0.986, phi 0.108, lnd 0.699; significant: lnQ

## Mesh-check trigger

- hits: beta@70-150mT × lnQ, beta_all × lnQ, lnw@10mT × lnQ, beta@35-70mT × r_p, lnw@10mT × phi, lnw@100mT × phi, beta@9-35mT × lnd, beta@35-70mT × lnd, beta@70-150mT × lnd, beta_all × lnd, lnw@10mT × lnd, lnw@50mT × lnd, lnw@100mT × lnd
- corners to check: K2, K3

## D3 check (fit A residuals, d/L_eff < 10 vs ≥ 10)

- beta_all: n(<10) = 35, Δ = +0.001, p = 0.962
- lnw@50mT: n(<10) = 35, Δ = +0.004, p = 0.947

## α 0.01 centre test (pairs: 8)

- ln(w α0.01 / w α0.1) per amplitude: 9 mT +nan ± nan, 20 mT +nan ± nan, 35 mT +nan ± nan, 50 mT +nan ± nan, 70 mT +nan ± nan, 100 mT +nan ± nan, 150 mT +nan ± nan
- lnR@9mT: R = nan (95 % CI nan … nan, n = 8): mixed
- lnR@50mT: R = nan (95 % CI nan … nan, n = 8): mixed
