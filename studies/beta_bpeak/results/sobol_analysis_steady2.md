# Sobol design v2: analysis (PLAN 4.0)

Runs in the fits: 116 (flagged 0); bootstrap 1000 resamples (0 failed). Range effect = output at factor max minus at min, other factors averaged over the design points. φ: the effect on ln E[w] contains the dilution 0.227 (value without it in brackets).

**Declared deviation (PROTOKOLL 7.54): steady-tail rule.** Only the measured cycles after the last cycle with |w_dis / w_loop − 1| > 0.20 are used. Dropped: 232 of 4060 measured cycles; amplitudes with fewer than 2 kept cycles (missing): 13 in 13 runs.

## Primary outputs (8 tests, Holm 5 %)

| output | factor | main fit (95 % CI) | 90 % CI | class | fit A | fit B |
|---|---|---|---|---|---|---|
| beta_all | lnQ | +0.152 (+0.070 … +0.228) | +0.084 … +0.211 | lever | +0.137 ± 0.038 | +0.132 ± 0.039 |
| beta_all | r_p | +0.157 (+0.033 … +0.268) | +0.054 … +0.248 | lever | +0.140 ± 0.057 | +0.141 ± 0.051 |
| beta_all | phi | -0.058 (-0.156 … +0.032) | -0.144 … +0.019 | not a lever | -0.061 ± 0.049 | -0.056 ± 0.049 |
| beta_all | lnd | -0.577 (-0.677 … -0.481) | -0.663 … -0.497 | lever | -0.588 ± 0.046 | -0.598 ± 0.047 |
| lnw@50mT | lnQ | -0.066 (-0.254 … +0.117) | -0.225 … +0.088 | not determined | -0.104 ± 0.102 | -0.120 ± 0.073 |
| lnw@50mT | r_p | +0.037 (-0.200 … +0.264) | -0.156 … +0.235 | not determined | +0.027 ± 0.126 | +0.063 ± 0.096 |
| lnw@50mT | phi | +0.062 [-0.165] (-0.144 … +0.284) | -0.114 … +0.234 | not determined | +0.030 ± 0.123 | -0.006 ± 0.092 |
| lnw@50mT | lnd | +0.308 (+0.093 … +0.506) | +0.140 … +0.470 | lever | +0.324 ± 0.111 | +0.312 ± 0.088 |

## Secondary outputs (hints only)

| output | factor | main fit (95 % CI) | fit A | fit B | without flagged runs | with w_dis |
|---|---|---|---|---|---|---|
| beta@9-35mT | lnQ | +0.114 (-0.036 … +0.264) | +0.132 ± 0.078 | +0.133 ± 0.065 | +0.121 | +0.114 |
| beta@9-35mT | r_p | +0.167 (-0.054 … +0.392) | +0.134 ± 0.115 | +0.055 ± 0.085 | +0.148 | +0.136 |
| beta@9-35mT | phi | -0.098 (-0.275 … +0.101) | -0.061 ± 0.090 | -0.061 ± 0.082 | -0.100 | -0.103 |
| beta@9-35mT | lnd | -0.307 (-0.513 … -0.121) | -0.292 ± 0.098 | -0.306 ± 0.078 | -0.316 | -0.290 |
| beta@35-70mT | lnQ | +0.159 (-0.095 … +0.386) | +0.001 ± 0.130 | +0.078 ± 0.095 | +0.124 | +0.119 |
| beta@35-70mT | r_p | +0.307 (+0.032 … +0.557) | +0.264 ± 0.134 | +0.307 ± 0.125 | +0.285 | +0.283 |
| beta@35-70mT | phi | -0.049 (-0.348 … +0.223) | -0.203 ± 0.159 | -0.139 ± 0.120 | -0.066 | -0.058 |
| beta@35-70mT | lnd | -1.035 (-1.320 … -0.752) | -1.056 ± 0.131 | -1.044 ± 0.114 | -1.061 | -1.070 |
| beta@70-150mT | lnQ | +0.307 (+0.102 … +0.528) | +0.317 ± 0.122 | +0.313 ± 0.081 | +0.270 | +0.264 |
| beta@70-150mT | r_p | +0.016 (-0.240 … +0.254) | +0.060 ± 0.145 | +0.112 ± 0.106 | +0.028 | +0.031 |
| beta@70-150mT | phi | -0.001 (-0.263 … +0.234) | +0.012 ± 0.155 | +0.073 ± 0.102 | -0.017 | -0.021 |
| beta@70-150mT | lnd | -0.350 (-0.593 … -0.148) | -0.401 ± 0.136 | -0.352 ± 0.097 | -0.385 | -0.383 |
| lnw@10mT | lnQ | -0.309 (-0.475 … -0.149) | -0.307 ± 0.097 | -0.285 ± 0.090 | -0.298 | -0.285 |
| lnw@10mT | r_p | -0.261 (-0.475 … -0.037) | -0.242 ± 0.120 | -0.239 ± 0.116 | -0.224 | -0.207 |
| lnw@10mT | phi | +0.275 (+0.076 … +0.506) | +0.220 ± 0.121 | +0.239 ± 0.112 | +0.274 | +0.275 |
| lnw@10mT | lnd | +1.123 (+0.943 … +1.332) | +1.150 ± 0.113 | +1.196 ± 0.108 | +1.157 | +1.126 |
| lnw@100mT | lnQ | -0.013 (-0.157 … +0.096) | -0.042 ± 0.069 | -0.047 ± 0.059 | -0.025 | -0.023 |
| lnw@100mT | r_p | +0.102 (-0.072 … +0.250) | +0.091 ± 0.087 | +0.085 ± 0.077 | +0.096 | +0.103 |
| lnw@100mT | phi | +0.159 (+0.019 … +0.329) | +0.146 ± 0.088 | +0.182 ± 0.074 | +0.140 | +0.144 |
| lnw@100mT | lnd | -0.185 (-0.336 … -0.008) | -0.176 ± 0.085 | -0.146 ± 0.071 | -0.184 | -0.179 |

## Scatter model (fit A residuals, Holm over 4 slopes)

- beta_all: slopes lnQ +0.08, r_p -0.05, phi +0.44, lnd +0.40; p lnQ 0.795, r_p 0.900, phi 0.236, lnd 0.249; significant: none
- lnw@50mT: slopes lnQ +0.82, r_p -0.01, phi +0.62, lnd -0.14; p lnQ 0.008, r_p 0.987, phi 0.108, lnd 0.699; significant: lnQ

## Mesh-check trigger

- hits: beta@70-150mT × lnQ, beta_all × lnQ, lnw@10mT × lnQ, beta@35-70mT × r_p, beta_all × r_p, lnw@10mT × r_p, lnw@10mT × phi, lnw@100mT × phi, beta@9-35mT × lnd, beta@35-70mT × lnd, beta@70-150mT × lnd, beta_all × lnd, lnw@10mT × lnd, lnw@50mT × lnd, lnw@100mT × lnd
- corners to check: K2, K3

## D3 check (fit A residuals, d/L_eff < 10 vs ≥ 10)

- beta_all: n(<10) = 35, Δ = +0.004, p = 0.874
- lnw@50mT: n(<10) = 35, Δ = +0.004, p = 0.947

## α 0.01 centre test (pairs: 8)

- ln(w α0.01 / w α0.1) per amplitude: 9 mT +nan ± nan, 20 mT +nan ± nan, 35 mT +nan ± nan, 50 mT +nan ± nan, 70 mT +nan ± nan, 100 mT +nan ± nan, 150 mT +nan ± nan
- lnR@9mT: R = nan (95 % CI nan … nan, n = 8): mixed
- lnR@50mT: R = nan (95 % CI nan … nan, n = 8): mixed
