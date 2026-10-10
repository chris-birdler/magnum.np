# Stage 1: jumps in the dissipation time series (exploratory)

n = 115 α 0.1 runs (Sobol, centre, corners), last 3 cycles per amplitude, 256 samples per cycle. Medians over the runs, 95 % CI of the median by bootstrap over the runs (2000 resamples).

## A. Measures against the amplitude

| B (mT) | RF non-smooth | RN non-repeating | RH sharp repeating | dS5 burst excess | rho event coupling (runs > 0) | samples with an n60 change | R_np snapshots | balance off (runs) |
|---|---|---|---|---|---|---|---|---|
| 9 | 0.041 (0.027 … 0.057) | 0.006 (0.005 … 0.008) | 0.031 (0.023 … 0.045) | 0.005 | 0.01 (0.00 … 0.02), 69/108 | 17 % | 0.000 | 3 |
| 20 | 0.139 (0.117 … 0.166) | 0.022 (0.009 … 0.051) | 0.119 (0.107 … 0.133) | 0.028 | 0.05 (0.02 … 0.07), 81/106 | 33 % | – | 5 |
| 35 | 0.230 (0.206 … 0.247) | 0.109 (0.056 … 0.153) | 0.172 (0.165 … 0.190) | 0.046 | 0.10 (0.08 … 0.11), 91/114 | 49 % | – | 8 |
| 50 | 0.231 (0.211 … 0.242) | 0.122 (0.103 … 0.144) | 0.181 (0.166 … 0.193) | 0.045 | 0.13 (0.11 … 0.14), 100/112 | 60 % | 0.029 | 0 |
| 70 | 0.217 (0.206 … 0.226) | 0.117 (0.103 … 0.127) | 0.165 (0.156 … 0.178) | 0.036 | 0.13 (0.11 … 0.14), 104/115 | 69 % | – | 2 |
| 100 | 0.185 (0.176 … 0.191) | 0.119 (0.110 … 0.132) | 0.139 (0.131 … 0.144) | 0.027 | 0.11 (0.08 … 0.11), 108/115 | 76 % | – | 2 |
| 150 | 0.153 (0.146 … 0.162) | 0.101 (0.093 … 0.111) | 0.109 (0.103 … 0.113) | 0.018 | 0.07 (0.06 … 0.08), 109/115 | 81 % | 0.025 | 0 |

balance off: |w_dis / w_loop − 1| > 0.20 over the kept cycles (state still changes, or a jump falls between the samples). Check without these run-amplitude pairs:

| B (mT) | RF median, all | RF median, balance ok |
|---|---|---|
| 9 | 0.041 (n 115) | 0.040 (n 112) |
| 20 | 0.139 (n 115) | 0.134 (n 110) |
| 35 | 0.230 (n 115) | 0.220 (n 107) |
| 50 | 0.231 (n 115) | 0.231 (n 115) |
| 70 | 0.217 (n 115) | 0.217 (n 113) |
| 100 | 0.185 (n 115) | 0.185 (n 113) |
| 150 | 0.153 (n 115) | 0.153 (n 115) |

## B. Amplitude B_pk of the largest RF per run

All runs: median B_pk 42 mT (95 % CI 36 … 51), n 115. B_pk at 9 or 150 mT (end of the range): 8 runs.

Linear model ln B_pk = c0 + Σ c_j x_j (x_j coded −1 … +1, OLS, HC3 standard errors). Range effect = factor on B_pk from the low end to the high end of the factor, with 95 % CI:

| factor (low → high end) | factor on B_pk (95 % CI) | p |
|---|---|---|
| L_eff 30 → 12 l_ex (K_eff up) | 1.41 (1.07 … 1.87) | 0.015 |
| r_p 0 → 3 | 1.32 (0.91 … 1.93) | 0.147 |
| φ 0.55 → 0.69 | 0.99 (0.72 … 1.36) | 0.957 |
| d/l_ex 212 → 300 (A 20 → 10 pJ/m) | 0.46 (0.33 … 0.62) | 0.000 |

R² = 0.21. Groups (median B_pk, 95 % CI):

| group | n | median B_pk (mT) |
|---|---|---|
| A low third (9.7 … 12.5 pJ/m) | 39 | 32 (28 … 39) |
| A mid third | 38 | 36 (32 … 46) |
| A high third (15.7 … 20.0 pJ/m) | 38 | 56 (51 … 65) |
| L_eff 12 l_ex | 26 | 54 (35 … 70) |
| L_eff 18 l_ex | 38 | 45 (35 … 53) |
| L_eff 24 l_ex | 24 | 36 (29 … 49) |
| L_eff 30 l_ex | 27 | 41 (32 … 51) |

## C. Mesh check: the same measures at dx 3 and dx 1.5 l_ex (paired seeds, one job per amplitude)

Jobs with 5 cycles, the last 3 kept. Medians over the seeds. n60 max = largest number of pairs > 60° in the kept cycles. Sign test (two-sided, exact) over the seeds of both points.

| B (mT) | point | n pairs | RF dx 3 | RF dx 1.5 | seeds with RF(dx 3) > RF(dx 1.5) | RN dx 3 | RN dx 1.5 | n60 max dx 3 | n60 max dx 1.5 | balance off dx 3 / dx 1.5 |
|---|---|---|---|---|---|---|---|---|---|---|
| 9 | K3 (A 20 pJ/m) | 6 | 0.037 | 0.110 | 0/6 | 0.036 | 0.040 | 23 | 11 | 0 / 0 |
| 9 | base (A 10 pJ/m) | 4 | 0.071 | 0.114 | 1/4 | 0.059 | 0.098 | 28 | 6 | 1 / 1 |
| 9 | both, sign test | 10 | | | 1/10, p = 0.021 | | | | | |
| 50 | K3 (A 20 pJ/m) | 6 | 0.237 | 0.113 | 5/6 | 0.154 | 0.083 | 24 | 2 | 0 / 1 |
| 50 | base (A 10 pJ/m) | 4 | 0.230 | 0.114 | 3/4 | 0.198 | 0.109 | 22 | 2 | 0 / 1 |
| 50 | both, sign test | 10 | | | 8/10, p = 0.109 | | | | | |
| 100 | K3 (A 20 pJ/m) | 6 | 0.215 | 0.073 | 6/6 | 0.138 | 0.042 | 26 | 2 | 0 / 0 |
| 100 | base (A 10 pJ/m) | 4 | 0.166 | 0.082 | 4/4 | 0.138 | 0.076 | 24 | 3 | 1 / 0 |
| 100 | both, sign test | 10 | | | 10/10, p = 0.002 | | | | | |
