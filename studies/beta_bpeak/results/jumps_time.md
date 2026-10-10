# Stage 1: jumps in the dissipation time series (exploratory)

n = 115 α 0.1 runs (Sobol, centre, corners), the last 3 cycles per amplitude if they are steady (7.54), 256 samples per cycle. Medians over the runs, 95 % CI of the median by bootstrap over the runs (2000 resamples).

## A. Measures against the amplitude

| B (mT) | RF non-smooth | RN non-repeating | RH sharp repeating | dS5 burst excess | rho event coupling (runs > 0) | samples with an n60 change | R_np snapshots | left out, not steady (runs) |
|---|---|---|---|---|---|---|---|---|
| 9 | 0.040 (0.025 … 0.055) | 0.006 (0.005 … 0.007) | 0.028 (0.021 … 0.042) | 0.004 | 0.01 (0.00 … 0.02), 65/104 | 17 % | 0.000 | 4 |
| 20 | 0.133 (0.111 … 0.164) | 0.020 (0.006 … 0.038) | 0.112 (0.096 … 0.132) | 0.027 | 0.04 (0.02 … 0.06), 75/100 | 33 % | – | 6 |
| 35 | 0.218 (0.193 … 0.242) | 0.091 (0.023 … 0.135) | 0.167 (0.159 … 0.179) | 0.044 | 0.08 (0.06 … 0.10), 80/103 | 49 % | – | 11 |
| 50 | 0.230 (0.210 … 0.240) | 0.121 (0.100 … 0.143) | 0.181 (0.168 … 0.193) | 0.044 | 0.13 (0.10 … 0.14), 98/110 | 60 % | 0.029 | 2 |
| 70 | 0.216 (0.202 … 0.224) | 0.116 (0.103 … 0.126) | 0.165 (0.155 … 0.176) | 0.035 | 0.13 (0.11 … 0.14), 101/112 | 69 % | – | 3 |
| 100 | 0.184 (0.175 … 0.190) | 0.117 (0.109 … 0.130) | 0.139 (0.130 … 0.144) | 0.027 | 0.10 (0.08 … 0.11), 105/112 | 76 % | – | 3 |
| 150 | 0.153 (0.146 … 0.162) | 0.101 (0.093 … 0.111) | 0.109 (0.103 … 0.113) | 0.018 | 0.07 (0.06 … 0.08), 109/115 | 81 % | 0.025 | 0 |

Check: RF median with the steady rule against all run-amplitudes (the last 3 cycles, also when not steady):

| B (mT) | RF median, steady | RF median, all |
|---|---|---|
| 9 | 0.040 (n 111) | 0.041 (n 115) |
| 20 | 0.133 (n 109) | 0.139 (n 115) |
| 35 | 0.218 (n 104) | 0.230 (n 115) |
| 50 | 0.230 (n 113) | 0.231 (n 115) |
| 70 | 0.216 (n 112) | 0.217 (n 115) |
| 100 | 0.184 (n 112) | 0.185 (n 115) |
| 150 | 0.153 (n 115) | 0.153 (n 115) |

## B. Amplitude B_pk of the largest RF per run

All runs: median B_pk 50 mT (95 % CI 42 … 52), n 115. B_pk at 9 or 150 mT (end of the range): 6 runs.

Linear model ln B_pk = c0 + Σ c_j x_j (x_j coded −1 … +1, OLS, HC3 standard errors). Range effect = factor on B_pk from the low end to the high end of the factor, with 95 % CI:

| factor (low → high end) | factor on B_pk (95 % CI) | p |
|---|---|---|
| L_eff 30 → 12 l_ex (K_eff up) | 1.12 (0.85 … 1.49) | 0.423 |
| r_p 0 → 3 | 1.20 (0.86 … 1.67) | 0.283 |
| φ 0.55 → 0.69 | 0.86 (0.61 … 1.22) | 0.403 |
| d/l_ex 212 → 300 (A 20 → 10 pJ/m) | 0.67 (0.51 … 0.88) | 0.004 |

R² = 0.07. Groups (median B_pk, 95 % CI):

| group | n | median B_pk (mT) |
|---|---|---|
| A low third (9.7 … 12.5 pJ/m) | 39 | 39 (34 … 46) |
| A mid third | 38 | 47 (36 … 52) |
| A high third (15.7 … 20.0 pJ/m) | 38 | 54 (51 … 60) |
| L_eff 12 l_ex | 26 | 52 (34 … 70) |
| L_eff 18 l_ex | 38 | 50 (38 … 53) |
| L_eff 24 l_ex | 24 | 39 (31 … 52) |
| L_eff 30 l_ex | 27 | 50 (41 … 53) |

## C. Mesh check: the same measures at dx 3 and dx 1.5 l_ex (paired seeds, one job per amplitude)

Jobs with 5 cycles, the last 3 used. A pair counts only if both jobs are steady in these 3 cycles (7.54). Medians over the seeds. n60 max = largest number of pairs > 60° in the used cycles. Sign test (two-sided, exact) over the seeds of both points.

| B (mT) | point | n pairs | RF dx 3 | RF dx 1.5 | seeds with RF(dx 3) > RF(dx 1.5) | RN dx 3 | RN dx 1.5 | n60 max dx 3 | n60 max dx 1.5 | pairs left out (not steady) |
|---|---|---|---|---|---|---|---|---|---|---|
| 9 | K3 (A 20 pJ/m) | 5 | 0.035 | 0.114 | 0/5 | 0.035 | 0.031 | 20 | 20 | 1 |
| 9 | base (A 10 pJ/m) | 2 | 0.071 | 0.098 | 0/2 | 0.059 | 0.080 | 30 | 12 | 2 |
| 9 | both, sign test | 7 | | | 0/7, p = 0.016 | | | | | |
| 50 | K3 (A 20 pJ/m) | 4 | 0.239 | 0.113 | 4/4 | 0.163 | 0.082 | 24 | 2 | 2 |
| 50 | base (A 10 pJ/m) | 1 | 0.189 | 0.114 | 1/1 | 0.145 | 0.109 | 22 | 2 | 3 |
| 50 | both, sign test | 5 | | | 5/5, p = 0.062 | | | | | |
| 100 | K3 (A 20 pJ/m) | 6 | 0.215 | 0.073 | 6/6 | 0.138 | 0.042 | 26 | 2 | 0 |
| 100 | base (A 10 pJ/m) | 3 | 0.165 | 0.085 | 3/3 | 0.137 | 0.081 | 23 | 3 | 1 |
| 100 | both, sign test | 9 | | | 9/9, p = 0.004 | | | | | |
