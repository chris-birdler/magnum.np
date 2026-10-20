# Noise floor of the energy in the converge relaxations (eval_relax_noise.py)

Threshold 3e-10 K_d V per chunk (PLAN 4.13 S7, D1 b). A rise of E is noise (alpha 1, no field).

| run | relaxation | iterations | E/(K_d V) | rises | rise median | rise 90 % | tail median \|dE\| | quiet share on plateau | cells | BPs end |
|---|---|---|---|---|---|---|---|---|---|---|
| K1_0_ad | virgin | 8000 | 1.695e-03 | 0 | nan | nan | 3.3e-10 | nan | 1.26e+06 | 0 |
| K1_0_ad | relax2_p0 | 1500 | 1.684e-03 | 0 | nan | nan | 1.5e-10 | nan | 1.26e+06 | 0 |
| K1_1_ad | virgin | 12000 | 1.815e-03 | 3 | 3.7e-10 | 5.2e-10 | 3.7e-10 | 0.40 | 1.26e+06 | 0 |
| K1_1_ad | relax2_p0 | 4500 | 1.710e-03 | 2 | 3.0e-10 | 4.6e-10 | 5.1e-10 | 0.40 | 1.26e+06 | 0 |
| K1_2_ad | virgin | 14000 | 1.748e-03 | 5 | 6.7e-10 | 1.4e-09 | 3.7e-10 | 0.36 | 1.26e+06 | 1 |
| K1_2_ad | relax2_p0 | 5000 | 1.726e-03 | 5 | 4.1e-10 | 1.3e-09 | 2.2e-10 | 0.25 | 1.26e+06 | 0 |
| K2_0_ad | virgin | 32500 | -7.431e-03 | 2 | 1.9e-10 | 3.3e-10 | 1.3e-09 | 0.67 | 3.51e+06 | None |
| K2_0_nx_sa200 | relax2_p0 | 5000 | -7.478e-03 | 4 | 7.6e-10 | 1.5e-09 | 1.1e-10 | 0.33 | 3.51e+06 | 0 |
| K2_1_ad | virgin | 33000 | -7.435e-03 | 8 | 4.2e-10 | 4.7e-10 | 4.2e-10 | 0.29 | 3.51e+06 | None |
| K2_2_ad | virgin | 25000 | -7.346e-03 | 1 | 6.1e-10 | 6.1e-10 | 6.6e-10 | 0.33 | 3.51e+06 | 6 |
| K2_2_ad | relax2_p0 | 1500 | -7.518e-03 | 2 | 1.3e-10 | 1.9e-10 | 1.3e-10 | 1.00 | 3.51e+06 | 0 |
| K2_2_nx | virgin | 43000 | -7.486e-03 | 4 | 4.7e-10 | 8.5e-10 | 3.1e-09 | 0.25 | 3.51e+06 | None |
| K3_0_ad | virgin | 22500 | 4.657e-04 | 12 | 4.0e-10 | 1.4e-09 | 4.9e-10 | 0.28 | 1.12e+06 | 0 |
| K3_0_ad | relax2_p0 | 7500 | 4.562e-04 | 8 | 5.3e-10 | 1.1e-09 | 2.7e-10 | 0.23 | 1.12e+06 | 0 |
| K3_1_ad | virgin | 19000 | 5.242e-04 | 13 | 7.1e-10 | 1.7e-09 | 5.5e-10 | 0.17 | 1.12e+06 | 2 |
| K3_1_nx | virgin | 14500 | 7.247e-04 | 4 | 5.2e-10 | 1.1e-09 | 9.0e-10 | 0.33 | 1.12e+06 | None |
| K3_2_ad | virgin | 25500 | 4.360e-04 | 9 | 4.4e-10 | 1.1e-09 | 3.5e-10 | 0.35 | 1.12e+06 | 1 |
| K3_2_ad | relax2_p0 | 6000 | 4.110e-04 | 6 | 4.3e-10 | 1.4e-09 | 2.7e-10 | 0.20 | 1.12e+06 | 0 |
| K4_0_ad | virgin | 16500 | -3.385e-04 | 2 | 1.4e-10 | 1.5e-10 | 2.2e-09 | 0.75 | 3.80e+06 | 1 |
| K4_0_ad | relax2_p0 | 11500 | -3.718e-04 | 9 | 4.7e-10 | 9.5e-10 | 1.7e-10 | 0.26 | 3.80e+06 | 0 |
| K4_0_nx | virgin | 19500 | -2.617e-04 | 5 | 2.2e-10 | 3.1e-10 | 3.3e-10 | 0.67 | 3.80e+06 | None |
| K4_1_ad | virgin | 13000 | -3.291e-04 | 4 | 2.4e-10 | 4.8e-10 | 4.4e-10 | 0.40 | 3.80e+06 | None |
| K4_1_nx_sa200 | relax2_p0 | 9500 | -3.703e-04 | 10 | 2.8e-10 | 5.2e-10 | 3.1e-10 | 0.35 | 3.80e+06 | 0 |
| K4_1_nx_sad | relax2_p0 | 6500 | -3.703e-04 | 3 | 4.7e-10 | 9.9e-10 | 1.5e-10 | 0.44 | 3.80e+06 | 0 |
| K4_1_nx_sb200 | relax2_p0 | 7500 | -3.715e-04 | 5 | 1.1e-09 | 2.0e-09 | 3.4e-10 | 0.30 | 3.80e+06 | 0 |
| K4_1_nx_sbd | relax2_p0 | 3000 | -3.632e-04 | 2 | 4.3e-10 | 6.3e-10 | 3.9e-10 | 0.67 | 3.80e+06 | 0 |
| K4_2_ad | virgin | 48000 | -3.486e-04 | 34 | 7.9e-10 | 1.7e-09 | 6.2e-10 | 0.13 | 3.80e+06 | 2 |
| K4_2_ad | relax2_p0 | 10000 | -3.853e-04 | 8 | 1.0e-09 | 1.8e-09 | 6.4e-10 | 0.18 | 3.80e+06 | 0 |

Relaxations: 28 (with at least 3 rises: 20).
Rise median: min 2.2e-10, median 4.7e-10, max 1.1e-09; rise 90 %: median 1.1e-09, max 2.0e-09.
Rise median / threshold: median 1.58, max 3.62; relaxations with the rise median above the threshold: 17.
Quiet share on the plateau: min 0.13, median 0.29.
Correlation of log rise median with log |E|: 0.11.
Correlation of log rise median with log cells: 0.01.
Correlation of log rise median with log iterations: -0.09.
With BPs at the end: n 4, rise median 6.9e-10; without: n 11, rise median 4.7e-10.
