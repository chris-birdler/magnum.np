# Tolerance twins of the main run at 30 MHz (reading fixed in PROTOKOLL 7.101)

c: median closure |M(end) − M(start)| / M_peak over the measured cycles (2 … 6) of all 7 stages.

| seed | c atol 1e-5 (K1_k) | c atol 1e-6 (K1_k_L6) | r | start E/K_d rel. diff. | start max abs. diff. of M/M_s | w_loop cycle 0 at 9 mT, L6 / 1e-5 − 1 | control |
|---|---|---|---|---|---|---|---|
| 0 | 0.000547 | 0.000545 | 1.00 | 0.00e+00 | 1.82e-09 | +4.83e-06 | same start |
| 1 | 0.000570 | 0.000570 | 1.00 | 0.00e+00 | 4.92e-09 | +2.62e-05 | same start |
| 2 | 0.000539 | 0.000538 | 1.00 | 0.00e+00 | 4.38e-09 | -2.01e-05 | same start |

**Reading: physical** (r ≥ 0.5 physical, r < 0.5 numerical; at least 2 of 3 seeds)

Per B, mean over the 3 complete seeds (reported, not part of the reading):

| B (mT) | median closure 1e-5 | median closure 1e-6 | median w_dis / w_loop 1e-5 | median w_dis / w_loop 1e-6 | P 1e-5 (W/cm³) | P 1e-6 (W/cm³) | P ratio 1e-6 / 1e-5 |
|---|---|---|---|---|---|---|---|
| 9 | 0.000012 | 0.000012 | 1.0209 | 1.0008 | 2.29 | 2.29 | 1.000 |
| 20 | 0.000114 | 0.000114 | 0.9932 | 0.9904 | 18.6 | 18.6 | 1.000 |
| 35 | 0.000244 | 0.000244 | 1.0110 | 1.0103 | 69.2 | 69.2 | 1.000 |
| 50 | 0.000416 | 0.000416 | 1.0060 | 1.0056 | 136 | 136 | 1.000 |
| 70 | 0.001175 | 0.001175 | 0.9983 | 0.9982 | 359 | 359 | 1.000 |
| 100 | 0.001227 | 0.001229 | 1.0374 | 1.0373 | 972 | 972 | 1.000 |
| 150 | 0.001990 | 0.001863 | 0.9966 | 1.0008 | 2.2e+03 | 2.18e+03 | 0.992 |

## Exploratory (not pre-registered; added 2026-10-15 after the reading)

Event cycles per stage (measured cycles): the same with both tolerances in 21 of 21 stages.

β per seed (main-fit definition: ln w at the nominal amplitudes, steady-tail rule; n 3, no effect claim):

| seed | β 9 … 150 mT, 1e-5 / 1e-6 | β 9 … 35 mT | β 35 … 70 mT | β 70 … 150 mT |
|---|---|---|---|---|
| 0 | 2.3111 / 2.3062 | 2.5706 / 2.5706 | 2.0665 / 2.0665 | 2.2770 / 2.2540 |
| 1 | 2.5532 / 2.5532 | 2.4089 / 2.4089 | 2.8585 / 2.8585 | 2.5804 / 2.5803 |
| 2 | 2.4816 / 2.4797 | 2.4647 / 2.4647 | 2.4728 / 2.4728 | 2.2941 / 2.2821 |

R1 (share of the loss at the drive frequency), lag of the fundamental against the drive and loss shares per zone (acc.json, calm cycles), mean over the seeds; 1e-5 / 1e-6:

| B (mT) | R1 | lag (deg) | largest abs. change of a zone share (zones r<0.5, 0.5<=r<0.8, r>=0.8, texture top 10%, p top 10%) |
|---|---|---|---|
| 9 | 0.9710 / 0.9904 | -0.17 / -0.17 | 0.0135 |
| 20 | 0.8039 / 0.8058 | -0.29 / -0.29 | 0.0020 |
| 35 | 0.6157 / 0.6162 | -0.35 / -0.35 | 0.0007 |
| 50 | 0.5591 / 0.5593 | -0.34 / -0.34 | 0.0003 |
| 70 | 0.4980 / 0.4981 | -0.45 / -0.45 | 0.0001 |
| 100 | 0.3606 / 0.3606 | -0.60 / -0.60 | 0.0001 |
| 150 | 0.3657 / 0.3655 | -0.61 / -0.60 | 0.0025 |
