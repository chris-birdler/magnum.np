# PLAN 4.12 V1: Newell candidates against the alias sum (check_demag_alias.py v1_report N0b)

Energies per box volume in K_d. Reference: fit c0 + c1 / Gm + c2 / Gm² over Gm = 8, 16, 32; u: difference to the linear extrapolation from 16 and 32. dE = (E − E_ref) / E_ref; dH = RMS over the magnetic cells of |H − H_ref| / RMS |H_ref|. Largest trace error of the alias kernels: 9.2e-15.

## Calibration

| case | E(8) | E(16) | E(32) | E_ref | u | exact | E_ref / exact − 1 | dH vs exact |
|---|---|---|---|---|---|---|---|---|
| cube f=1/8 (111) | 0.0362724 | 0.0363625 | 0.0364096 | 0.0364582 | 3.9e-05 | 0.0364583 | -0.0002 % | – |
| checkerboard (1,-1,0) | 0.4763894 | 0.4877866 | 0.4937833 | 0.4999788 | 4.0e-04 | 0.5000000 | -0.0042 % | 0.0042 % |
| B0b geometry uniform (111) | 0.0757233 | 0.0758929 | 0.0759821 | 0.0760742 | 3.8e-05 | 0.0760744 | -0.0003 % | – |

**Calibration: reference calibrated.**

## Candidates (12 states)

| candidate | max dE | max (dE + 2u) | max dH | max (dH + 2u_H) | reading |
|---|---|---|---|---|---|
| 1 | -14.9142 % | 15.0048 % | 5.3022 % | 5.3772 % | not valid |
| 2 | -5.2463 % | 5.3369 % | 1.8666 % | 1.9416 % | not valid |
| 4 | -1.6189 % | 1.7095 % | 0.5767 % | 0.6517 % | not valid |
| 8 | -0.4515 % | 0.5421 % | 0.1615 % | 0.2364 % | not valid |
| 2x | +0.1918 % | 0.2824 % | 0.0666 % | 0.1416 % | not valid |
| 4x | +0.0063 % | 0.0967 % | 0.0039 % | 0.0791 % | not determined |

**Final setting: none (stop and report).**
