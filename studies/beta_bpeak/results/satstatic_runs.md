# Static saturation test: per-run checks (night 2026-10-19/20, working log for PROTOKOLL 7.128)

Code ebcea74, spool 25_satstatic.jobs. w_l = mean w_loop of the last 5 cycles of the 9 mT stage (as 7.125).
Reference K4_1 (7.125): nx 1.73e-7, AC from saturation 1.35e-7, AC 300 mT (option C) 5.24e-7.

| run | done (UTC) | static relax (it, conv) | BPs after static | passes | BPs end pass 0 | relax2 it | E/(K_d V) prepared | 9 mT w_l | cycles, calm tail | BPs in cycles |
|---|---|---|---|---|---|---|---|---|---|---|
| K4_1_nx_sbd | 00:47 | 2500, yes (ramp max 88 BPs) | 0 | 1 | 0 (max 94 in AC) | 3000 | -3.6317e-4 | 1.310e-7 | 8, 5 | 0 |
| K4_1_nx_sad | 00:51 | 2000, yes (no ramp, uniform start) | 0 | 1 | 0 (max 94 in AC) | 6500 | -3.7029e-4 | 2.508e-7 | 8, 5 | 0 |
| K4_1_nx_sa200 | 01:23 | 2500, yes (no ramp, uniform start) | 0 | 1 (41 half cycles) | 0 (max 996 in AC; 115 at the end of half 0) | 9500 | -3.7034e-4 | 2.251e-7 | 8, 5 | 0 |

| K4_1_nx_sb200 | 01:25 | 3500, yes (ramp max 363 BPs) | 0 | 1 (41 half cycles) | 0 (max 996 in AC) | 7500 | -3.7150e-4 | 4.015e-7 | 8, 5 | 0 |
Notes K4_1, h_sat 1 - φ (00:56):
- Static states: E/(K_d V) -0.4166043424 (A) and -0.4166043405 (B), difference 1.9e-9 (relative 5e-9); alignment identical
  to 3 digits.
- AC pass 0 (35 half cycles from h 0.45): M_par peak and BP count per half cycle are identical for A and B in halves 0 ... 5
  (BPs at the end 19 18 6 7 0 0). They separate at half 6 (h 0.118, the first large drop of M_par 0.536 -> 0.469):
  BPs at the end 2 (A) and 1 (B); afterwards A 0 1 1 0 ..., B 1 1 2 2 1 0 ...
- 9 mT: A 2.508e-7, B 1.310e-7, ln(A/B) +0.649 -> rule of 7.127: "not erased" (BPs after pass 1: 0 and 0, same).
- The energy difference 1.9e-9 is below the fp32 resolution (eps 1.2e-7; the runs use --precision single). Thus the two
  saturated states are the same state within the precision of the model. Known rounding-level differences between A and
  B: other GPU (A gpu 2, B gpu 1) and another integrator history before the AC pass (B had the 64-step ramp, A not; the
  adaptive time step differs).
- Marker: A has max_angle 85 deg in cycles 0 ... 2 and 93 deg from cycle 3 on, and w_loop still rises by about 0.3 % per
  cycle (dW 0.0026 at c7, within the calm rule); B stays at 85 deg. A changed its state during the measurement, B not.
- Candidate mechanism (inference, n 1 pair; hold until the 2 M_s pair is done): the AC pass selects one of several states
  at the first irreversible switching (half 6), from rounding-level differences. Discriminator: sa200 against sb200
  (0 cells > 15 deg at 2 M_s). If that pair also separates, the start does not matter and the spread belongs to the AC
  path; if it agrees, the 1 - φ field left a memory in the 1.4 % tilted cells.

Notes K4_1, h_sat 2 M_s (01:30):
- Static states: E/(K_d V) -2.121323355 (A) and -2.121323335 (B), difference 2.0e-8 (relative 1e-8, below fp32 eps).
- AC pass 0 (41 half cycles from h 2): M_par peak and BP count per half cycle identical for A and B in halves 0 ... 12
  (BPs at the end 115 83 56 51 23 26 29 17 18 2 2 6 1); they separate at half 13 (h 0.110, the first large drop of
  M_par 0.531 -> 0.437): 1 (A) and 0 (B). Thus at both fields the paths separate at the same h (0.11 ... 0.12).
- 9 mT: A 2.251e-7, B 4.015e-7, ln(A/B) -0.579 -> "not erased" (BPs after pass 1: 0 and 0).
- The sign of ln(A/B) changes between the two fields (+0.649 at 1 - φ, -0.579 at 2 M_s): no systematic memory of the start.
- Question 2 (field), same start: A ln(2 M_s / 1 - φ) -0.108, B +1.120 -> both "matters" by the rule (A just above 0.10).
- K4_1 overall: four 9 mT losses 1.31 / 2.25 / 2.51 / 4.02 e-7 from two saturated states; with 7.125: 1.31 ... 5.24 e-7.

K2_0 (h_sat 1 - φ = 0.38, 2 M_s):

| run | done (UTC) | static relax (it, conv) | BPs after static | passes | BPs end pass 0 | relax2 it | E/(K_d V) prepared | 9 mT w_l | cycles, calm tail | BPs in cycles |
|---|---|---|---|---|---|---|---|---|---|---|
| K2_0_nx_sad | FAIL 04:50 (exit code 3) | 80000, no (absolute criterion below the fp32 noise of E: E flat within 1e-9 K_d V, torque rms 0.03 A/m) | 0 | 1 | 0 | – | – | – | – | – |
| K2_0_nx_sa200 | 04:5x | 8500, yes | 0 | 1 (41 half cycles) | 0 (max 1292 in AC) | 5000 | -7.4780e-3 | 1.350e-7 | 7, 5 | 0 |

- 04:50 the dispatcher started the gate retry of K2_0_nx_sad with the old code ebcea74 and relax_maxiter 200000. The
  patch 4bb78f1 (static relaxation with the relative criterion, Chris "Ok" 04:3x) is committed and pushed; the switch
  of the dispatcher worktree on the box was blocked by the permission check. Stop of the retry and restart: needs the
  go of Chris.
- 04:5x Chris: "Ok" to the patch 4bb78f1 (the static relaxation in the field uses the relative criterion relax_tol_E
  1e-6, 2 chunks in a row; the zero-field relaxations keep the absolute threshold 3e-10), then "Go" to stop the old
  retry and restart. 04:54 dispatcher worktree switched to 4bb78f1 (/root/disp_commit); 04:55 the retry (PID 109449,
  old code, in its setup) stopped with SIGTERM; the dispatcher put it back (exit code -15, attempt 1 of 1) and started
  it at once on gpu 3: fresh start (no checkpoint), code 4bb78f1, relax_maxiter 200000 (the gate-retry value; it
  matters only if a zero-field relaxation needs more than 80000 iterations).
