# PLAN: β(B_peak) production study (draft for approval)

Status 2026-10-06: plan v2 (4.0) prepared and tested (PROTOKOLL 7.30,
7.33); waiting for Chris's go and a credit top-up. The target is the total
loss at 30 MHz with the model α (Chris 2026-10-03, PROTOKOLL §12); the
history of steps 1 … 1d is in §3 … 3f.
Basis: PROTOKOLL.md (model, tests, mesh decision).

## 1. Question

Which parameters control the Steinmetz exponent β of a powder core of 1 µm
particles at B_peak = 10 … 150 mT, and how can β stay small?

Output: β_core of the mean loss curve over all amplitudes and in three
windows (centres ≈ 18, 50 and 102 mT; the lowest window starts at 9 mT, so
β at 10 mT itself is not delivered) and the loss density E[w] (P in W/cm³ at
f = 30 MHz) at 10, 50 and 100 mT as functions of four factors, with errors;
a classification of the factors (PLAN 4.0, decision rules).

## 2. Fixed settings (from PROTOKOLL)

| Setting | Value | Reason |
|---|---|---|
| model | Herzer cubes (Q_eff, L_eff) + particle-scale stress r_p; FCC cell of 4 spheres | PROTOKOLL §2 |
| mesh | dx = 3 l_ex (10 nm at A = 10 pJ/m; 10 … 14.2 nm in reading 1 of PLAN 4.0) | no mesh effect on β detected at 1 µm (95 % CI ≈ ±0.2 … 0.4); level at dx 3 higher than at dx 1.5 in 3 of 4 seeds (+21 … +43 %, not significant), direction to the converged value unknown (PROTOKOLL 7.28, 7.32) |
| initial state | virgin state (random m on 6 l_ex blocks, relaxed), ascending amplitudes | PROTOKOLL §7.4 … 7.6 |
| damping | α = 0.1 (8 centre runs at α 0.01) | model parameter (PROTOKOLL 7.17, 7.31) |
| precision | fp32, atol 10⁻⁵ | numerical floor test PASS |
| amplitudes | B_peak = 9 / 20 / 35 / 50 / 70 / 100 / 150 mT (D2), 7 cycles at α 0.1 (9 at α 0.01), cycles 0 and 1 discarded | study range; β_core in the windows 9-20-35, 35-50-70, 70-100-150 mT and over all amplitudes (4.0) |
| frequency | f = 30 MHz (f/f_M = 7.14·10⁻⁴) | target is the total loss with the model α (Chris 2026-10-03, PROTOKOLL §12) |

## 3. Step 1: frequency test (go 1)

Question: which part of the loss at 30 MHz is dynamic (damping), and is
β ≈ 2 a damping artefact?

| Runs | Content |
|---|---|
| F1: 4 seeds at f | base point: d/l_ex = 300, L_eff = 12, r_p = 0, φ = 0.65, dx = 3 |
| F2: the same 4 seeds at 2 f | same arguments, f/f_M doubled |
| F3, F4 (option D4): 4 seeds at f and 2 f | low-loss corner L_eff = 30, r_p = 3 (loss ×0.24 … 0.31): the dynamic share is largest there |

Evaluation: w(f) = w_h + c·f gives w_h = 2 w(f) − w(2 f) and the dynamic
share 1 − w_h/w(f) per amplitude (mean ± SE over 4 seeds); β from w and
from w_h.

Decision rule (fixed before the runs):
- the 90 % upper bound of the dynamic share is < 30 % at all amplitudes,
  and the 90 % upper bound of |β(w) − β(w_h)| is < 0.2: step 2 at f only.
  With σ(ln w) ≈ 0.2 and 4 seeds the SE of the dynamic share is ≈ 0.14;
  if the test is not decisive, more seeds (≈ 1 $ each pair) before step 2.
- else: options for Chris (pairs f / 2 f for every design point, +50 %
  cost; or a lower f, more cost per cycle).

Cost (measured: 0.059 h per cycle at d/l_ex = 300, dx = 3, V100; the time
per cycle does not depend on the amplitude): 7 amplitudes × 7 cycles = 49
cycles: ≈ 3.0 h per run at f, ≈ 1.6 h at 2 f (assumption: half the steps
per cycle, not measured). 16 runs (F1 … F4, 4 seeds each) ≈ 37 V100-h
≈ 5 $. Wall time ≈ 5 h on 8 GPUs, ≈ 9 h on 4 GPUs.
Seeds 901 … 904 (base) and 905 … 908 (corner). The F1 and F3 runs are also
used as design points of step 2.

## 3b. Step 1b: mesh error at the small end of the d/l_ex range

Question: the staircase surface of the voxelized spheres has a relative
roughness dx/d = 1.0 % at d/l_ex = 300 and 1.4 % at 212. At fixed dx this
trend is confounded with the factor d/l_ex. If the mesh error (dx 3 vs
dx 1.5) is the same at both ends, it cancels in the d/l_ex effect.

| Runs | Content |
|---|---|
| R212_dx15 (4 seeds 911 … 914) | d/l_ex = 212, dx = 1.5, base point, 9 and 50 mT |
| R212_dx3 (same 4 seeds) | d/l_ex = 212, dx = 3 |

Reference: d/l_ex = 300, dx 3 vs dx 1.5 (planning pilot): level +0.23 ±
0.07 at 50 mT (ln w), β −0.20 ± 0.23.
Rule: if the level difference dx 3 − dx 1.5 at d/l_ex = 212 agrees with the
value at 300 within ± 0.10 in ln w (90 %), the d/l_ex effect on the level is
clean; else it is corrected by the measured mesh share or Chris decides on a
smoother surface. β (9 → 50 mT) only coarse (≈ ± 0.3).
Cost: ≈ 31 V100-h ≈ 4 … 5 $ (9 M cells at dx = 1.5), ≈ 8 h on 4 GPUs.

## 3c. Step 1c: surface roughness (smooth surface at dx = 3)

Step 1b showed that the mesh error is not the same at both ends of the
d/l_ex range, but not why. This test changes only the surface: surface cells
carry their magnetic volume fraction f (`--surface fraction`; Ms, A, K_eff
and K_p scaled with f; φ exact). Runs: d/l_ex 212 (seeds 911 … 914) and 300
(seeds 1 … 3, one job per amplitude), 9 and 50 mT, all else as the
staircase references (checked: only `surface` differs).
- If the smooth dx = 3 moves to dx = 1.5 and the d-dependent difference
  disappears: the staircase is the cause, and the smooth surface removes it
  without extra cost (a model change: decision by Chris).
- Else: the cause is in the interior (unresolved structures, walls).
Note: the partial cells need a gap of ≈ 3 cells; with this surface φ ≤ 0.66
at d/l_ex = 212 … 300 (φ = 0.69 gives contact).
Cost ≈ 4.5 V100-h ≈ 0.5 $, chained after the α test.

## 3d. Step 1d: relaxation of the virgin state to convergence

Finding (2026-10-03): at d/l_ex = 300 the 9 mT stage is not steady in all
earlier runs: w_dis/w_loop = 2 … 76 (α = 0.1) and 17 … 44 (α = 0.01); at
d/l_ex = 96 it is 1.1 … 1.3. Cause: the relaxation of the virgin state stops
after 5000 iterations without convergence; the stop criterion (max torque <
100 A/m) cannot be met at dx = 3 because the unresolved structures keep a
large local torque. The particle keeps relaxing during the 9 mT cycles.
Consequence: all 9 mT values at 1 µm so far are not reliable (and β at
10 mT). From 20 mT the balance is ≈ 1.
Fix: `--relax_mode converge` (chunks of 500 iterations at α = 1 until the
relative energy change per chunk < 10⁻⁶; RMS torque logged). CPU check:
converged, then w_dis/w_loop = 1.01 at the first amplitude.
Test: d/l_ex 300, dx 3, base point, α 0.1 and 0.01, seeds 901 … 903, 9 and
20 mT; rule: w_dis/w_loop at 9 mT within 0.9 … 1.1 in the kept cycles.
Cost ≲ 13 V100-h ≈ 1.5 $, chained after step 1c.

## 3e. Repeats after the audit (PROTOKOLL 7.13, go 2026-10-03, option B)

All repeats use `--relax_mode converge` with the hard gate (a run stops with
exit code 3 if the relaxation does not converge; analyze.py excludes such
runs) and save the relaxed virgin state as `init.pt`.
- `alphatest2`: α 0.01, seeds 904 … 906, now with the converged relaxation
  (A010c_*).
- `baseline`: realisation scatter at 1 µm: PROD (d/l_ex 300, dx 3,
  7 amplitudes), base point, α 0.1, seeds 921 … 928 (≈ 27 V100-h). It replaces
  σ at 9 mT in the planning of the number of realisations.
- `mesh300`: dx 1.5 vs dx 3 at d/l_ex 300, seeds 921 … 924, 9 / 50 / 100 mT,
  5 cycles; the 9 mT job relaxes and writes init.pt, the 50 and 100 mT jobs
  start from it (≈ 183 V100-h). Compared with the valid 50 mT value at
  d/l_ex 212 (PROTOKOLL 7.10, rule of step 1b).
- Later decision: the same at d/l_ex 212 (≈ 47 V100-h) if the 9 mT mesh error
  at 1 µm is to be transferred to the small end of the range.
Total ≈ 230 V100-h ≈ 22 $ on one 4×V100 instance (0.41 $/h), ≈ 2.5 days.

## 3f. Domain analysis (prepared 2026-10-03, no GPU cost)

Question (Chris): does β depend on the total wall area? Tool `domains.py`,
tested on synthetic states with known walls (PROTOKOLL 7.14):
- wall measure T(2 L_eff) (rotation of m on scales above 2 L_eff);
- switched volume between two phases of a cycle.
A segmentation into domains was tested and archived (it over-counts walls
2 … 8×). Use: first on the snapshots (snaptest) and on init.pt of `baseline`;
only if this gives a clear signal, the measures go into the runs as a per-cycle
time series (decision by Chris).

## 4. Step 2: randomized design (go 2)

### 4.0 Plan v2 (2026-10-06, for approval by Chris)

Supersedes 4.1 … 4.4 and §5 where they differ. Basis: PROTOKOLL 7.13 … 7.33
(7.33: decisions on the third audit, one by one with Chris).

**Terms.**
- Ms = Js/μ₀ saturation magnetisation; l_ex = √(2A/(μ₀Ms²)) exchange length; K_d = μ₀Ms²/2; reduced lengths are
  in units of l_ex. d particle diameter; L_eff edge of the Herzer cubes
  (random-anisotropy cells); Q_eff = (l_ex/L_eff)² = K_eff/K_d; r_p = K_p/K_eff
  particle anisotropy (residual stress) relative to K_eff; φ packing fraction.
- Run: one simulation of one design point with one seed (its own cubes,
  axes and initial state). Job: one line of a job file = one call of
  run_loops.py; in the design one job = one run; in the mesh300 protocol one
  run = 3 jobs (one per amplitude). DONE / POSTPROC_DONE / POSTPROC_FAILED:
  marker files in the run folder (run finished / analysed on the instance /
  analysis failed). Realisation scatter: the scatter of a result over
  seeds at the same factors.
- init.pt: the relaxed virgin state of a run (start of the amplitude
  ramp); checkpoint.pt: the state after the last finished amplitude stage.
- w: loss per cycle and core volume (J/m³), w = mean of w_loop = ∮H dB over
  the kept cycles; P = f·w (W/cm³). E[w]: mean of w over the realisations
  (what a core of ≈ 10⁹ particles shows). B̂ = B_peak, the peak flux
  density of the core (mT). β_core = d ln E[w]/d ln B̂. w_dis: the same loss
  from the LLG dissipation (check of w_loop). μ: relative permeability of
  the core (≈ 8, PROTOKOLL §0).
  ln E[w] is the loss output in all rules below.
- Window: 3 neighbouring amplitudes (9-20-35, 35-50-70, 70-100-150 mT;
  centres ≈ 18 / 50 / 102 mT); β in a window = slope of a power law through
  its 3 points. "β all": over all 7 amplitudes.
- Range effect of a factor: prediction of the main fit at the factor
  maximum minus at its minimum, the other factors averaged over the design
  points.
- r: distance of a cell from its particle centre, R = d/2 particle radius.
- x = (ln Q_eff, r_p, φ, ln d/l_ex): the regression variables.
- p: local LLG dissipation per cell; p_dis: its box mean; p_fund: the part
  of the fundamental (drive frequency); R1 = p_fund/p_dis; lag: phase of the
  local fundamental against the drive; texture: |∇m| of the cycle-mean state
  smoothed over 0.5 L_eff (top 10 % cells = "texture top 10 %"); F90:
  fraction of the volume rotated by > 90° (PROTOKOLL 7.14); l_C:
  correlation length of m.
- Gamma GLM: regression of w with a log link and a Gamma error model (it
  estimates ln E[w]); OLS: ordinary least squares; HC3: heteroscedasticity-
  consistent standard error; Huber fit: robust regression that down-weights
  outliers; Holm: Holm–Bonferroni step-down correction for several tests;
  TOST: two one-sided tests (equivalence test).
- Sobol sequence: a quasi-random set of points that fills the factor space
  evenly (scrambled, scipy.stats.qmc, design seed 20261006).
- V100-h: GPU hours on one V100 (cost unit; ≈ 0.1 $ per V100-h on the
  4×V100 instance at 0.41 $/h).

**Design.** 96 Sobol points (seed 1001 + point index, mapping fixed, so the
design can be extended) + 8 centre replicates C0 … C7 (d/l_ex 252, L_eff 18,
r_p 1.5, φ 0.62; seeds 1901 … 1908) + 4 corners × 3 seeds (seeds 2011 …
2013, 2021 … 2023, 2031 … 2033, 2041 … 2043), all at α 0.1; at α 0.01 only the
8 centre replicates (same seeds, started from the init.pt of the α 0.1
centre run, 9 cycles per amplitude): 124 runs. The α 0.01 Sobol subset is
dropped (Chris 2026-10-06); a later own α study can start from the kept
init.pt files. The corners check whether the realisation scatter is the
same over the factor space (the centre replicates measure it only at the
centre):

| Corner | d/l_ex | L_eff/l_ex | r_p | φ | why |
|---|---|---|---|---|---|
| K1 | 212 | 30 | 1.5 | 0.62 | fewest cubes per particle (nominal 185, realised 195; the smallest Sobol point has 183) |
| K2 | 300 | 12 | 3 | 0.62 | narrowest local wall (2 cells) |
| K3 | 212 | 18 | 1.5 | 0.69 | smallest gap between particles (1.7 cells) |
| K4 | 300 | 18 | 1.5 | 0.55 | largest box, many cubes |

**Fixed settings.** Js = 1.5 T (K_d = 895 kJ/m³), dx = 3 l_ex, f = 30 MHz,
B̂ = 9 / 20 / 35 / 50 / 70 / 100 / 150 mT, 7 cycles per amplitude at α 0.1
(9 at α 0.01), cycles 0 and 1 discarded, converged virgin state (hard
gate), uniaxial anisotropy, staircase surface, T = 0, no eddy currents.

**Two readings of the same runs.** The simulation uses reduced units: A
enters only through l_ex = √(2A/(μ₀Ms²)); at fixed Js the time scale and
K_d do not depend on A. Thus a run with d/l_ex = 212 is the same
calculation for "d = 1 µm, A = 19.9 pJ/m" and for "A = 10 pJ/m, d = 709 nm";
w (J/m³), f and β are the same. The effects of d and A cannot be
separated: it is one factor. **Reading 1 (main, d = 1 µm fixed, A varied)**
and reading 2 (A = 10 pJ/m fixed, d varied) below. Nominal ranges; the
realised values are in results/sobol_design.csv (d/l_ex 209.7 … 303.9,
thus A 9.7 … 20.4 pJ/m in reading 1).

**Factors and their real values** (Js = 1.5 T):

| Factor (reduced), sampling | Range | Reading 1: d = 1 µm, A varied | Reading 2: A = 10 pJ/m, d varied |
|---|---|---|---|
| d/l_ex, log-uniform | 212 … 300 | A = 19.9 … 9.9 pJ/m (l_ex 4.72 … 3.33 nm) | d = 709 … 1003 nm (l_ex 3.34 nm) |
| L_eff/l_ex, 4 levels (level = ⌊4u⌋, u = Sobol coordinate in [0, 1)), linear in L | 12 / 18 / 24 / 30 | Herzer length = cube edge 40 … 57 / 60 … 85 / 80 … 113 / 100 … 142 nm (at A 10 … 20 pJ/m) | 40 / 60 / 80 / 100 nm |
| ↳ Q_eff = (l_ex/L_eff)², K_eff = Q_eff K_d = A/L_eff² | 6.9 / 3.1 / 1.7 / 1.1·10⁻³ | K_eff 6.22 / 2.76 / 1.55 / 0.99 kJ/m³ (independent of A) | same |
| ↳ K1 for grains D = 10 nm (Herzer K_eff = K1⁴D⁶/A³, K1 ∝ A^¾) | | 50 … 84 / 41 … 69 / 35 … 59 / 32 … 53 kJ/m³ (at A 10 … 20 pJ/m) | 50 / 41 / 35 / 32 kJ/m³ |
| r_p, uniform (linear) | 0 … 3 | K_p = r_p K_eff: 0 … 18.7 kJ/m³ (L 12) … 0 … 3.0 kJ/m³ (L 30); residual stress, one uniaxial axis per particle | same |
| φ, uniform (linear) | 0.55 … 0.69 | packing fraction 55 … 69 % | same |

Regression variables: ln Q_eff, r_p, φ, ln d/l_ex. Note: the factor d/l_ex
alone is "A changed at fixed Q_eff and fixed r_p", not "A changed, material
otherwise the same"; the full A effect at 1 µm follows from the regression
(4.3: dβ/d ln A = −3 ∂β/∂ ln Q_eff − ½ ∂β/∂ ln(d/l_ex), plus a r_p term if
K_p is fixed). The grain size D = L_eff with K1 = K_eff (x = 1, PROTOKOLL 2.1)
is an equivalent reading of the cubes.

**Derived quantities and limits** (nominal; reading 1 / reading 2):

| Quantity | Min | Max | Limit / rule |
|---|---|---|---|
| mesh dx = 3 l_ex | 10.0 nm (A 10) / 10.0 nm | 14.2 nm (A 20) / 10.0 nm | fixed in reduced units |
| box edge a | 1.47 µm / 1.04 µm (d/l_ex 212, φ 0.69) | 1.56 µm / 1.56 µm (d/l_ex 300, φ 0.55) | a multiple of 12 l_ex |
| cells | 1.1 M | 3.8 M | cost ∝ cells |
| smallest gap between particles | 1.7 cells = 24 nm / 17 nm (d/l_ex 212, φ 0.69) | 10.4 cells = 104 nm (d/l_ex 300, φ 0.55) | ≥ 1 void cell over 26 neighbours (contact rule) |
| d/L_eff | 7.1 (d/l_ex 212, L 30) | 25 (d/l_ex 300, L 12) | ≥ 7 (decision D3, §7) |
| cubes per particle | 185 | 8 200 | residual check vs d/L_eff (decision D3, §7) |
| local wall width L_eff/√(1 + r_p) | 2.0 cells = 20 nm (L 12, r_p 3, A 10) … 28 nm (A 20) / 20 nm | 10 cells = 100 … 142 nm (L 30, r_p 0) / 100 nm | ≥ 2 cells (resolved) |
| wall width π L_eff | 126 nm | 445 nm / 315 nm | |
| drive H amplitude (μ ≈ 8, relative permeability of the core) | 0.9 kA/m (9 mT) | 14.9 kA/m (150 mT) | |
| α | 0.01 (8 centre runs) | 0.1 (main) | model parameter (PROTOKOLL 7.17) |

**Validity of the mesh** (PROTOKOLL 7.28, 7.32): no mesh effect on β
detected at the base point (1 µm, n = 4; 95 % CI up to ≈ ±0.2 … 0.4); the
small end d/l_ex 212 and the corners are not checked. The absolute loss at
dx 3 was higher than at dx 1.5 in 3 of 4 seeds (mean Δ ln w +0.19 … +0.36,
i.e. +21 … +43 %; not significant, n = 4); the direction to the converged
value is not known.
**Mesh check (rule fixed before the start, Chris 2026-10-06, audit 3
item 6, option A):** no further mesh tests before the design. After the
main run and its regression: **trigger** = in the main fit (Gamma GLM), a
factor has |range effect| ≥ 0.15 on β_core (any window or all) or on
ln E[w] (10 / 50 / 100 mT), with its 95 % CI (not corrected, on purpose
generous) excluding 0. L_eff or r_p (they set the wall width in cells) →
corner K2 (wall 2 cells); d/l_ex or φ (they set the gap in cells) → corner
K3 (gap 1.7 cells). **Check** = 6 seeds at dx 1.5 and dx 3, paired, mesh300
protocol (9 / 50 / 100 mT, 5 cycles), run without the 12 h job timeout
(`--job_timeout_h 30`: a dx 1.5 job needs ≈ 21 h); compared: the mesh
difference at the corner against the one at the base point (mesh300).
Equal → the effect is physical; different → the difference is the artefact
share. Expected 95 % CI of the difference of the mesh differences: ±0.13 … 0.21 (in β and in ln w). **Cost (estimates from mesh300
rates, cost ∝ cells):** K2 ≈ 250 V100-h ≈ 25 $ (≈ 41 V100-h per seed:
dx 1.5 and dx 3, 3 jobs each; ≈ 3.5 days on 4 GPUs), K3 ≈ 80 V100-h ≈ 8 $
(≈ 13 V100-h per seed, ≈ 1 day). Run only with a
separate go after the regression. Without trigger: no check; every effect
of L_eff, r_p, d/l_ex, φ is reported as a relative effect at dx 3 (a mesh
error of β up to ≈ 0.3 not excluded, in both directions); every d/l_ex effect on
ln E[w] carries the measured mesh caveat (Δ ln w at 50 mT: +0.26 ± 0.16 at
d/l_ex 300, −0.06 ± 0.14 at 212, 7.28).

**Outputs per run** (`analyze.features`): β in the 3 windows and β all,
ln w at 10 / 50 / 100 mT, and the flags per amplitude (definitions in the
`analyze.py` docstring): `c` loop not closed (drift), `o` loop centre
shifted (|offset| > 0.10), `*` energy balance |w_dis/w_loop − 1| > 2 %,
`#` jump of n60 (number of neighbour-cell pairs with > 60° between their m),
`!` cycle-to-cycle scatter above the target.
**Snapshots (decided 2026-10-06, Chris, option c):** all 124 runs: m at 16
phases in each of the last 2 cycles at 9, 50 and 150 mT (`--snap_phases 16
--snap_cycles 2`; 96 files of 48 MB at 148³ cells, measured; 1.6 … 5.5 GB
per run, mean 3.1 GB). Analysed on the instance right after each run:
`phasemap_batch.py`, `phase_extra.py` (R1, lag spread, where the loss sits,
texture, F90, slices) and `domains_batch.py` (wall measure, switched
volume, l_C of the last cycle at 90/180/270/360°, of init.pt and of the end
state checkpoint.pt; ≈ 14 states × 38 s CPU per run).
**Accumulators (decided Chris 2026-10-06, audit 3 item 9b, option B):** all
124 runs, all 7 amplitudes (`--acc`): per cell over the cycles ≥ 2 the mean
p, the mean m and the fundamental m₁ from all 256 samples per cycle (no
aliasing). Reduced on the instance by `acc_reduce.py` → acc.json (sum check
p_box/p_dis = box mean of p over the measured p_dis, must be 1; R1; shares
of p and p_fund and the lag in the classes core r/R < 0.5 / middle /
shell r/R ≥ 0.8 / texture top 10 % / p top 10 %) and acc_slices.npz.
Instance test (d/l_ex 96, 2 amplitudes): sum check 1.000000, R1 0.960 /
0.206 against 0.956 / 0.207 from the 16-phase snapshots, run time +4 … 5 %,
loss unchanged.

**Deletion and fetch rules.**
- Snapshots (*.vti) and accumulator files (acc_*.npz) (standing go, Chris
  2026-10-06): deleted on the instance by `instance_postproc.py` as soon as
  the analysis files are written and checked there and their md5 is logged;
  the results are fetched and md5-checked on the laptop later. If the
  post-processing of a run fails (POSTPROC_FAILED), its files stay; they are
  removed only by a manual step after the failure is understood.
- init.pt: stays on the instance until the whole design is finished (the
  α 0.01 centre runs and a later α study start from it); then deletion only
  with a separate go after the dependency check of `cleanup_instance.py
  --include_init`.
- checkpoint.pt (decision Chris 2026-10-06): deleted only by Claude during
  the daily fetch, after a dry run of `cleanup_instance.py`; conditions: the
  run is DONE and POSTPROC_DONE, domains.json holds the end state ("final"),
  and summary.json, domains.json, checkpoint.pt and init.pt are listed in the
  local MANIFEST.md5 (fetched and verified).
- Fetch (decided Chris 2026-10-06, audit 3 item 9a): during the daily fetch
  the result files and the new init.pt and checkpoint.pt are downloaded
  (state files of all runs: init.pt 2.98 GB of the 116 α 0.1 runs +
  checkpoint.pt 3.18 GB of all 124 runs = 6.15 GB, ≈ 2 GB per day),
  md5-verified and kept local as raw data (MANIFEST.md5, not in git).
- Disk (40 GB): used now 4.4 GB; kept per run init.pt + checkpoint.pt + CSV
  ≈ 0.05 GB (12 B per cell per state file; α 0.01 runs write no init.pt);
  expected bytes per run (snapshots + accumulators + state files)
  2.0 … 6.6 GB, mean 3.7 GB. **Disk guard with reservation** (`run_queue.py
  --min_free_gb 3`): a job starts only if the free space minus the expected
  bytes still to come of all running jobs minus its own expected bytes is
  ≥ 3 GB; the queue waits instead of filling the disk (tested on the
  instance). If checkpoints are not deleted for a day, fewer jobs run in
  parallel near the end.

**Unattended operation** (PROTOKOLL 7.33 items 2, 4, 5, C2; all tested).
The chain (`instance/chain_sobol.sh`) runs the code of one pinned commit,
restarts itself after an instance restart (onstart, lock), and appends all
logs. Order: K1_0, K2_0, K3_0, K4_0 first (relaxation, time, disk and
post-processing at the extremes in the first ≈ 4 h), then the other main
and corner runs in random order (seed 20261006), then the 8 α 0.01 runs.
`run_queue.py`: a failed job runs once more (it resumes from its
checkpoint, also after a gate retry); a relaxation-gate failure (exit code 3, deterministic) runs
once more with relax_maxiter 200000 and a fresh start, then it is excluded;
a NaN state stops a run (exit code 4); job timeout 12 h; failures of 2
different jobs on one GPU within 10 min pause it 30 min, the 3rd such burst
disables it (ALERT); the retry state of each job is kept in
QUEUE_STATE.json, so a restart gives no new attempts; status file
status_sobol.json every 5 min. Post-processing runs under a restart loop.

**Analysis (estimand decided Chris 2026-10-06, audit 3 item 7, option C).**
- Runs in the fits: the 104 main runs and the 12 corner runs at α 0.1;
  excluded runs (4.3a) are left out.
- Main fit: for each of the 7 amplitudes, w of each run is regressed on
  (ln Q_eff, r_p, φ, ln d/l_ex) with a Gamma GLM, linear + quadratic +
  2-factor interactions (15 terms); ln E[w] at 10 / 50 / 100 mT (10 mT:
  linear interpolation in ln B̂ between the predicted means at 9 and 20 mT,
  as in `analyze.features`) and β_core
  of the windows and over all amplitudes follow from the predicted means; the
  SE and the range effects come from a bootstrap over the runs (the runs are
  the independent samples).
- Cross-checks, reported next to the main fit: (A) OLS of ln w and of the
  per-run β with HC3 SE (it estimates E[ln w], the typical run); (B) Huber
  fit of the same (the bulk without rare states); the main fit without
  flagged runs and the main fit with w_dis instead of w_loop. **Flagged
  run** (rule fixed 2026-10-07 from the earlier data, PROTOKOLL 7.39):
  |w_dis/w_loop − 1| > 0.20 at any amplitude, i.e. a state change during the
  measured cycles (7 of 20 runs of the stopped design, 3 of 8 base-point
  runs); the flags `c`, `o`, `#`, `!` are reported per run only. A clear
  difference between the main fit and fit B
  means that the rate of rare states depends on the factors: reported as a
  finding.
- φ: w is per core volume, so the φ effect on ln E[w] contains the trivial
  dilution ln(0.69/0.55) = 0.23; it is reported with and without it.
- Scatter model (decided Chris 2026-10-06, replaces the F-test with
  n = 3): for the 2 primary outputs, the log of the squared residuals of
  fit A is regressed linearly on the 4 factors (all 116 runs, centre and
  corner replicates included); if a slope is significant (Holm over the 4
  slopes of one output, 5 %), the fits are repeated with the weights
  1/σ²(x) from this model. The corner and centre replicates give the pure
  scatter at 5 points as a check (only a change of σ by ≈ 2× is detectable).
- α 0.01: per seed the ratio to the α 0.1 run of the same seed at the
  centre (n = 8, paired): w and R = |m₁|²(α 0.01)/|m₁|²(α 0.1), with |m₁|²
  the mean over the magnetic cells from the accumulators; plus the loss map
  (acc.json) and the snapshots.

**Decision rules (fixed before the start, Chris 2026-10-06, audit 3 item 8).**
- Primary outputs: β_core over all amplitudes and ln E[w] at 50 mT, each
  against the 4 factors: 8 tests of the range effect, Holm correction,
  family level 5 %.
- Class (only for the 2 primary outputs): "lever" = range effect
  significant after the correction; "not a lever" = 90 % CI of the range
  effect inside ±0.15 (β; ln E[w]: ±0.15 = −14 … +16 % loss; TOST);
  otherwise "not determined". All other outputs (β windows, ln E[w] at 10 /
  100 mT, snapshot and accumulator measures), interaction and quadratic
  terms: reported with CI and called a hint, not a finding.
- α 0.01 centre test (predictions): viscous loss (∝ α) R ≈ 1; α-independent
  loss (micro-avalanches) R ≈ 10; resonance R ≈ 100. The hypothesis whose
  value lies in the 95 % CI of ln R is assigned; none → "mixed". Note: the
  paired loss ratio 0.59 at 50 mT (7.31) gives R ≈ 6 if the loss is mainly
  in the fundamental (estimate), so "mixed" is the likely outcome; the test
  then gives the share.
- Limits of the design (stated in the report): α and f are fixed and cannot
  be ranked as levers; r_p ≤ 3 is probably in the linear range (open option
  for later: pre-test r_p 3 and 6, PROTOKOLL §12 option B, ≈ 2 $); β at
  10 mT itself is not delivered (lowest window centre ≈ 18 mT).
- Expected SE of the range effects on β at N = 96 (σ from the baseline,
  for OLS of the per-run β, fit A): 0.06 (β all), 0.10 (9-35 mT), 0.06
  (35-70 mT), 0.07 (70-150 mT); the main fit is checked against these with
  synthetic data (below).
- Analysis code (audit 3 C4): the code for the main fit, the bootstrap,
  Holm, TOST and the scatter model is written and tested on synthetic data
  with a known answer during the runs, before any result is looked at.

**Restart after the stop (PROTOKOLL 7.36 … 7.38).** The first start
(2026-10-06) used a non-isotropic stress-axis set and was stopped; its 20
finished runs are archived (runs/_archive_sobol_v2_oldaxes, not used in the
analysis). The code since bc6d5bb has the isotropic tetrahedral axes
(7.37). The job files are unchanged (the axes are set in the code). Steps:
(1) Chris's go; (2) on the instance: no Sobol run folder may exist (else
run_queue skips the old DONE runs); fetch the pinned commit, write
/root/sobol_commit, copy the chain, add the onstart line, start (see the
head of instance/chain_sobol.sh); (3) after ≈ 4 h check K1_0 … K4_0.
Decided before the restart (Chris 2026-10-07): r_p stays relative
(K_p = r_p K_eff, r_p 0 … 3; the ratio of coherent to random anisotropy is
the transferable quantity; r_p ≫ 3 is a stated limit, pre-test option);
flag rule as above (7.39).

**Cost and time (estimates from the measured 3.4 V100-h per run at 3.2 M
cells, cost ∝ cells, + 4 … 5 % for the accumulators):** per run 1.2 … 4.2
V100-h (largest box 4.2), mean ≈ 2.4. Main design 104 runs ≈ 250 V100-h; corners 12 runs
≈ 30 V100-h; α 0.01 centre 8 runs ≈ 20 V100-h (per cycle the same rate as
α 0.1, measured); total ≈ 300 V100-h. Wall time ≈ 3.2 days on 4 GPUs. The
instance is paid per hour (0.41 $/h), also for idle GPUs at the end of each
queue: rental ≈ 32 … 35 $ (estimate). Credit now ≈ 11 $: top up to ≈ 45 $
(an empty credit stops the instance). The mesh check (if triggered) needs
its own budget: K3 ≈ 8 $, K2 ≈ 25 $ (estimates).

### 4.1 Factors and ranges

(4.1 … 4.4 are the earlier plan; where they differ from 4.0, 4.0 holds.)

| Factor | Range | Meaning |
|---|---|---|
| L_eff/l_ex | 12, 18, 24, 30 (discrete: the cube edge must be a multiple of 6 l_ex) | Q_eff = (l_ex/L_eff)² = 6.9·10⁻³ … 1.1·10⁻³ (Herzer anisotropy) |
| r_p | 0 … 3 | particle-scale stress anisotropy K_p/K_eff |
| φ | 0.55 … 0.69 | packing fraction |
| d/l_ex | 212 … 300 | at 1 µm and Js = 1.5 T: A = 20 … 10 pJ/m (exchange); together with Q_eff (∝ A⁻³ by Herzer) this gives the effect of A |

All corners are valid (checked 2026-10-02 and by the review): no particle
contact (smallest gap 1.7 cells at d/l_ex = 212, φ = 0.69), local wall width
≥ 2 cells (L_eff = 12, r_p = 3, exactly 2), d/L_eff ≥ 7.

Note on PROTOKOLL: PROTOKOLL §3 fixed d/l_ex = 300 and the window
d/L_eff ≥ 10 (≥ 500 cubes per particle). This plan varies d/l_ex (exchange
A, Chris 2026-10-02) and has d/L_eff = 7 … 25: at d/l_ex = 212 and
L_eff = 30 there are ≈ 185 cubes per particle (decision D3).

### 4.2 Design

- **On hold (Chris 2026-10-04)** until mesh300, a frequency test f/2, f, 2f at
  9 … 20 mT with the converged state, and phase-resolved snapshots plus a
  ring-down test are done (PROTOKOLL 7.24).
- **Decided 2026-10-04 (Chris, option 3):** main design at α 0.1 with
  N = 96 + 8; a subset at α 0.01 with the first 32 Sobol points + 8 centre
  replicates, the same points and seeds; extensible to the full 96 + 8 at
  α 0.01 if needed (decision after the evaluation of the subset).
  **Superseded 2026-10-06 (Chris):** at α 0.01 only the 8 centre
  replicates, with snapshots (PLAN 4.0); a later own α study is possible.
- N points of a scrambled Sobol sequence in 4 dimensions (scipy.stats.qmc,
  fixed design seed). The design must stay extensible (a prefix of the
  sequence is a valid design): the mapping of each point does not depend
  on N. L_eff: level = ⌊4u⌋ of the Sobol coordinate u (fixed; for N = 32,
  64, 96 the Sobol net gives exactly N/4 points per level); r_p and φ
  linear in u, d/l_ex linear in ln. The regression uses the realised values
  (`d_lex_used`, `phi_vox`), not the nominal ones.
- One realisation per point (own seed: own cubes and own initial state).
  Seed = 1001 + point index, the same at α 0.1 and 0.01 (no overlap with
  the pilot seeds 1 … 8 and the test seeds 901 … 928).
- α 0.01 runs start from the relaxed virgin state of the α 0.1 run of the
  same point (`--init_from`; the relaxation runs at α = 1 and does not
  depend on α) and use 9 cycles per amplitude instead of 7 (at α 0.01 the
  9 mT loops jump from cycle to cycle, PROTOKOLL 7.15).
- Random order of the runs over the GPUs (no correlation of factor and GPU).
- 8 replicates at the centre point (L_eff = 18, r_p = 1.5, φ = 0.62,
  d/l_ex = 252): pure error σ independent of the model (with the corner
  replicates a check of the scatter model, 4.0).
- Disk, deletion and fetch: see PLAN 4.0 (superseded here).
- Extension: the next Sobol points of the same sequence (up to 128 / 256 in
  total at α 0.1) if the errors of step 2 are too large (decision by
  Chris); at α 0.01: a later own α study (separate plan).

### 4.3 Analysis

- Per run (`analyze.features`, changed 2026-10-04): local Steinmetz exponent
  from a power law w = c B^β fitted over 3 neighbouring amplitudes
  (9-20-35, 35-50-70, 70-100-150 mT; centres ≈ 18 / 50 / 102 mT), β over
  all amplitudes (β_pow), and ln P at 10 / 50 / 100 mT (interpolated in
  ln b). 3 points leave 1 degree of freedom: χ² checks the local power law.
  No quadratic fit in ln-ln (it forces a shape; the curvature was not
  significant, −0.07 ± 0.09).
- Regression of each output on x = (ln Q_eff, r_p, φ, ln d/l_ex): linear,
  quadratic and two-factor interaction terms (15 coefficients), ordinary
  least squares (superseded by 4.0: Gamma GLM on w as main fit, OLS and Huber
  as cross-checks, scatter model instead of the F-test).
  The centre is not exactly at the coded zero for ln Q_eff (+0.11); the
  regression does not need it there.
- Ranking: change of the output over the full range of each factor (with
  SE), main effects and the largest interactions.
- Fixed check (D3): residuals vs d/L_eff, separately for d/L_eff < 10 and
  ≥ 10 (points with fewer than ≈ 500 cubes per particle).
- Effect of A at fixed 1 µm and fixed r_p: dβ/d ln A = −3 ∂β/∂ ln Q_eff
  − ½ ∂β/∂ ln(d/l_ex). If instead the stress anisotropy K_p is fixed, then
  r_p = K_p/K_eff ∝ A³ and the term +3 r_p ∂β/∂r_p is added. A factor 2 in A
  changes Q_eff by 8×, more than the design range (6.25×): the formula is a
  local derivative (review estimate: SE ≈ 0.27 per unit ln A at N = 128).
- Target error: SE ≤ 0.1 for the range effect of each factor on β.
  **Update 2026-10-04 (PROTOKOLL 7.17, baseline n = 8, converged virgin
  state):** σ = 0.15 (β_pow), 0.27 (9-35 mT), 0.15 (35-70 mT), 0.19
  (70-150 mT). With the scaling below: N = 29 / 97 / 30 / 47 for SE ≤ 0.1;
  at N = 64: SE 0.07 / 0.12 / 0.07 / 0.09; at N = 96: SE 0.06 / 0.10 / 0.06 /
  0.07. Decision D1 again by Chris. The text below is the original basis.
  Basis: σ_β ≈ 0.45 per realisation (pilot, β 9 → 50 mT, dx 1.5; the value
  0.11 at dx 3 comes from n = 3 only). Review Monte Carlo with this design:
  N = 128 + 8 gives SE ≈ 0.11 (Q_eff) … 0.14 (r_p, φ, d/l_ex); about
  256 + 8 runs are necessary for SE ≤ 0.1. β at 10 mT comes from 9 and
  20 mT only (Δ ln b = 0.8) and has a larger σ (≈ 0.4 … 1.0).

### 4.3a Exclusion and failure rules (fixed before the runs)

- A run is excluded if it crashed, gives NaN, or misses an amplitude.
  A failed point is run again once with the same seed (it resumes from its
  checkpoint); if it fails again, it is reported and left out (no
  replacement point, to keep the Sobol balance). **Relaxation gate not met
  (exit code 3, decided 2026-10-06, audit 3 item 5):** the relaxation is
  deterministic (same seed, same iterations), so the one retry uses the
  same seed with relax_maxiter 200000 (instead of 80000) and starts again
  (--no_resume); a converged state does not depend on the limit. If it
  fails again, the point is excluded; the factor values of all excluded
  points are reported, to show a cluster in the factor space.
- Order: one seed of each corner K1 … K4 first (tests relaxation, time,
  disk and post-processing at the extremes in the first ≈ 4 h), then the
  shuffled rest.
- Flags `#`, `o`, `c`, `*`, `!` are reported per run but do not exclude it.
  A sensitivity analysis repeats the regression without the flagged runs
  (flagged: |w_dis/w_loop − 1| > 0.20 at any amplitude, PLAN 4.0) and with
  w_dis instead of w_loop.
- Raw data stay local with MANIFEST.md5 (no commit); the per-run feature
  CSV and the regression results are committed.
- Results go to PROTOKOLL (tables, response plots, recommendation with
  validity limits).

### 4.4 Cost and time (estimates)

| | Value |
|---|---|
| per run, 7 amplitudes, converged relaxation (measured, d/l_ex 300) | 3.4 V100-h (α 0.1); α 0.01 per cycle the same rate (measured, audit 3) |
| main design 96 + 8 + corners 12 at α 0.1 | ≈ 280 V100-h (incl. accumulators, PLAN 4.0) |
| α 0.01: 8 centre runs (decided 2026-10-06) | ≈ 20 V100-h |
| rental (0.41 $/h, incl. idle GPUs) | ≈ 32 … 35 $ (PLAN 4.0) |
| a later own α study | separate plan and go |

GPU selection: offers are shown to Chris before rent (low hourly rate).

### 4.5 Completion to 2 calm cycles, then the conditioning test (2026-10-11, for approval)

Reason: PROTOKOLL 7.53 … 7.55. The events in the measured cycles are the end
of a transient. The steady analysis (7.54, main result by decision of
Chris) keeps only the measured cycles after the last event. Chris wants
every run-amplitude with at least 2 such cycles, then the conditioning test
(option 1 of 7.55). Reviewed by a context-free agent (2026-10-11); its
findings are worked in.

Terms. A **stage** is one amplitude of a run. The **measured cycles** of a
stage are the cycles from 2 on (cycle 0: drive ramp, cycle 1: drive
correction). A **calm cycle** has w_loop > 0 and |w_dis / w_loop − 1| ≤ 0.20.
The **calm tail** of a stage is the number of calm cycles at the end of its
measured cycles; it equals the number of cycles that the steady analysis
(7.54) keeps. "Calm" is not the old flag `steady` of run_loops (closure and
dW rule), which this plan does not use. Cost unit: V100-h (1 GPU on the
4 × V100 instance).

Is 2 enough? Probability of an event in the next cycle after k calm cycles
in a row (116 runs, cycles 1 … 6; eval_events.py, numbers in PROTOKOLL
7.55 / chat 2026-10-11): k = 1 / 2 / 3 / 4: 1.6 / 0.8 / 0.4 / 0.1 % (20 … 35
mT: 2.6 / 1.7 / 0.7 / 0.3 %). For the value of the loss density, 2 cycles
are enough: the scatter from cycle to cycle is 2.5 %, the scatter between
realisations 16 … 35 %. A rule with an extra closure limit (10⁻³) would need
a rerun of 105 of 116 runs; not proposed. The conditioning test (D)
measures the remaining drift directly.

A. Code (no GPU).
   - `run_loops.py`: `--min_calm N` (default 0: off, the old rule) and
     `--calm_tol 0.20`. With N > 0 a measured stage runs at least
     `cycles_per_amp` cycles and then continues until its calm tail is ≥ N,
     at most `max_cycles_per_amp` cycles. With N > 0 this replaces the
     closure / dW stop rule. `--cycles_last_amp M`: the last stage runs M
     cycles (for D). Summary: calm tail per stage; a stage that reaches the
     maximum without N calm cycles is marked.
   - Analysis: rerun folders have the suffix `_mc2` (min calm 2). Option
     `--reruns` of analyze_design (with --steady_tail): a `_mc2` run replaces
     its original. The pre-registered analysis stays on the originals.
     PATTERN and the post-processing pattern accept the suffix. eval_k3 and
     eval_jumps_time get the same mapping. eval_k3 also gets the steady-tail
     rule (the verdict of 7.52 used all measured cycles, events included).
   - Determinism gate `check_reruns.py`: up to the first extended stage, a
     `_mc2` run must repeat its original cycle by cycle (w_loop and w_dis
     equal to 6 significant digits). Same GPU model (Tesla V100-SXM2-16GB)
     and torch (2.7.1+cu126) as the originals: true, same instance.
   - Tests: unit test of the stop decision; a small CPU run with N = 2 gives
     the same cycles 0 … 6 as without the option; mapping test.
B. Reruns. Arguments as in the original job line, with `--init_from
   <run>/init.pt` (always explicit), `--min_calm 2`, `--max_cycles_per_amp`
   replaced (15 for the Sobol runs, 12 for the mesh jobs), all `--snap_*`
   removed, accumulators on. All init.pt files are on the instance (checked
   2026-10-11).
   - 13 Sobol runs with a calm tail < 2 (7.54): S033, S035, S055 (first
     extended stage 9 mT), K4_1, S076, S091 (20 mT), K4_0, S020, S026, S047,
     S067 (35 mT), S034 (70 mT), K1_0 (100 mT). They are not a random
     sample: 11 of 13 have A ≤ 15 pJ/m, 6 of 13 have L_eff 30.
   - Consequence: an extended stage changes the start state and the drive
     correction of all later stages of that run. 56 of the 91 stages of the
     13 runs change, although they were calm before. After an event the
     loss density is lower (7.55), so the change has a direction. Thus the
     steady fit is reported with the originals and with the `_mc2` runs
     (sensitivity row).
   - The other 103 runs have a calm tail ≥ 2 at all amplitudes, so the rule
     would stop them at 7 cycles, as run (true by construction). The rerun
     is the continuation of the original only if the GPU run is
     deterministic: the gate in A checks this. On a mismatch the chain
     starts no new job and reports.
   - 4 mesh / K3 jobs with a calm tail < 2: M300_dx15_s922_b9,
     M300_dx3_s921_b50, M300_dx3_s921_b100, MK3_dx15_s2034_b50.
   - K3 fill-in jobs (running): the list of jobs with a calm tail < 2 is
     made after the fill-in ends. Cap: 15 V100-h.
   - A stage that reaches the maximum without 2 calm cycles is missing in
     the steady analysis (the original value is not used) and is listed.
   - Limits that stay: the accumulators of a `_mc2` run sum all cycles from
     2 on, more cycles than in the other runs, events included. The
     snapshot measures of the 13 runs stay from the originals; at 9 mT the
     snapshots of S033, S035 and S055 lie in event cycles. Both are marked in
     the mechanism evaluations.
C. Analysis after B: acceptance = every stage of the steady analysis has a
   calm tail ≥ 2; every exception is listed. Results: steady analysis with
   `_mc2` (main), the same with the originals (sensitivity), pre-registered
   (check). K3 mesh check again with the steady-tail rule. The time-series
   measures keep their rule (≥ 3 kept cycles) and list the gaps.
D. Conditioning test (option 1 of 7.55), centre point only, seeds C0 … C5
   (n 6). History as in the design up to the tested stage, then 30 cycles:
   per seed 3 runs from init.pt, `--min_calm 2 --cycles_last_amp 30`:
   (a) 9 mT; (b) 9 → 20 mT; (c) 9 → 20 → 35 mT. Cycles are numbered 0 … 29.
   - Drift d = ln P (mean of cycles 28 … 29) − ln P (mean of the calm tail
     of cycles 2 … 6, the design estimator), per seed and amplitude.
   - Decision per amplitude, mean of d over the 6 seeds, t(5): the design
     value stands if the 90 % CI lies inside ±0.05 (5 %, 2 × the cycle
     scatter); systematic error if the 95 % CI excludes 0 and |mean| ≥ 0.05;
     else not determined. Also β 9 → 20 and 20 → 35 mT from the design
     estimator and from cycles 28 … 29.
   - Report d and the number of events per seed (rare late events make the
     distribution two-peaked; n 6 is enough for a slow continuous drift,
     not for rare events).
   - The runs (b) and (c) repeat the 9 / 20 mT stages of the design: an
     extra check of the determinism against C0 … C5.
E. Order and cost (estimates; measured cycle times; 0.41 $/h):
   - fill-in K3 (running, end ≈ 11 UTC 2026-10-11) → one queue: the long
     mesh job first, then the Sobol reruns, the other mesh jobs, the
     fill-in reruns, then D.
   - Sobol reruns: 29.6 V100-h without extension (measured), ≈ 35 with it.
     Mesh jobs: 16 … 28 V100-h (1.72 / 0.49 / 0.06 h per cycle, 7 … 12
     cycles; the dx 1.5 job at d/l_ex 300 alone 12 … 21 h). Fill-in reruns:
     ≤ 15. D: 6 seeds × 111 cycles × 0.038 h ≈ 25 V100-h.
   - Total ≈ 90 … 105 V100-h, wall ≈ 26 … 30 h, ≈ 11 … 13 $ plus idle time.
     The credit after the fill-in is ≈ 0.7 $: a top-up of ≈ 20 $ is needed
     before the start.

## 5. Optional step 3: mesh control (separate go)

**Superseded 2026-10-06 by the mesh-check rule of PLAN 4.0** (trigger,
corners K2 / K3, 6 seeds). The text below is the earlier idea.

One factor change (L_eff 12 → 30 at the base point) at dx = 3 and dx = 1.5:
is the effect the same on both meshes? A run at dx = 1.5 and d/l_ex = 300
costs 2.8 h relaxation + 1.87 h per cycle (measured). With 3 realisations
per arm the SE of the difference of the effects is ≈ 0.45: this cannot find
anything. A useful version: 2 amplitudes (9 and 50 mT), ≥ 8 realisations
per arm: ≈ 16 × 16 h ≈ 260 V100-h ≈ 35 $ (dx 1.5 part) + small dx 3 part.
Decision by Chris after step 2.

## 6. Code work before step 2 (no GPU cost)

1. `sobol_design.py`: Sobol design (validity check of every point, seed
   per point, extension by N) → jobs/sobol_*.txt (done).
2. `analyze_design.py` (not written yet; PLAN 4.0: written and tested on
   synthetic data with a known answer during the runs, before any result is
   looked at): per-run features to a CSV, Gamma GLM with bootstrap, fits A
   and B, scatter model, Holm, TOST, plots.
3. Tests: design generation (reproducible, valid, extensible: done),
   regression on synthetic data with a known answer (with item 2).

## 7. Open decisions for Chris

| | Question | Options (cost estimates) |
|---|---|---|
| D1 | number of design points | **decided: 96 + 8 at α 0.1 + 4 corners × 3 seeds, extensible; at α 0.01 only the 8 centre replicates (2026-10-06; the α 0.01 Sobol subset is dropped)** |
| D2 | amplitudes | **decided: 7 (9 / 20 / 35 / 50 / 70 / 100 / 150 mT)** |
| D3 | d/l_ex as factor and the cube window | **decided: (a) d/l_ex 212 … 300, d/L_eff ≥ 7; check of the residuals vs d/L_eff (below / above 10) in the analysis** |
| D4 | frequency test corner | **decided: base point + corner L_eff 30, r_p 3** |

## 8. Risks

| Risk | Consequence | Measure |
|---|---|---|
| large dynamic share at 30 MHz | β contains the damping loss | accepted: the target is the total loss with the model α (PROTOKOLL §12) |
| regression model too simple | biased effects | residuals vs factors, scatter model, centre and corner replicates; extend points |
| σ larger at some corners (small d, 9 mT) | larger errors there | scatter model (4.0) with weights 1/σ²(x); corner replicates; extension |
| level bias of dx = 3 (higher in 3 of 4 seeds, +21 … +43 %, not significant; direction unknown) | absolute losses uncertain; d effect on ln E[w] mesh-confounded | state as validity limit; mesh-check rule (4.0) |
| d/l_ex = 212 is not "1 µm" for A = 10 pJ/m | — | interpret d/l_ex as exchange at fixed 1 µm (Js fixed) |
