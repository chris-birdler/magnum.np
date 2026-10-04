# PLAN: β(B_peak) production study (draft for approval)

Status 2026-10-03: step 1 (frequency test) FAILED: from ≈ 35 mT the loss
at 30 MHz is mainly damping loss (PROTOKOLL 7.9). Step 1b: rule not met
(PROTOKOLL 7.10). Step 2 is on hold until Chris decides the method for the
hysteresis loss (PROTOKOLL §12).
Basis: PROTOKOLL.md (model, tests, mesh decision).

## 1. Question

Which parameters control the Steinmetz exponent β of a powder core of 1 µm
particles at B_peak = 10 … 150 mT, and how can β stay small?

Output: β at 10, 50 and 100 mT and the loss density P (W/cm³ at f = 30 MHz)
as functions of four factors, with errors; a ranking of the factors.

## 2. Fixed settings (from PROTOKOLL)

| Setting | Value | Reason |
|---|---|---|
| model | Herzer cubes (Q_eff, L_eff) + particle-scale stress r_p; FCC cell of 4 spheres | PROTOKOLL §2 |
| mesh | dx = 3 l_ex (10 nm) | β not mesh-dependent at 1 µm; level known to ≈ 30 % (PROTOKOLL §4) |
| initial state | virgin state (random m on 6 l_ex blocks, relaxed), ascending amplitudes | PROTOKOLL §7.4 … 7.6 |
| damping | α = 0.1 | steady cycles (PROTOKOLL §7.2) |
| precision | fp32, atol 10⁻⁵ | numerical floor test PASS |
| amplitudes | B_peak = 9 / 20 / 35 / 50 / 70 / 100 / 150 mT (D2), 7 cycles each, cycles 0 and 1 discarded | study range; β at 10 / 50 / 100 mT by local fit (two-sided at 50 and 100 mT) |
| frequency | f = 30 MHz (f/f_M = 7.14·10⁻⁴) | unless step 1 shows a large dynamic share |

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

### 4.1 Factors and ranges

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

- **Decided 2026-10-04 (Chris, option 3):** main design at α 0.1 with
  N = 96 + 8; a subset at α 0.01 with the first 32 Sobol points + 8 centre
  replicates, the same points and seeds; extensible to the full 96 + 8 at
  α 0.01 if needed (decision after the evaluation of the subset).
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
  d/l_ex = 252): pure error σ independent of the model, lack-of-fit test.
- Disk of the instance (40 GB): large files are deleted when no longer needed
  (Chris, 2026-10-04), with `cleanup_instance.py`: dry run first, deletion
  only after the go for the printed list. Snapshots only after domains.json
  and dissmap.json are fetched and md5-verified; checkpoints only after the
  run is DONE and fetched; **init.pt never by default** (the α 0.01 runs
  start from the init.pt of the α 0.1 run of the same point).
- Extension: the next Sobol points of the same sequence (up to 128 / 256 in
  total at α 0.1; up to 96 at α 0.01) if the errors of step 2 are too large
  (decision by Chris).

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
  least squares (the realisation scatter dominates the cycle statistics);
  residual σ compared with the centre replicates (lack-of-fit F-test).
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
  A failed point is run again once with the same seed; if it fails again,
  it is reported and left out (no replacement point, to keep the Sobol
  balance).
- Flags `#`, `o`, `c`, `*` are reported per run but do not exclude it.
  A sensitivity analysis repeats the regression without the flagged runs.
- Raw data stay local with MANIFEST.md5 (no commit); the per-run feature
  CSV and the regression results are committed.
- Results go to PROTOKOLL (tables, response plots, recommendation with
  validity limits).

### 4.4 Cost and time (estimates)

| | Value |
|---|---|
| per run (V100, dx = 3), 5 amplitudes | 0.9 h (d/l_ex = 212) … 2.2 h (d/l_ex = 300); mean ≈ 1.42 h (review) |
| per run, 7 amplitudes, converged relaxation (measured, d/l_ex 300) | 3.4 V100-h (α 0.1); α 0.01 with 9 cycles ≈ 4 V100-h (estimate) |
| main design 96 + 8 at α 0.1 (decided 2026-10-04) | ≈ 240 … 350 V100-h ≈ 25 … 35 $ (estimate; d/l_ex 212 … 300 is cheaper than 300; exact from jobs.py when built) |
| subset 32 + 8 at α 0.01 (decided) | ≈ 110 … 160 V100-h ≈ 11 … 16 $ (estimate) |
| later extension of α 0.01 to 96 + 8 | + ≈ 64 runs ≈ + 20 … 26 $ |
| if step 1 needs f / 2 f pairs | × 1.5 |

GPU selection: offers are shown to Chris before rent (low hourly rate).

## 5. Optional step 3: mesh control (separate go)

One factor change (L_eff 12 → 30 at the base point) at dx = 3 and dx = 1.5:
is the effect the same on both meshes? A run at dx = 1.5 and d/l_ex = 300
costs 2.8 h relaxation + 1.87 h per cycle (measured). With 3 realisations
per arm the SE of the difference of the effects is ≈ 0.45: this cannot find
anything. A useful version: 2 amplitudes (9 and 50 mT), ≥ 8 realisations
per arm: ≈ 16 × 16 h ≈ 260 V100-h ≈ 35 $ (dx 1.5 part) + small dx 3 part.
Decision by Chris after step 2.

## 6. Code work before step 2 (no GPU cost)

1. `jobs.py`: group `freqtest` (step 1) and group `sobol` (design, validity
   check of every point, seed per point, extension by m).
2. `analyze_design.py`: per-run features to a CSV, regression, lack-of-fit,
   ranking, plots.
3. Tests: design generation (reproducible, valid, extensible), regression on
   synthetic data with a known answer.

## 7. Open decisions for Chris

| | Question | Options (cost estimates) |
|---|---|---|
| D1 | number of design points | **decided 2026-10-04: 96 + 8 at α 0.1 plus 32 + 8 at α 0.01 (same points), both extensible** (was: 128 + 8, based on σ_β ≈ 0.45 before the converged relaxation) |
| D2 | amplitudes | **decided: 7 (9 / 20 / 35 / 50 / 70 / 100 / 150 mT)** |
| D3 | d/l_ex as factor and the cube window | **decided: (a) d/l_ex 212 … 300, d/L_eff ≥ 7; check of the residuals vs d/L_eff (below / above 10) in the analysis** |
| D4 | frequency test corner | **decided: base point + corner L_eff 30, r_p 3** |

## 8. Risks

| Risk | Consequence | Measure |
|---|---|---|
| large dynamic share at 30 MHz | β measures damping | step 1 first; options for Chris |
| regression model too simple | biased effects | lack-of-fit test with centre replicates; extend points |
| σ larger at some corners (small d, 9 mT) | larger errors there | weights from the cycle statistics; extension |
| level bias of dx = 3 (≈ 30 %) | absolute losses uncertain | state as validity limit; optional step 3 |
| d/l_ex = 212 is not "1 µm" for A = 10 pJ/m | — | interpret d/l_ex as exchange at fixed 1 µm (Js fixed) |
