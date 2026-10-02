# PLAN: β(B_peak) production study (draft for approval)

Status: approved settings 2026-10-02 (D1 … D4, section 7). Step 1 is set
up; it starts after Chris approves the GPU selection. Step 2 needs go 2.
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

- N points of a scrambled Sobol sequence in 4 dimensions (scipy.stats.qmc,
  fixed design seed; N = 128 or 256, decision D1). Mapping: L_eff by
  quantiles to the 4 levels (exactly N/4 points per level), r_p and φ
  linear, d/l_ex linear in ln. The regression uses the realised values
  (`d_lex_used`, `phi_vox`), not the nominal ones.
- One realisation per point (own seed: own cubes and own initial state).
  Seeds start at 1001 (no overlap with the pilot seeds 1 … 8).
- Random order of the runs over the GPUs (no correlation of factor and GPU).
- 8 replicates at the centre point (L_eff = 18, r_p = 1.5, φ = 0.62,
  d/l_ex = 252): pure error σ independent of the model, lack-of-fit test.
- Extension: the next 128 Sobol points of the same sequence (256 in total)
  if the errors of step 2 are too large (decision by Chris).

### 4.3 Analysis

- Per run: ln P and β at 10 / 50 / 100 mT (local fit, PROTOKOLL §8).
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
| N = 128 + 8, 7 amplitudes (decided) | ≈ 270 V100-h ≈ 35 … 43 $, ≈ 34 h on 8 GPUs |
| extension to 256 + 8 | + ≈ 255 V100-h ≈ + 33 … 40 $ |
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
| D1 | number of design points | **decided: 128 + 8 first, then extension to 256 if the errors are too large** |
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
