# Handover: β(B_peak) study with magnum.np (fp32)

Owner: Christoph Vogler (Chris). Branch: `study/beta-bpeak-pinning` (based on
`perf/memory-optimizations`, the fp32 fork). Folder: `studies/beta_bpeak/`.

Language rules for this study: documents in ASD-STE100 Simplified Technical
English, SI units only. Chris writes in German or English.

**Decision rule.** Chris makes the decisions. Before you change the physics
model, the parameter matrix, the material values or the budget, give him the
options with pros, cons, chances and risks. Then wait for his decision.
Engineering fixes (bugs, logging, speed) do not need a decision, but tell him.

---

## 1. Question

Steinmetz: P = k f^α B_peak^β. In soft magnetic powder cores (FeSiCr, particle
diameter 1 µm to 5 µm), how do the particle size, the fill factor φ and short-range
internal stress (dislocations after pressing) change β as a function of
B_peak? Practical goal: keep β small up to high B_peak without a large
increase of the total loss.

## 2. Physics background (result of the discussion before this handover)

1. **Shape of β(B_peak).** Small B: linear response, all losses ∝ B² → β ≈ 2.
   Medium B: Rayleigh hysteresis (W ∝ B³) grows faster → β rises towards 3.
   Large B: the loop becomes square, H_c is almost constant → β falls (to 1).
   So β_eff(B) = Σ (P_i/P) β_i is a weighted mean of the loss parts.
2. **Wall stiffness in 1 µm to 5 µm particles is magnetostatic.** A wall shift x
   in a particle gives a net moment. Stiffness k_ms ≈ 4 μ0 N Ms²/d
   (≈ 7·10¹¹ N/m³ at 5 µm, 3·10¹² N/m³ at 1 µm, N = 1/3). The pinning
   stiffness from K1 at grain boundaries (≈ 7 γ_w/D²) is 10² to 10³ times
   smaller. Consequence: A and K1 are almost unimportant for β.
3. **What pins against k_ms.** Only fluctuations of the wall energy on a
   short length scale ξ (10 nm to 100 nm). Curvature of the pinning potential
   ≈ Δγ_w/ξ² is (D/ξ)² ≈ 10³ to 10⁴ times larger than that of grain-boundary
   pinning. Dislocation stress fields do this through
   K_σ = (3/2) λs σ. A HOMOGENEOUS stress only adds a uniform anisotropy
   (like K1) and does not pin. Barkhausen criterion: a jump happens where
   the pinning curvature is more negative than −k_total.
4. **Expected effects.**
   - σ_rms·λs large → more Barkhausen jumps → β ≈ 2.5 to 3 over a wide B range.
   - Smaller d → larger k_ms → fewer jumps → β ≈ 2 up to higher B.
   - Fill factor φ: higher φ lowers the particle magnetization at given B
     (good, nucleation later) but also lowers N_eff and k_ms (bad, jumps
     earlier). There is an optimum φ that depends on the pinning strength.
   - Annealing helps only if it removes dislocations (recrystallisation),
     not only long-range residual stress (recovery).

### Hypotheses to test

| ID | Hypothesis | Signature in the results |
|---|---|---|
| H1 | Short-range stress (σ_rms, ξ) causes the β ≈ 3 plateau | β_max and plateau width increase with σ_rms; σ = 0 gives β ≈ 2 up to the nucleation range |
| H2 | k_ms ∝ 1/d delays the jump onset | onset B of β > 2 moves up when d goes down (0.5 → 0.75 → 1 µm) |
| H3 | φ trade-off | onset B of β > 2 is not monotonic in φ, or the direction changes with σ_rms |
| H4 | A and K1 are secondary | K1 = 20 kJ/m³ vs 30 kJ/m³ changes W and β less than the seed scatter |

## 3. Decisions already made (with reasons)

| Item | Decision | Reason |
|---|---|---|
| Model level | Full micromagnetics, FCC unit cell, 4 spheres, periodic | only this model gives k_ms, jump statistics and nucleation without assumptions |
| Demag | `DemagFieldPBC` (true periodic, k = 0 mode removed) | infinite core without macroscopic demag field, like a toroid: H_macro = H_ext |
| Cell size | 10.2569 nm for all production runs (N = 72/108/144/140/160) | cost ∝ Δx⁻⁵; wall parameter Δ = √(A/K1) ≈ 22 nm → 2.2 cells; ξ = 50 nm → 5 cells. Validate with 7.03 nm |
| Particle size | d = 0.5, 0.75, 1.0 µm (not 1 µm to 5 µm) | V100/5090 memory and time. Use the d-scaling to extrapolate |
| Damping | α = 0.1 | near critical damping of the wall oscillator → fastest relaxation; α = 1 is 10× slower |
| Frequency | 30 MHz, check with 15 MHz | quasi-static for the walls (relaxation ≈ 0.1 ns) and for rotation (≈ 1 ns ≪ 33 ns) |
| Protocol | 1 saturating cycle, then descending amplitudes (factor 0.5, 9 steps), ≥ 3 cycles each | the descending series is an AC demagnetization; wiping-out and congruency; check once with ascending order |
| Precision | fp32 (`set_precision("single")`), time in fp64 | 2× faster, validated in `bench/fp32_validation.pdf`; check fp32 vs fp64 at small B |
| Stress model | Gaussian random field, 6 components, σ_rms per component, correlation length ξ | 2 parameters, can be calibrated with XRD micro-strain |
| Drive | homogeneous, along x ([100] of the cube) | crystal axes are random per particle |
| Crystal axes | one single crystal per particle; 3 fixed orientation sets `AXES_SEEDS = {1: 16295, 2: 13903, 3: 6982}` (jobs.py); all main series use set 1 | the sets are selected so that the statistics along the drive are those of an isotropic powder (E_a/K1 = 0.20, nearest-easy-axis cos = 0.83, see note N9). A random set of 4 is too unreliable: the first choice (seed 1) was a 2.7σ outlier |
| No supercell | 1 FCC cell (4 particles), no 2×2×2 supercell | per cost, a supercell gives the same statistics as more cycles/seeds; it does not remove the systematic limitations (section 10). Decision Chris, 2026-09-29 |
| GPU | RTX 5090 (verified hosts on vast.ai) for fp32; fp64 runs only on V100/A100/H100 | best price per V100-hour; consumer GPUs have 1/64 FP64 |

## 4. Material values

**Status: assumptions. Replace them with measurements when Chris gives them.**

| Parameter | Value | Source / status |
|---|---|---|
| μ0 Ms | 1.80 T | assumption for Fe-3.5Si-4.5Cr (mass %). **Most important parameter (k_ms ∝ Ms²). Ask for the VSM value.** |
| A | 15 pJ/m | estimate: Fe D = 281 meV Å² → 19.8 pJ/m; Si and Cr lower D (Fe-4 at% Si 260 meV Å², Fe-30Cr 220 meV Å²) → 13 pJ/m to 17 pJ/m |
| K1 | 30 kJ/m³ (cubic) | Fe-Si torque data K1 = 52.9 − 2.79·A kJ/m³ (A in at% Si), scaled to room temperature; Cr effect not known → 25 kJ/m³ to 32 kJ/m³ |
| λs | 5·10⁻⁶ (isotropic) | assumption; Fe-3Si gives ≈ 6.5·10⁻⁶. Single crystals have λ100 ≠ λ111 (not modelled) |
| σ_rms | 0, 150, 400 MPa | 400 MPa ≈ yield stress (cold-worked), 150 MPa partially recovered. XRD: σ_rms ≈ E·ε_rms |
| ξ | 50 nm | dislocation spacing scale; must be ≥ 3 cells |
| Void | μ0 Ms = 10⁻⁸ T, A = 0, K1 = 0 | as in the fp32 validation run; Ms must be > 0 (division by Ms) |

## 5. Code map

| File | Content |
|---|---|
| `geometry.py` | FCC box with 4 voxelized spheres, packing fraction, contact check (no face contact of different particles → no exchange coupling), random rotations |
| `stress.py` | Gaussian random stress field (FFT filter exp(−k²ξ²/4) → covariance exp(−r²/2ξ²)), normalised to σ_rms on the magnetic cells |
| `field_terms_extra.py` | `StressAnisotropyField` (h = 3λs σ·m/(μ0 Ms)), `SinusoidalDrive` (H(t) from state.t, no device sync) |
| `run_loops.py` | one run: setup, relaxation, protocol, diagnostics, checkpoints, resume |
| `analyze.py` | W(B_peak), P, β_eff = d ln W/d ln B, error bars, plots, frequency separation (`--pair`) |
| `jobs.py` | job matrix and cost model → `jobs/*.txt` |
| `run_queue.py` | runs a job file on several GPUs, skips `DONE`, resumes, `--shard i/n` for several instances |
| `setup_vast.sh` | instance setup (keeps the image torch, checks sm_120 support for RTX 5090, runs the tests) |
| `test_study.py` | CPU tests: geometry, rotations, stress statistics and correlation length, stress term = −∂E/∂m (autograd), drive |

### Outputs of one run (`runs/<name>/`)

- `config.json`: all parameters, φ_vox, l_ex, git commit, GPU.
- `summary.json`: per stage and cycle: B_peak, H_peak, M_peak, W_loop, W_dis,
  closure, dW_rel, steps, wall time, mean dt.
- `samples_<k>_<stage>.csv`: 257 samples per cycle: t, H, M_par, Mx, My, Mz,
  B, p_dis, M per particle.
- `checkpoint.pt` (m after the last complete stage), `DONE`, `stdout.log`.

Definitions: B = μ0 (H_ext + ⟨M⟩_box) (core average, voids included).
W_loop = ∮ H dB = μ0 ∮ H dM. W_dis = ∫⟨p_dis⟩ dt with
p_dis = α γ μ0 Ms/(1+α²) |m × H_eff|². For a closed cycle W_loop = W_dis
exactly (all other energies are state functions). This is the energy
balance check.

## 6. Procedure

Do the steps in this sequence. Do not start step 4 before the gates G1 to G3
are passed or Chris accepts a deviation.

1. **Setup** (each instance): `bash setup_vast.sh`. It must end with
   "Setup complete". On an RTX 5090, the torch build must support sm_120
   (CUDA ≥ 12.8).
2. **Benchmark** (≈ 0.2 GPU-h): `python run_queue.py jobs/bench.txt --gpus 0`.
   Read `mean_dt` and `wall_s/steps` from `runs/bench_*/summary.json`.
   Put the measured values into `jobs.py` (`--s_per_step_Mcell`, `--dt_ps`,
   `--speedup`) and run `python jobs.py` again. Report the new total cost to
   Chris. **If the cost is more than 2× the estimate below, stop and ask.**
3. **Validation**: `jobs/validation.txt` on the fp32 GPUs and
   `jobs/validation_fp64.txt` on a V100/A100. Evaluate the gates (section 7).
4. **Production**: `python run_queue.py jobs/production.txt --gpus 0,...,7`.
   With several instances use `--shard i/n`. Copy `runs/` back (rsync) before
   you destroy an instance.
5. **Analysis**: `python analyze.py runs/P_* runs/V_*`, and
   `python analyze.py --pair runs/<30 MHz run> runs/V_f15_...` for the
   frequency separation. Write a short report (tables and plots) for Chris.

### Cost estimate (before the benchmark)

`python jobs.py` (V100 fp32 = 1): benchmark 0.2 h, validation fp32 ≈ 20 h,
validation fp64 ≈ 2.3 h, production ≈ 81 h, total ≈ 104 V100-hours.
Assumptions: 0.011 s per RKF45 step per 10⁶ cells (V100 fp32, from the fp32
validation run), mean dt = 2.0 ps at 10.26 nm (CPU smoke test: 2.4 ps to
3.0 ps), 3.3 cycles per amplitude.
RTX 5090: expected speed-up ≈ 2 to 2.5 → ≈ 45 GPU-hours → 8 GPUs ≈ 6 h wall
time, ≈ 15 $ to 30 $ on vast.ai (marketplace prices change).

## 7. Validation gates

| Gate | Test | Pass criterion | If it fails |
|---|---|---|---|
| G0 | `pytest test_study.py` | all pass | fix the code |
| G1 | benchmark | cost ≤ 2× estimate | ask Chris (options: fewer amplitudes, smaller d, fewer seeds) |
| G2 | energy balance, all steady cycles | \|W_dis/W_loop − 1\| < 2 % | increase `--samples`; check closure |
| G3 | grid pinning: `V_dx7_*` vs `P_d050_N072_*` | W(B) agrees within 10 % to 20 % (take the cycle scatter into account); the σ = 0 run gives the lower limit of numerical hysteresis | if not: tell Chris; options: 7 nm for all runs (≈ 6× cost) or report the grid limit |
| G4 | fp32 vs fp64 (`V_fp64_*` vs `P_d050_N072_*`) | W at the 2 smallest amplitudes within 5 % (or within the cycle scatter) | repeat the smallest amplitudes in fp64 from the fp32 state (see note N6) |
| G5 | quasi-static: 15 MHz vs 30 MHz (`--pair`) | dynamic share at 30 MHz is known; if > 30 % at small B, use W_h from the separation for β | add a 15 MHz twin for the affected runs (ask Chris, cost) |
| G6 | protocol: ascending vs descending | W(B) agrees within the cycle and seed scatter | report; use the descending results |
| G7 | steadiness | most measured stages have `steady = true` | increase `--max_cycles_per_amp` |

## 8. Results of the CPU smoke tests (for orientation)

Tiny box (N = 20, dx = 10 nm, d = 130 nm, φ = 0.54, σ_rms = 150 MPa, 300 MHz).
Not physically relevant, only a code test.

- The code runs in fp32 and fp64; resume from a checkpoint reproduces the
  cycle exactly.
- Energy balance of a closed cycle: 64 samples → 4 % error, 256 samples →
  0.2 %, 1024 samples → 0.1 %. Keep 256.
- Mean RKF45 step: 2.4 ps to 3.0 ps at dx = 10 nm (A = 15 pJ/m, Js = 1.8 T).
- fp32 and fp64 agree to 4 digits until the first switching event. After it,
  the paths separate (chaotic), but W stays in the same range (3 %).
- **Cycle-to-cycle scatter is large** for 4 small particles: W changes by
  10 % to 18 % from cycle to cycle at medium amplitude, also when the
  closure is small (the loop takes different metastable paths). Ascending and
  descending protocols gave 10 % to 15 % different W. For this reason
  `analyze.py` averages all cycles except the first and gives a standard
  error, and the matrix has 3 seeds for the key points. Expect less scatter
  with larger particles, but do not assume it.

## 9. Numerical notes (gotchas)

- N1. `set_precision()` must come before the first `Mesh`. `run_loops.py` does it.
- N2. **Time in fp64.** `state.t` is a float64 tensor on purpose. In fp32,
  t + dt with t ≈ 100 ns and dt ≈ 1 ps has a rounding error of ≈ 1 % per
  step. The drive reads `state.t`; the logged H uses the same `state.t`.
  `state.t` is reset to 0 at each stage start.
- N3. Void cells: Ms must be > 0 (field terms divide by Ms), A = 0 stops the
  exchange coupling (harmonic mean of A). The contact check in
  `geometry.py` makes sure that no two particles touch through a face.
- N4. `DemagFieldPBC` removes the k = 0 mode: the mean demag field is 0, so
  H_ext is the macroscopic H. Do not use `DemagField` with pseudo-PBC images
  here (it gives a macroscopic shape demag field and has the known
  `zero_origin` far-field defect).
- N5. The first call of each term triggers `torch.compile` (some minutes on
  a new instance). This is in the wall time of the first cycle.
- N6. A run resumes automatically from `checkpoint.pt` (m and the next stage
  index). To repeat single amplitudes from a saved state (e.g. in fp64 for
  G4), use a NEW output folder:
  `run_loops.py <same args> --init_from runs/X/checkpoint.pt --n_sat 0
  --H_list <amplitudes> --precision double --out runs/X_fp64tail`.
  Use the checkpoint of the stage BEFORE the amplitude: runs with
  `--keep_stage_ckpt` (the N = 72 production runs) keep `ckpt_<k>_<stage>.pt`
  after every stage.
- N7. Material parameters use copy-on-write (fork feature). Do not write
  in-place to a homogeneous parameter without `Material.materialize()`.
- N8. The stress field is set to 0 in the void and normalised on the
  magnetic cells only. It does not satisfy div σ = 0.
- N9. Orientation statistics of 4 crystals (Monte Carlo): E_a/K1 along the
  field = 0.20 ± 22 % (1 crystal: ± 44 %), nearest-easy-axis cos = 0.83 ± 6 %.
  `config.json` → `axes_stats` records the values of each run. The
  orientation is important mainly at large B (Q = K1/K_d ≈ 0.02). A factor
  that does not depend on B cancels in β = d ln W/d ln B.
- N10. Statistical error of β: with a relative scatter δ of W per amplitude
  and an amplitude ratio r, δβ ≈ √2 δ / ln r. δ = 10 %, r = 2 → δβ ≈ 0.2.
  This is large compared to the β range 2 … 3. Watch `W_err` in the
  validation runs; if δβ > 0.1, add cycles at the important amplitudes (not
  a supercell, see section 3).

## 10. Known limitations

- 4 particles per box: small ensemble, large scatter (see 8).
- d ≤ 1 µm instead of 1 µm to 5 µm: the size effect is an extrapolation of the d-scaling.
- T = 0: no thermal activation over barriers (W_h at small B slightly too high).
- No eddy currents (small at 1 µm; add analytically: P = π² σ_el d² f² B²/20).
- Isotropic λs, Gaussian stress model, identical spheres, no particle contacts.
- 10 nm cells: vortex cores (≈ 2 to 3 l_ex = 7 nm to 10 nm) are under-resolved →
  nucleation fields at large B are only qualitative.

## 11. Open points for Chris (do not decide alone)

1. Actual FeSiCr composition and VSM value of Ms (replace 1.80 T).
2. XRD micro-strain before/after annealing → σ_rms values.
3. λs of the alloy (or λ100, λ111).
4. After the benchmark: accept the budget, or reduce (fewer amplitudes,
   factor 0.7 only in the interesting B range, fewer seeds).
5. Optional: Py reference (λs ≈ 0, K1 = 0) with d up to 2 µm at 10 nm cells, to
   see β from magnetostatics and nucleation only.
6. **Positive control (proposed, not decided):** with λs = 5·10⁻⁶ the stress
   anisotropy is only K_σ,rms = 1.1 kJ/m³ (150 MPa) … 3 kJ/m³ (400 MPa) ≪ K1.
   Estimated pinning curvature Δγ_w/ξ² ≈ 2 … 5·10¹⁰ N/m³ < k_ms ≳ 10¹¹ N/m³.
   The model thus predicts only a small σ effect for FeSiCr. Proposal: 1 run
   with λs = 25·10⁻⁶ (FeNi50-like), σ_rms = 400 MPa, d = 1 µm (≈ 4 V100-h) to
   show that the mechanism works in the model at all.

## 12. Expected deliverables from Claude Code

- Benchmark numbers and the updated cost table (before production).
- Gate table G1 to G7 with the measured values.
- `analyze.py` outputs: `all_runs.csv`, `beta_vs_B.png`, frequency separation tables.
- A short report (tables + plots) that answers H1 to H4, with error bars.
- Commit code changes and small result files (csv, png, json) to this branch.
  Do not commit `samples_*.csv` or `checkpoint.pt` of all runs (large).
