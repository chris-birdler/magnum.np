# Handover: β(B_peak) study with magnum.np (fp32, reduced units)

Owner: Christoph Vogler (Chris). Branch: `study/beta-bpeak-pinning` (based on
`perf/memory-optimizations`, the fp32 fork). Folder: `studies/beta_bpeak/`.

Rules: documents in ASD-STE100 Simplified Technical English, SI units only.
**Chris makes the decisions.** Before you change the model, the parameter
matrix or the budget, give him the options with pros, cons, chances and risks,
and wait. Engineering fixes (bugs, logging, speed) need no decision, but tell him.

## 1. Goal

Find which parameters change the Steinmetz exponent β(B_peak) most, in a
**fictitious** soft magnetic powder material. No specific alloy. Practical aim:
understand how to keep β small up to high B_peak.

## 2. Model (reduced units)

Periodic FCC unit cell, 4 spheres, `DemagFieldPBC` (mean demag field = 0,
like a toroid core: H_ext = H_macro), exchange, anisotropy, random stress
anisotropy, homogeneous sinusoidal drive along x. T = 0, LLG, RKF45, fp32.

The physics depends only on dimensionless numbers. Js and A only select the SI
unit system (`units.py`):

| Unit | Formula | Value (Js = 1.5 T, A = 10 pJ/m) |
|---|---|---|
| length l_ex | √(2A/(μ0 Ms²)) | 3.34 nm |
| energy density K_d | μ0 Ms²/2 | 895 kJ/m³ |
| field | Ms | 1.19 MA/m |
| frequency f_M | γ Ms/2π (28 GHz/T · Js) | 42.0 GHz |

| Reduced parameter | Meaning | Base value |
|---|---|---|
| d/l_ex | particle diameter | 300 (≈ 1 µm) |
| φ | packing fraction | 0.65 |
| Q = K/K_d | anisotropy, uniaxial | 0.02 (δ0 = √(A/K) = 7.1 l_ex) |
| κ = K_σ,rms/K_d | random stress anisotropy (Gaussian tensor field, energy −K_σ Σ s_ij m_i m_j, s_ij unit rms) | 0.008 |
| ξ/l_ex | correlation length of the stress field | 15 |
| dx/l_ex | cell size (numerical) | 3 |
| f/f_M | drive frequency | 7.14·10⁻⁴ (= 30 MHz) |
| α | damping (numerical choice) | 0.1 |

Output: w = W/K_d per cycle versus b = B_peak/Js, β_eff = d ln w/d ln b.

Physics background (short): at small b all losses are linear → β ≈ 2; Rayleigh
hysteresis (∝ b³) raises β towards 3; at large b the loop saturates → β falls.
In 1 µm particles the wall stiffness is magnetostatic, k_ms ≈ 4 μ0 N Ms²/d.
Only short-range fluctuations of the wall energy (scale ξ) can pin against
k_ms (Barkhausen criterion). The κ scan must cross this threshold; the old
FeSiCr values (κ ≈ 0.001) were below it.

## 3. Decisions (with reasons)

| Item | Decision | Reason |
|---|---|---|
| Material | fictitious, reduced units, SI unit system Js = 1.5 T, 30 MHz | general answer, not one alloy |
| Size | real regime: d/l_ex = 150 … 300 (0.5 … 1 µm) | small particles (d/l_ex ≤ 100) have other physics (vortex states, high surface fraction) |
| Anisotropy | Q = 0.02, **uniaxial**, deterministic easy axes cos θ = 1/8, 3/8, 5/8, 7/8 to the drive | K > 0 gives real walls with width δ0 ≫ l_ex; deterministic axes remove the orientation statistics problem (4 random axes: ±22 % scatter). Cubic is possible with `--aniso cubic` |
| Mesh | dx = 3 l_ex | walls resolved: δ0/dx = 2.4, grid (Peierls) pinning ∝ exp(−π² δ0/dx) ≈ 10⁻¹⁰. l_ex-scale structures (vortex cores, Bloch points) are NOT resolved; accepted (Bloch points are singular in micromagnetics anyway) |
| No supercells, no K = 0, no reduced size | rejected | same statistics per cost as more cycles; K = 0 or small d change the regime |
| Protocol | 1 saturating cycle, 9 descending amplitudes (factor 0.5), 3 … 4 cycles each, first cycle discarded | AC demagnetization; cycle-to-cycle scatter is large (10 … 18 % in a CPU test) |
| GPU | any fp32 GPU; choose RTX 3090 or RTX 5090 after the benchmark (cost per cycle) | "V100-h" is only a cost unit |

## 4. Mesh validity (what we accept and how we see it)

- Walls: resolved (see above).
- Unresolved structures exist anyway. CPU test at dx = 3 l_ex: largest
  neighbour angle ≈ 180° (a vortex core is narrower than one cell).
- Diagnostics in every run: `frac_pairs_gt60` per cycle and `n_pairs_gt60` per
  sample (csv). Jumps show nucleation/annihilation of cores or Bloch points.
  `analyze.py` flags such amplitudes with `#`. Interpret these points with care.
- The κ = 0 run gives the floor (numerical + magnetostatic hysteresis without
  pinning). The dx = 2 l_ex controls show how much of the result depends on
  the mesh. These are the only mesh checks.

## 5. Code map

| File | Content |
|---|---|
| `units.py` | reduced ↔ SI units, default f/f_M |
| `run_loops.py` | one run (reduced inputs), protocol, diagnostics, checkpoints, resume |
| `geometry.py` | FCC box, voxelized spheres, contact check, random rotations |
| `stress.py` | Gaussian random tensor field with correlation length ξ |
| `field_terms_extra.py` | `StressAnisotropyField`, `SinusoidalDrive` |
| `jobs.py` | job matrix and cost model → `jobs/*.txt` |
| `run_queue.py` | runs a job file on several GPUs, skips `DONE`, resumes, `--shard i/n` |
| `analyze.py` | w(b), β_eff(b), flags, features for the ranking, `--pair` frequency separation |
| `setup_vast.sh` | instance setup and tests |
| `test_study.py` | CPU tests |

Outputs of one run (`runs/<name>/`): `config.json` (all parameters, units,
N, φ_vox), `summary.json` (per cycle: b_peak, w_loop, w_dis, closure, dW_rel,
frac_pairs_gt60, steps, wall time, mean dt, and SI values),
`samples_<k>_<stage>.csv` (257 samples per cycle, SI), `checkpoint.pt`, `DONE`.

Energy balance: for a closed cycle w_loop = w_dis exactly
(w_dis = ∫ α γ μ0 Ms/(1+α²) |m × H_eff|² dt / K_d). CPU test: 0.2 % with 256 samples.

## 6. Procedure

1. `bash setup_vast.sh` on each instance (RTX 5090 needs torch with CUDA ≥ 12.8).
2. Benchmark: run `jobs/bench.txt` (1 cycle of the base point, ≈ 0.15 V100-h)
   on one RTX 3090 and one RTX 5090. Compare `wall_s` per cycle × price.
   Update `jobs.py` (`--s_per_step_Mcell`, `--dt_tau`) with the measured
   values and run `python jobs.py --speedup <x>`. Report the cost to Chris.
   **Stop and ask if the cost is more than 2× the estimate below.**
3. Run `jobs/production.txt` and `jobs/control.txt` together:
   `python run_queue.py jobs/production.txt --gpus 0,1,2,3` (and the control
   file on other GPUs, or `--shard i/n` on several instances). Copy `runs/`
   back before you destroy an instance.
4. `python analyze.py runs/P_* runs/C_*` and
   `python analyze.py --pair runs/P_kappa00800 runs/C_fhalf`.
5. Report to Chris: ranking of the parameters (change of β@b and β_max per
   parameter step, with error bars), plots, checks (section 8).

## 7. Matrix and cost

Base point: d/l_ex = 300, φ = 0.65, Q = 0.02, κ = 0.008, ξ/l_ex = 15, dx/l_ex = 3.

| Block | Runs |
|---|---|
| κ = 0 / 0.0025 / 0.008 / 0.025 | 4 |
| d/l_ex = 150 / 225 | 2 |
| φ = 0.47 / 0.71 | 2 |
| Q = 0.01 | 1 |
| seeds 2, 3 (stress field, initial m) | 2 |
| control: κ = 0 at d/l_ex = 150 (dx 3 and 2), κ = 0.008 at d/l_ex = 150 with dx 2, f/f_M halved | 4 |

Estimate (`python jobs.py`, before the benchmark): production ≈ 44 V100-h,
control ≈ 18 V100-h, total ≈ 62 V100-h. RTX 3090 ≈ 1× V100, RTX 5090 ≈ 2 … 2.5×.
On vast.ai ≈ 7 $ to 18 $ in total.

## 8. Checks (no extra runs)

| Check | Pass criterion | If it fails |
|---|---|---|
| energy balance (`*` flag) | \|w_dis/w_loop − 1\| < 2 % | more `--samples` |
| steadiness (`!` flag) | most stages steady | more `--max_cycles_per_amp` |
| floor | w(κ = 0) ≪ w(κ > 0) at small and medium b | the pinning signal is below the numerical floor: report, ask Chris |
| mesh | C_d150_dx2 vs P_d150 and the two κ = 0 runs: w(b) within the cycle scatter | report the mesh dependence; ask Chris |
| quasi-static | `--pair`: dynamic share at 30 MHz; if > 30 % at small b, use β from w_h | ask Chris |

## 9. Known limitations

4 particles per box (large scatter); T = 0; no eddy currents (add
analytically); identical spheres, no contacts; Gaussian stress model (no
div σ = 0); isotropic magnetostriction model; l_ex-scale structures and
Bloch points not resolved (flagged); 5 µm particles not reachable (only the
scaling from d/l_ex = 150 … 300).

## 10. Numerical notes

- `set_precision()` before the first `Mesh` (done in `run_loops.py`).
- `state.t` is float64 on purpose (fp32 time gives ≈ 1 % error per step at
  t ≈ 100 ns); reset to 0 at each stage.
- Void: Ms_void = 10⁻⁸ Ms (> 0, field terms divide by Ms), A = 0 (no exchange
  across the gap); `geometry.py` rejects face contact between particles.
- First call of each term triggers `torch.compile` (minutes).
- Material parameters are copy-on-write (fork feature): use
  `Material.materialize()` before in-place writes to homogeneous parameters.

## 11. Open points for Chris

1. After the benchmark: GPU type and budget.
2. After the first results: refine the κ range or add ξ/l_ex values
   (3 … 30) if the ranking shows that pinning dominates.
3. Optional: cubic anisotropy comparison (`--aniso cubic`) at the base point.
