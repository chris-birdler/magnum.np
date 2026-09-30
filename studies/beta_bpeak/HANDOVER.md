# Handover: β(B_peak) study with magnum.np (fp32, reduced units)

Owner: Christoph Vogler (Chris). Branch: `study/beta-bpeak-pinning` (based on
`perf/memory-optimizations`, the fp32 fork). Folder: `studies/beta_bpeak/`.

Rules: documents in ASD-STE100 Simplified Technical English, SI units only.
**Chris makes the decisions.** Before you change the model, the parameter
matrix or the budget, give him the options with pros, cons, chances and risks,
and wait. Engineering fixes (bugs, logging, speed) need no decision, but tell him.

## 1. Goal

Find the parameters that control the Steinmetz exponent β(B_peak) of a soft
magnetic powder core. The material is a **fictitious nanocrystalline**
powder. It lets us separate the influences cleanly. Practical aim: understand
how to keep β small up to high B_peak.

Fixed by Chris:
- particle diameter ≈ 1 µm,
- Bloch points and vortex cores (l_ex scale) are not resolved (accepted),
  on the condition that the mesh test (section 7) passes,
- no induced anisotropy (a pressed powder core has no common texture).

## 2. Model (reduced units)

Periodic FCC unit cell, 4 spheres, `DemagFieldPBC` (mean demag field = 0,
like a toroid core: H_ext = H_macro), exchange, anisotropy (below),
homogeneous sinusoidal drive along x. T = 0, LLG, RKF45, fp32.

Js and A only select the SI unit system (`units.py`):

| Unit | Formula | Value (Js = 1.5 T, A = 10 pJ/m) |
|---|---|---|
| length l_ex | √(2A/(μ0 Ms²)) | 3.34 nm |
| energy density K_d | μ0 Ms²/2 | 895 kJ/m³ |
| field | Ms | 1.19 MA/m |
| frequency f_M | γ Ms/2π (28 GHz/T · Js) | 42.0 GHz |

### 2.1 Anisotropy (`anisotropy_model.py`)

There are two sources. Both are uniaxial. There is no other anisotropy.

1. **Herzer residual anisotropy.** The grains (nm scale) are not resolved.
   Exchange averages their anisotropy (Herzer): K_eff = K1⁴D⁶/A³. The
   averaged anisotropy has a random direction on the scale of the exchange
   length L_eff = √(A/K_eff). The model uses cubes with edge L_eff. Each cube
   has one random isotropic easy axis and the amplitude K_eff. In reduced
   units, Q_eff = K_eff/K_d = (l_ex/L_eff)². The model has one free
   parameter here: L_eff (or Q_eff).
   - Translation to a material: Q_eff = Q1 (D/δ1)⁶ with Q1 = K1/K_d and
     δ1 = l_ex/√Q1. Example: K1 = 48 kJ/m³, D = 10 nm → Q_eff = 6·10⁻³.
   - Cube edge = L_eff is a convention. Herzer fixes L_eff only up to a
     factor of order 1. This factor moves all points in the same way.
   - Anisotropy fluctuations on scales below L_eff (grains, dislocation
     stresses) are averaged into Q_eff. They are not a separate factor.
2. **Residual stress on the particle scale.** One uniaxial anisotropy
   K_p = r_p K_eff per particle. The easy axes are deterministic:
   cos θ = 1/8, 3/8, 5/8, 7/8 to the drive (midpoint rule for an isotropic
   powder). There is no orientation statistics.

### 2.2 Parameters

| Reduced parameter | Meaning | Values |
|---|---|---|
| L_eff/l_ex | Herzer length = cube edge; Q_eff = (l_ex/L_eff)² | 12 / 18 / 24 / 30 → Q_eff = 6.9 / 3.1 / 1.7 / 1.1 ·10⁻³ |
| r_p = K_p/K_eff | particle-scale stress anisotropy | 0 / level from the pilot |
| φ | packing fraction | 0.55 … 0.69 |
| d/l_ex | particle diameter | 300 (≈ 1 µm), fixed |
| dx/l_ex | cell size (numerical) | 3 (mesh test: 2 and 1.5) |
| f/f_M | drive frequency (numerical, quasi-static aim) | 7.14·10⁻⁴ (= 30 MHz) |
| α | damping (numerical) | 0.02 (pilot compares 0.1) |

Output: w = W/K_d per cycle versus b = B_peak/Js, β_eff = d ln w/d ln b.

Physics background (short):
- At small b, reversible rotation and damping give w ∝ b², thus β ≈ 2.
  The damping part is proportional to α f. It is an artifact of the
  numerical α and the frequency pair controls it.
- Walls pin at the anisotropy fluctuations on the scale L_eff. This gives
  Rayleigh hysteresis (w ∝ b³), thus β rises towards 3.
- At large b the loop saturates, thus β falls.
- The Herzer parameter x = D/δ1 of the real material sets Q_eff = Q1 x⁶.
  For x ≪ 1 the anisotropy disappears. For d/L_eff < 5 a 1 µm particle has
  no walls, and magnetostatics (vortex states) controls the loss.

## 3. Decisions (with reasons)

| Item | Decision | Reason |
|---|---|---|
| Material | fictitious nanocrystalline powder, reduced units, SI unit system Js = 1.5 T, 30 MHz | general answer, clean separation |
| Anisotropy | only Herzer cubes (Q_eff) + particle-scale stress (r_p) | a global K per particle mixes with the pinning source; short-range stress merges into Q_eff |
| Cubes | equal cubes of edge L_eff, not Voronoi | Herzer assumes one grain size; no size distribution as a hidden parameter; cube faces lie on cell faces (no staircase pinning) |
| Q_eff window | 1.1·10⁻³ … 6.9·10⁻³ | L_eff ≥ 4 cells (resolved); d/L_eff ≥ 10 (walls in the particle, ≥ 500 cubes per particle, small realisation scatter) |
| Mesh-identical cubes | box edge a = multiple of 12 l_ex, L_eff = multiple of 6 l_ex | dx = 3, 2 and 1.5 see exactly the same cubes (tested on a common fine grid, with a negative control); d is adjusted by ≤ 2.5 % to keep φ exact |
| Contacts | at least one void cell between particles (no face, edge or corner contact) | no exchange coupling, no extreme stray fields at contacts; at d/l_ex = 150 and dx = 3, φ = 0.69 and 0.71 are not possible |
| Size | d/l_ex = 300 (1 µm), fixed | Chris |
| Scatter target | 10 % scatter of w from cycle to cycle (`--dW_tol 0.10`) | gives σ(β) ≈ 0.07 per point (section 8) |
| Protocol | 1 saturating cycle, descending amplitudes, 3 … 4 cycles each, first cycle discarded | AC demagnetization |
| Damping | α = 0.02 | 5× less artificial damping loss than α = 0.1 at the same cost per step; the pilot checks it |
| GPU | any fp32 GPU; choose RTX 3090 or RTX 5090 after the benchmark (cost per cycle) | "V100-h" is only a cost unit |

## 4. Mesh validity

| Structure | Status |
|---|---|
| walls (width ≥ L_eff ≥ 12 l_ex = 4 cells) | resolved; grid pinning ∝ exp(−π² δ/dx) ≈ 10⁻¹⁷ |
| anisotropy cubes (≥ 4 cells) | resolved; faces on cell faces |
| particle-scale stress | uniform per particle; local wall width 1/√(Q_eff (1 + r_p)) must stay ≥ 2 cells (`run_loops.py` warns) |
| vortex cores, Bloch points (l_ex scale) | NOT resolved. An audit test (2D disk with a vortex, Q = 0, quasi-static) found grid pinning of the core at dx = 3 l_ex: branch gap up to 0.09 in m, h_c,eff = 1.4·10⁻³. At dx = 1.5 l_ex the loop has no hysteresis. This can be larger than the Herzer pinning (estimate ≈ 5·10⁻⁵). The mesh test (section 7) decides |

Diagnostics in every run: `frac_pairs_gt60` per cycle and `n_pairs_gt60` per
sample (csv). `analyze.py` flags amplitudes with jumps (`#`). The control
T_L30_dx2 (same cubes as T_L30) shows how much the result depends on the mesh.

## 5. Code map

| File | Content |
|---|---|
| `units.py` | reduced ↔ SI units, default f/f_M |
| `anisotropy_model.py` | box quantization, Herzer cubes (mesh independent), particle axes |
| `geometry.py` | FCC box, voxelized spheres, contact check (26-neighbourhood) |
| `field_terms_extra.py` | `SinusoidalDrive` |
| `run_loops.py` | one run (reduced inputs), protocol, diagnostics, atomic checkpoints, resume with argument check |
| `jobs.py` | job matrix (bench, meshtest, pilot) and cost model → `jobs/*.txt` |
| `run_queue.py` | runs a job file on several GPUs, skips `DONE`, resumes, `--shard i/n` |
| `analyze.py` | w(b), β_eff(b) with errors, flags, features, ranking against the base point, `--compare` (two runs at the same b), `--pair` frequency separation |
| `setup_vast.sh` | instance setup (clone or update to the pushed branch, print the commit) and tests |
| `test_study.py` | CPU tests |

Outputs of one run (`runs/<name>/`): `config.json` (all parameters, units,
N, φ_vox, anisotropy data), `summary.json` (per cycle: b_peak, w_loop, w_dis,
closure, dW_rel, frac_pairs_gt60, steps, wall time, mean dt, SI values; per
stage: steady, closed, git commit), `samples_<k>_<stage>.csv` (257 samples per
cycle, SI), `checkpoint.pt`, `DONE`.

Resume: `run_loops.py` stops with an error if the physics or protocol
arguments differ from the run that wrote the checkpoint. Checkpoint and
summary are written atomically.

Energy balance: for a closed cycle w_loop = w_dis
(w_dis = ∫ α γ μ0 Ms/(1+α²) |m × H_eff|² dt / K_d). CPU smoke test: 1 %
with 32 samples. w_dis comes from point samples and can miss short
Barkhausen jumps. Thus the `*` flag is information only.

## 6. Procedure

1. Commit and push the branch. `setup_vast.sh` takes the pushed state.
2. `bash setup_vast.sh` on each instance (RTX 5090 needs torch with CUDA ≥ 12.8).
3. Benchmark: run `jobs/bench.txt` (1 cycle at d/l_ex = 300, ≈ 0.15 V100-h)
   on one RTX 3090 and one RTX 5090. Compare `wall_s` per cycle × price.
   Update `jobs.py` (`--s_per_step_Mcell`, `--dt_tau`) with the measured
   values. Report the cost to Chris.
4. Mesh test FIRST (done, section 7.1), then the steady test:
   `python run_queue.py jobs/steady.txt --gpus 0,1`. Evaluate it (section 7.2)
   and report to Chris. Stop if the mesh check fails.
5. Pilot: `python run_queue.py jobs/pilot.txt --gpus 0,1,2,3`. Copy `runs/`
   back before you destroy an instance.
6. `python analyze.py runs/T_*` and `python analyze.py --pair runs/T_L18 runs/T_L18_f2`.
   Report the pilot checks (section 7.3) to Chris.
7. Chris decides the production matrix (section 9). Then add it to `jobs.py`.

Never edit tracked files on an instance: `setup_vast.sh` discards such edits,
and a resume stops if the code commit changed (`--allow_new_commit` overrides).
Commit and push instead. The commit of each run is in `config.json` and in
every stage of `summary.json` ("-dirty" = uncommitted changes).

## 7. Mesh test and pilot

### 7.1 Mesh test (d/l_ex = 96, ≈ 6 V100-h, before the pilot)

Question: does the grid pin the unresolved vortex cores in the powder model
(spheres, periodic demag, Herzer cubes)? A small particle is a conservative
test: core pinning grows as the particle becomes smaller.

| Runs | dx/l_ex | Content |
|---|---|---|
| M_floor_dx3 / dx2 / dx15 | 3 / 2 / 1.5 | Q_eff = 0: only magnetostatics and exchange |
| M_L12_dx3 / dx2 / dx15 | 3 / 2 / 1.5 | L_eff = 12 l_ex (same cubes on all meshes) |

a = 144 l_ex, d = 97.5 l_ex, φ = 0.65, 4 amplitudes, α = 0.02.

Evaluation (tolerances fixed before the runs):

    python analyze.py --compare runs/M_floor_dx3 runs/M_floor_dx15
    python analyze.py --compare runs/M_L12_dx3  runs/M_L12_dx15
    python analyze.py --compare runs/M_L12_dx2  runs/M_L12_dx15

| Result | Meaning | Next step |
|---|---|---|
| PASS for dx = 3 vs 1.5 (level: 90 % interval of ln w_3/w_1.5 inside ±0.15; slope: 90 % interval of Δβ inside ±0.2) for both floor and L12 | the grid does not change w(b) and β | pilot with dx = 3 |
| FAIL for dx = 3, PASS for dx = 2 | dx = 3 is too coarse, dx = 2 is sufficient | report the cost of dx = 2 (≈ 7.6× per run: 3.4× cells, 2.25× steps) to Chris |
| FAIL for both | the cores must be resolved for this observable. Chris's condition "cores not resolved" does not hold | stop, report options to Chris (finer mesh at smaller d/l_ex, or a fictitious material with a larger l_ex) |

Also compare `frac_pairs_gt60` and the per-sample `n_pairs_gt60` jumps
between the meshes. Steps in m(t) at dx = 3 that disappear at dx = 1.5 are
grid pinning of cores.

**Result (2026-09-30, 2× V100, ≈ 0.40 $; raw data local with MANIFEST.md5):**
- b ≈ 0.74: steady and clean. L12: dx = 3 vs 1.5 differs by +3 %. Floor: +24 %.
- b ≤ 0.25: NOT decided. At α = 0.02 the cycles are not steady within 3 … 4
  cycles (closure up to 22 %, negative loop areas, w_dis/w_loop = 1.3 … 4.5,
  scatter 20 … 80 %). Thus `--compare` fails, but this is no proof of a mesh
  effect.
- dx = 3 keeps ≈ 24 neighbour pairs above 60° (unresolved cores); dx = 1.5
  has 0 … 3. At the smallest amplitude dx = 1.5 dissipates 5 … 10× more
  (w_dis). This points to cores that the coarse grid holds in place.
- Consequence: 3 … 4 cycles per amplitude are not sufficient at α = 0.02.

### 7.2 Steady test (d/l_ex = 96, ≈ 5.5 V100-h, before the pilot)

Same geometry as 7.1. Two amplitudes (h = 0.0292 and 0.0073, b ≈ 0.25 and
0.06), 12 cycles each, no early stop.

| Runs | dx/l_ex | Content |
|---|---|---|
| S_floor_dx3 / dx15 | 3 / 1.5 | Q_eff = 0, α = 0.02 |
| S_L12_dx3 / dx15 | 3 / 1.5 | L_eff = 12, α = 0.02 |
| S_L12_a010_dx3 / dx15 | 3 / 1.5 | L_eff = 12, α = 0.1 (faster relaxation, more damping loss) |

Evaluation (the last 6 of 12 cycles):

    python analyze.py --skip 6 --cycles runs/S_*
    python analyze.py --skip 6 --compare runs/S_L12_dx3 runs/S_L12_dx15
    python analyze.py --skip 6 --compare runs/S_floor_dx3 runs/S_floor_dx15
    python analyze.py --skip 6 --compare runs/S_L12_a010_dx3 runs/S_L12_a010_dx15

1. Steady? For the kept cycles: |drift of w| < 1 % per cycle (within 2
   errors), closure < 2 %, w_dis/w_loop within 0.9 … 1.1. If not: report the
   number of cycles that is necessary (cost) to Chris.
2. Mesh: the same tolerances as in 7.1 (level ±0.15, slope ±0.2, 90 %).
3. Damping: compare α = 0.02 and 0.1 (steady time and loss level).

### 7.3 Pilot (d/l_ex = 150, 4 amplitudes, ≈ 4.8 V100-h)

The pilot uses d/l_ex = 150 (8× cheaper). There are 60 … 940 cubes per
particle. Thus the realisation scatter is larger than in production.

| Run | Question | Pass criterion | If it fails |
|---|---|---|---|
| T_L12 … T_L30 vs T_floor (Q_eff = 0) | is there a pinning signal above the floor? | `--compare`: level of T_L above T_floor by more than 3 σ at small and medium b | the signal is below the floor: report, ask Chris |
| T_L30_dx2 vs T_L30 | mesh dependence (same cubes) | difference within the realisation scatter (T_L30 vs T_L30_seed2) | report the mesh dependence, ask Chris |
| T_L30_seed2 vs T_L30 | realisation scatter at the weakest point | — (input for the error budget) | — |
| T_L18_alpha010 vs T_L18 | size of the damping loss | — (shows the gain of α = 0.02) | — |
| T_L18_f2 vs T_L18 (`--pair`) | dynamic share at α = 0.02 | dynamic share < 30 % or less than 2 errors at small b | ask Chris (lower f or use β from w_h) |
| T_L18_d111 vs T_L18 | do the cube faces cause a preferred direction? | difference within the scatter | ask Chris (option: equal truncated octahedra) |
| T_L18_rp1, T_L18_rp3 | does the particle-scale stress change β? | — (Chris selects the r_p level) | — |
| all | scatter | s_cyc ≤ 0.10 | more cycles cost budget: ask Chris |

## 8. Error budget

| Quantity | Error (10 % scatter) | How `analyze.py` gets it |
|---|---|---|
| ln w at one amplitude | max(s_cyc, 0.10)/√n ≈ 0.07 (n = 2 … 3 measured cycles) | s_cyc pooled over all amplitudes of one run |
| β_eff at one amplitude | ≈ 0.07 (d/l_ex = 300, 9 amplitudes) | central difference over the two neighbour amplitudes |
| realisation | from the seeds of the base point | variance over the seeds minus the statistical part |
| change of a feature | √(stat.² + realisation² + base²) | `ranking()`; significance with Student-t of the seed scatter |
| dynamic share (`--pair`, 2f) | ≈ ±0.10 | two runs with 10 % scatter each |

Consequences:
- A change of β smaller than ≈ 0.2 is not a result.
- β_max from the central difference is smoothed (up to −0.3 in synthetic
  tests). Use β at fixed b for the ranking.
- With 5 features and many runs, some |d|/σ above the threshold are false
  hits. Report them as candidates.

## 9. Production matrix (open: Chris decides after the pilot)

Base: d/l_ex = 300, φ = 0.65, L_eff = 18 l_ex, r_p = 0. One d/l_ex = 300 run
costs ≈ 5 V100-h.

| Option | Runs | Cost | Pros | Cons |
|---|---|---|---|---|
| A: full factorial Q_eff (4) × r_p (2) × φ (2) + 2 seeds | 18 | ≈ 90 V100-h | all main effects and interactions | budget |
| B: Q_eff scan (4) at r_p = 0, φ = 0.65 + 2×2 (r_p, φ = 0.55 / 0.69) at L_eff = 18 + 2 seeds | 10 | ≈ 50 V100-h | Q_eff curve and both other factors with their interaction | interactions with Q_eff not visible |
| C: B + floor (Q_eff = 0) at d/l_ex = 300 | 11 | ≈ 55 V100-h | B + floor at production size | cost |

A dx = 2 control at d/l_ex = 300 costs ≈ 38 V100-h. Thus the mesh check stays
at d/l_ex = 150 (pilot).

## 10. Known limitations

- 4 particles per box, identical spheres, FCC order, no contacts.
- T = 0 (no thermal activation). No eddy currents (add them analytically).
- Herzer length only up to a factor of order 1.
- No anisotropy on intermediate scales (≈ 300 nm), by choice.
- Vortex cores and Bloch points not resolved (flagged).

## 11. Numerical notes

- `set_precision()` before the first `Mesh` (done in `run_loops.py`).
- `state.t` is float64 on purpose (fp32 time gives ≈ 1 % error per step at
  t ≈ 100 ns); reset to 0 at each stage.
- Void: Ms_void = 10⁻⁸ Ms (> 0, field terms divide by Ms), A = 0 (no exchange
  across the gap).
- Random streams: initial m uses [seed, 0], cube axes use [seed, 1000 + p].
- First call of each term triggers `torch.compile` (minutes).
- Material parameters are copy-on-write (fork feature): use
  `Material.materialize()` before in-place writes to homogeneous parameters.
- The RKF45 solver discards one trial step at each sample boundary (≈ 2 %
  cost). This is in the library, not changed.

## 12. Open points for Chris

1. After the benchmark: GPU type and budget.
2. After the pilot: the r_p level and the production matrix (section 9).
3. After the pilot: α = 0.02 confirmed or not.
4. Audit points that wait for the mesh test (then Chris decides):
   - separate the hysteresis loss w_h from the damping loss for every point
     (w_h = 2 w(f) − w(2f), +50 % cost), because the damping loss can move β by 0.2 … 0.5;
   - pilot controls at L18 instead of L30 (d/L_eff = 5 at L30 and d/l_ex = 150
     is the no-wall limit), and 2 more seeds at L18 (+0.75 V100-h);
   - pilot pass criteria as fixed tolerances with `--compare` (matched b);
   - ranking by a weighted regression over all runs (the present ranking
     compares only runs that differ in one parameter);
   - production option D: Q_eff (4) × r_p (2) at φ = 0.65 + φ = 0.55 / 0.69 at
     L18 for both r_p + 2 seeds = 14 runs, ≈ 70 V100-h (contains Q × r_p).
