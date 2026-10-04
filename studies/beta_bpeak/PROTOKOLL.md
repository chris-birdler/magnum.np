# Protokoll: β(B_peak) study with magnum.np (fp32, reduced units)

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
how to keep β small.

**Range of interest: B_peak = 10 … 150 mT** (b = 0.0067 … 0.1 with Js = 1.5 T).
The most frequent cases are 10, 50 and 100 mT. Higher B_peak is rare and not
part of the study. B_peak is the amplitude (half peak-to-peak), as in the
Steinmetz equation.

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
   - **K_eff is an input, not a result.** Each cube gets K_eff = Q_eff K_d
     (base point: L_eff = 12 l_ex = 40.1 nm, K_eff = 895 kJ/m³ / 144 =
     6.22 kJ/m³ = A/L_eff²). The simulation does not compute the averaging
     from the grains to K_eff; the Herzer formula does this before the run.
     The simulation computes the response of m to this cube landscape
     (walls, pinning, ripple, loss, β).
   - The cubes are Herzer's self-consistent fixed point: averaging cubes of
     edge L and K = A/L² over a larger edge L′ gives K′ = K (L/L′)^{3/2}, and
     L′ = √(A/K′) gives L′ = L. Thus a wall on the scale L_eff sees the full
     K_eff; there is no further averaging. The mean over a whole particle,
     K_eff/√N ≈ 64 J/m³ (N ≈ 9400 cubes), acts only on a uniform rotation of
     the particle and is small against K_d.
   - Herzer applied to the model itself (K1 = K_eff = 6.22 kJ/m³ per cube,
     D = cube edge = 40.1 nm, A = 10 pJ/m): δ1 = √(A/K1) = 40.1 nm, thus
     x = D/δ1 = 1.00 and K1⁴D⁶/A³ = 1.494·10¹⁵ · 4.164·10⁻⁴⁵ / 10⁻³³ =
     6.22 kJ/m³ = K1. There is no reduction: one cube per exchange volume.
     Herzer's averaging needs many grains per exchange volume (D ≪ L_ex).
     The model sits exactly at the edge of the Herzer regime (x = 1), at the
     transition to the classical grain regime; each cube acts with its full
     anisotropy. Two readings give the same simulation: (i) grains of 40 nm
     with K1 = 6.2 kJ/m³ and no averaging, (ii) grains of 10 nm with
     K1 ≈ 50 kJ/m³ after the averaging. For comparison, 10 nm grains with
     K1 = 6.2 kJ/m³ give x = 0.25, K_eff = 1.5 J/m³ and L_ex ≈ 2.6 µm > d.
   - **Assumption, not tested in the simulation:** the translation of Q_eff
     to a material (K1, D) holds only if Herzer's averaging holds. A check
     with explicit grains (one cell per 10 nm grain, K1 = 48 kJ/m³; Herzer
     predicts L_eff ≈ 43 nm ≈ base point) was offered (≈ 10 V100-h) and not
     run (decision Chris, 2026-10-03).
   - Range of the cubes in material terms: for K1 ≈ 20 … 50 kJ/m³ the window
     L_eff = 12 … 30 l_ex covers grains of D ≈ 10 nm; for a FINEMET-like
     K1 = 8 kJ/m³ and D = 10 nm, K_eff ≈ 4 J/m³ and L_eff ≈ 1.6 µm > d: the
     particle has almost no anisotropy, and only the particle-scale
     anisotropy (item 2) is left. This end is not covered by the cubes.
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
| α | damping (numerical) | 0.1 (steady cycles, section 7.2) |

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
| Amplitudes | B_peak = 150 / 100 / 70 / 50 / 35 / 20 / 9 mT (`--b_list`; 9 mT keeps 10 mT inside the range). The drive of each stage comes from b/h of the stage before and is corrected once after cycle 0 with b/h of cycle 0 (factor limited to 0.67 … 1.5) | study range 10 … 150 mT; no time on saturation. b of cycle 0 is within ±4 % of the later cycles (V100 data), thus each stage hits its target to ≈ 4 % |
| Drive | sinusoidal H (not controlled B) | below 150 mT the core is almost linear: B(t) has 0.2 … 3 % harmonics, b varies by 0.1 … 3 % from cycle to cycle (mesh test data). Thus sinusoidal H ≈ sinusoidal B (Steinmetz condition) |
| Initial state | open: decided by the initproto test (section 7.4). Candidates: A = virgin state (random m, full relaxation at H = 0 and α = 1: the T = 0 analogue of the anneal and the cooling without field; a real powder is never magnetized before use), B = AC demagnetization (decaying saturating cycles, as in IEC 60404-6), C = 1 saturating cycle (present) | C gave minor loops around a remanent state (offset up to 0.55 at 23 mT) |
| Protocol | A and B: ascending amplitudes (the classic Rayleigh procedure: every larger loop erases the smaller ones, the loops stay centred). C: descending. 7 cycles each (fixed), cycles 0 and 1 discarded (drive correction) | AC demagnetization. At T = 0 the soft powder does not lock into a periodic cycle (steady test): the mean over cycles is the measurement, not a single steady cycle |
| Loss per cycle | w = mean of w_loop (∮H dB) over the kept cycles; w_dis (LLG dissipation) as a check (`--w dis`) | w_loop of one cycle contains the change of the stored energy when the cycle is not closed; this part averages out over cycles. w_dis scatters less |
| Damping | α = 0.1 | with α = 0.02 the cycles at dx = 1.5 do not settle (section 7.2); at dx = 3 the loss changes only by −12 % between α = 0.02 and 0.1 |
| GPU | any fp32 GPU; choose RTX 3090 or RTX 5090 after the benchmark (cost per cycle) | "V100-h" is only a cost unit |

## 4. Mesh validity

| Structure | Status |
|---|---|
| walls (width ≥ L_eff ≥ 12 l_ex = 4 cells) | resolved; grid pinning ∝ exp(−π² δ/dx) ≈ 10⁻¹⁷ |
| anisotropy cubes (≥ 4 cells) | resolved; faces on cell faces |
| particle-scale stress | uniform per particle; local wall width 1/√(Q_eff (1 + r_p)) must stay ≥ 2 cells (`run_loops.py` warns) |
| vortex cores, Bloch points (l_ex scale) | NOT resolved at dx = 3 (≈ 17 … 23 neighbour pairs above 60°). Their effect depends on the particle size (below) |

**Mesh decision (planning pilot 7.8, α = 0.1, virgin state, L_eff = 12 l_ex,
mean ± SE over realisations):**

| | loss vs dx = 1.5 (9 mT / 50 mT) | Δβ | pairs > 60° |
|---|---|---|---|
| d/l_ex = 96, dx 1 | ×1.6 / ×1.22 | −0.03 ± 0.29 | 0 |
| d/l_ex = 96, dx 2 | ×0.57 / ×0.98 | +0.35 ± 0.26 | ≈ 8 |
| d/l_ex = 96, dx 3 | ×0.37 / ×0.26 | −0.22 ± 0.28 | ≈ 23 |
| d/l_ex = 300 (1 µm), dx 3 | ×1.4 / ×1.26 | −0.20 ± 0.23 | ≈ 17 (dx 1.5: ≈ 2 … 5) |

- **At d/l_ex = 300 (1 µm) dx = 3 is sufficient for β** (−0.20 ± 0.23). The
  loss level is ≈ 26 % high at 50 mT (± 5 %) and ≈ 40 % high at 9 mT
  (uncertain). dx = 3 costs 0.06 h per cycle, dx = 1.5 costs 1.87 h (V100).
- At d/l_ex = 96 (330 nm) dx = 3 is not sufficient: the loss is 3 … 4× too low.
- Cause: in a 330 nm particle a vortex carries the process. Its core (size
  l_ex) is not resolved at dx = 3; the grid holds it, it moves little at
  small fields, and its loss is missing. In a 1 µm particle many walls
  (width ≈ L_eff = 4 cells, resolved) carry the process; the few unresolved
  structures are a small part of a 30× larger volume.
- Comparison with the earlier assumptions: "walls resolved, unresolved cores
  accepted" holds for 1 µm, not for small particles. The audit prediction
  (grid-pinned cores add ≈ 10× hysteresis) is refuted: at d/l_ex = 96 the
  coarse grid gives less loss, at 1 µm only +26 %. The earlier statements
  "dx 3: Δβ +0.5, dx 2: Δβ +0.7" came from single pairs with cycle-only
  errors and are wrong.
- Use for effects: compare parameters at the same dx; then the level bias
  cancels if it is similar for all parameters. A control of one parameter
  change at dx = 3 and dx = 1.5 checks this (≈ 100 V100-h, decision by Chris).

**Convergence of the loss level (validity limit):**

- d/l_ex = 96, dx 1 vs dx 1.5: loss ×1.35 at 50 mT (Δ ln w = +0.30 ± 0.24,
  1.3 SE), ×1.45 at 9 mT (+0.37 ± 0.53); one pair of the convergence test:
  dx 1.5 ≈ 30 % low; energy of the virgin state 0.00166 vs 0.00186 K_d.
  Thus the level is not converged at dx = 1.5 for d/l_ex = 96 (a hint, not
  proven with n = 3). β does not change (−0.03 ± 0.29).
- Cause: the discretization error of the exchange energy grows
  approximately with (dx / size of the carrying structure)². At d/l_ex = 96
  a vortex core (size ≈ l_ex) carries the process: (1.5)² and (1)² are both
  large, thus even dx = 1 is coarse for the core (it needs dx ≲ 0.5 l_ex).
  At d/l_ex = 300 walls (≈ L_eff = 12 l_ex) carry it: (3/12)² ≈ 6 %,
  (1.5/12)² ≈ 1.6 % (estimate).
- Measured at d/l_ex = 300: dx 3 vs dx 1.5 +26 % at 50 mT, more than the
  6 % estimate; probably a part comes from the ≈ 17 unresolved structures at
  dx = 3. Whether dx = 1.5 itself is converged at d/l_ex = 300 is not known
  (dx = 1 costs ≈ 300 V100-h per realisation).
- **Validity:** absolute losses at 1 µm are known to ≈ 30 % (mesh). β and
  ratios between parameters at the same dx are not mesh-dependent within
  the errors of all comparisons.

Diagnostics in every run: `frac_pairs_gt60` per cycle and `n_pairs_gt60` per
sample (csv). `offset` per cycle = loop centre / amplitude: a value above
0.10 (flag `o`) is a minor loop around a remanent state, not around the
demagnetized state. The mesh test (factor-4 amplitude steps) gave offsets up
to 0.25 (L12) and 0.55 (floor) at 23 mT with dx = 1.5. `analyze.py` flags amplitudes with jumps (`#`). The control
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
   Then the numerical floor test and the initial-state test:
   `python run_queue.py jobs/numfloor.txt --gpus 0,1` and
   `python run_queue.py jobs/initproto.txt --gpus 0,1` (sections 7.3, 7.4),
   then `python run_queue.py jobs/demagtest.txt --gpus 0,1` (section 7.5).
5. Pilot: `python run_queue.py jobs/pilot.txt --gpus 0,1,2,3`. Copy `runs/`
   back before you destroy an instance.
6. `python analyze.py runs/T_*` and `python analyze.py --pair runs/T_L18 runs/T_L18_f2`.
   Report the pilot checks (section 7.7) to Chris.
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

**Result (2026-09-30, last 6 of 12 cycles):**

| Run pair | level ln(w_A/w_B) | Δβ | Note |
|---|---|---|---|
| α = 0.1: dx 3 vs 1.5 (L12) | −0.16 ± 0.04 | +0.08 ± 0.06 (PASS) | both steady: scatter 2 … 8 %, drift < 1 % per cycle |
| α = 0.02: dx 3 vs 1.5 (L12) | +0.56 ± 0.08 (w_dis) | −0.01 ± 0.11 (PASS) | dx 1.5 not steady: scatter of w_loop up to 92 % |
| α = 0.02: dx 3 vs 1.5 (floor) | +0.63 ± 0.06 | −0.02 ± 0.08 (PASS) | |
| dx 3: α 0.1 vs 0.02 | −0.12 ± 0.05 | +0.12 ± 0.07 | at dx = 3 the loss hardly depends on α |
| dx 1.5: α 0.1 vs 0.02 | +0.59 ± 0.07 | +0.04 ± 0.10 | at dx = 1.5 and α = 0.02 the cycles do not settle |

Conclusions:
- The slope β does not depend on the mesh (Δβ ≤ 0.1 in all pairs). The
  predicted 10× grid pinning of the cores does not occur.
- With α = 0.1 the model is steady and mesh-independent (level −16 %, Δβ
  +0.08). With α = 0.02 the fine mesh does not settle within 12 cycles and
  gives about half the loss.
- Thus α = 0.1 is the better choice: steady cycles (scatter below the 10 %
  target) and mesh independence. The damping loss at dx = 3 is small
  (α 0.02 → 0.1 changes w by −12 %).

### 7.3 Numerical floor test (d/l_ex = 96, dx = 3, virgin state, ascending, ≈ 4.2 V100-h)

The loops in the study range are small: at 23 mT the mesh test gave w ≈ 10⁻⁷
… 10⁻⁶ K_d per cycle. At 10 mT, w is smaller again. Question: do the solver
tolerance and fp32 make a part of this loss?

| Run | precision | atol |
|---|---|---|
| F_si_atol1e-05 (reference of the study) | fp32 | 10⁻⁵ |
| F_si_atol1e-06 | fp32 | 10⁻⁶ |
| F_do_atol1e-05 | fp64 | 10⁻⁵ |
| F_do_atol1e-06 (best) | fp64 | 10⁻⁶ |

L_eff = 12 l_ex, virgin state (protocol A), B_peak = 9 / 50 / 150 mT
ascending, 7 cycles each. Evaluation:

    python analyze.py --skip 2 --cycles runs/F_*
    python analyze.py --skip 2 --compare runs/F_si_atol1e-05 runs/F_do_atol1e-06
    python analyze.py --skip 2 --w dis --compare runs/F_si_atol1e-05 runs/F_do_atol1e-06

PASS: at 9, 50 and 150 mT the mean w of the reference is within ±15 % of the
best run (90 % interval), for w_loop and for w_dis. If it fails: the smallest
setting that passes sets atol and precision for the pilot (cost to Chris).

**Result (2026-09-30): PASS.** fp32 / atol 10⁻⁵ vs fp64 / atol 10⁻⁶: level
−0.012 ± 0.036, Δβ −0.012 ± 0.032. At 9 and 50 mT all four settings agree to
better than 1 %. Energy balance w_dis/w_loop within 2 %. fp32 and atol 10⁻⁵
stay. (fp64 on the small grid cost only ≈ 0.3 h, not the modelled 1.7 h.)

### 7.4 Initial-state test (d/l_ex = 96, L_eff = 12, ≈ 2.6 V100-h)

Question: does the loss in the study range depend on how the demagnetized
state is made? B_peak = 9 / 50 / 150 mT, 7 cycles each.

| Run | Protocol | dx |
|---|---|---|
| F_si_atol1e-05 (from 7.3) | A: virgin state, ascending | 3 |
| I_A_virgin_dx15 | A: virgin state, ascending | 1.5 |
| I_B_acdemag_dx3 | B: 17 decaying saturating cycles (−20 % per half cycle), ascending | 3 |
| I_C_reset_dx3 | C: 1 saturating cycle, descending (present) | 3 |

Evaluation (`--skip 2`, both `--w loop` and `--w dis`):
- `--compare` A (dx 3) vs B: PASS (level ±0.15, slope ±0.2) → the demagnetized
  state does not matter; use A (cheaper, no demagnetization cycles).
- A vs C: expected to differ at small B_peak (offset flag `o` in C).
- A dx = 3 vs A dx = 1.5: the virgin state from a random start contains many
  vortices and Bloch points; the grid can hold some of them. Compare
  `init.n_pairs_gt60` and `init.E_over_Kd` in summary.json and the loss.
- Offsets: `o` flags in A and B must be absent.

**Result (2026-09-30, α = 0.02, `--skip 2`):**

| Pair | level | Δβ | Offset at 9 / 50 / 150 mT |
|---|---|---|---|
| A dx 3 vs B dx 3 | +0.12 ± 0.04 | −0.10 ± 0.03 | A: −0.53 / −0.09 / −0.04; B: +0.11 / +0.02 / 0.00 |
| A dx 3 vs C dx 3 | −0.37 ± 0.07 | −0.43 ± 0.06 | C: −0.80 / −0.15 / −0.05 |
| A dx 3 vs A dx 1.5 | −0.61 ± 0.07 | +0.46 ± 0.06 | A dx 1.5: +0.03 / +0.02 / 0.00 |

Conclusions:
- C (present protocol) is out: minor loops around a remanent state (offset
  up to −0.80), Δβ = −0.4.
- A at dx = 3 is not a clean virgin state: the relaxation (5000 iterations,
  not converged) leaves 17 unresolved large-angle pairs and a net M of
  ≈ −0.003 Ms. This is comparable to the 9 mT loop amplitude (offset −0.53).
  At dx = 1.5 the virgin state has no such pairs and M ≈ 0. Thus A depends
  on the mesh through the initial state.
- B (AC demagnetization) is steady (scatter 0.7 … 11 %) and nearly centred.
- The tests ran with α = 0.02. After the steady test, α = 0.1 is preferred.
  Open: A and B with α = 0.1 at dx = 3 and 1.5 (section 12).

### 7.5 Demagnetization test (d/l_ex = 96, L_eff = 12, α = 0.1, ≈ 5.7 V100-h)

Question: gives the AC demagnetization (B) the same state and loss as a clean
virgin state, and can the demagnetizing cycles run at 4 f (only the end
state counts; 17 cycles at f cost ≈ +35 % per production run, at 4 f ≈ +8 %)?

| Run | Protocol | dx |
|---|---|---|
| D_B_dx3 / dx15 | B at f | 3 / 1.5 |
| D_B4f_dx3 / dx15 | B, demagnetizing cycles at 4 f (`--demag_f_factor 4`) | 3 / 1.5 |
| D_A_dx15 | A (virgin state; clean at dx = 1.5) | 1.5 |

Evaluation (`--skip 2`, `--compare`, tolerances as in 7.1, offsets):
B4f vs B (both meshes), B dx 3 vs B dx 1.5, B dx 1.5 vs A dx 1.5. If B4f
passes against B and against A: production uses B4f with α = 0.1 and dx = 3.

**Result (2026-09-30, α = 0.1, `--skip 2`):**

| Run | w at 9 / 50 / 150 mT | offset at 9 / 50 / 150 mT |
|---|---|---|
| B dx 3 | 5.4·10⁻⁸ / 5.1·10⁻⁶ / 1.8·10⁻⁴ | −0.59 / −0.11 / −0.04 |
| B4f dx 3 | 5.4·10⁻⁸ / 2.0·10⁻⁶ / 1.9·10⁻⁴ | −1.26 / −0.23 / −0.08 |
| B dx 1.5 | 1.6·10⁻⁷ / 6.5·10⁻⁶ / 1.4·10⁻⁴ | +0.71 / +0.13 / +0.05 |
| B4f dx 1.5 | 1.4·10⁻⁷ / 7.8·10⁻⁶ / 1.1·10⁻⁴ | +0.68 / +0.12 / +0.05 |
| A dx 1.5 | 1.5·10⁻⁷ / 4.8·10⁻⁶ / 1.0·10⁻⁴ | +0.02 / +0.01 / +0.01 |

| Pair | level | Δβ |
|---|---|---|
| B dx 3 vs B dx 1.5 | −0.37 ± 0.04 | **+0.46 ± 0.03 (FAIL)** |
| B4f dx 3 vs B4f dx 1.5 | −0.56 ± 0.07 | **+0.47 ± 0.06 (FAIL)** |
| B4f vs B, dx 3 / dx 1.5 | −0.28 / −0.09 | −0.02 / −0.03 (PASS) |
| B dx 1.5 vs A dx 1.5 | +0.21 ± 0.04 | +0.11 ± 0.03 (PASS) |

Conclusions:
- **In the study range the mesh dx = 3 is probably not sufficient.** At 9 mT
  dx = 3 gives ≈ 3× less loss than dx = 1.5 (2 pairs, same seed); β is
  ≈ 0.5 steeper (see the correction in 7.6: paired errors). The steady
  test (90 and 375 mT) did not show this: the problem is at small B_peak.
  At small fields the coarse grid holds structures in place that move on
  the fine grid (this agrees with the mesh test: at 23 mT dx = 1.5 dissipated
  5 … 10× more).
- The AC demagnetization of 4 particles ends with a net M of ≈ 0.35 % Ms
  (offset ±0.6 … 1.3 relative to the 9 mT loop). The virgin state A at
  dx = 1.5 has M ≈ 0.02 % Ms (offsets ≤ 0.02).
- Correction (Chris): this net M does not matter for the loss. At 9 mT, where
  the offset is largest, B and A give the same w (ratio 1.03). At 50 and
  150 mT, where the offset is small, B is 35 % higher. Thus A and B differ by
  the domain configuration of the demagnetized state, not by the net M. B is
  not out. Open: is this 35 % a real difference or the scatter between
  initial states (C_dx15_is2 in 7.6 answers this)?
- The AC demagnetization does not remove the structures that the coarse grid
  holds (dx = 3: ≈ 21 large-angle pairs after the demagnetization, as after
  the random start). At dx = 1.5 neither A nor B has such structures.
- Open: is dx = 1.5 converged at 9 mT? And the cost: a run at d/l_ex = 300
  costs 6 V100-h at dx = 3, 26 at dx = 2 and 72 at dx = 1.5 (section 12).

### 7.6 Mesh convergence test (d/l_ex = 96, L_eff = 12, α = 0.1, ≈ 7.7 V100-h)

Virgin state A with the same start on every mesh (`--init_block 6`: random m
on 6 l_ex blocks in physical coordinates), B_peak = 9 / 20 / 50 mT ascending,
7 cycles.

| Run | dx/l_ex | Purpose |
|---|---|---|
| C_dx2 / C_dx15 / C_dx1 | 2 / 1.5 / 1 | convergence of w and β at small B_peak |
| C_dx15_is2 | 1.5 | second initial state (`--init_seed 2`, same cubes): scatter of the virgin state |
| C_dx15_direct50 | 1.5 | 50 mT directly from the virgin state: if it equals the 50 mT stage of C_dx15, every amplitude can run as its own job (parallel production) |

Evaluation (`--skip 2`): `--compare` C_dx2 and C_dx15 against C_dx1 (level
±0.15, Δβ ±0.2); C_dx15 vs C_dx15_is2 gives the scale of the initial-state
scatter; the 50 mT w of C_dx15_direct50 vs C_dx15 (±15 %).
Also decide A vs B: if the difference A (D_A_dx15) vs B (D_B_dx15) of 7.5
(level +0.21, Δβ +0.11) lies within the initial-state scatter, A and B are
equivalent (then A: no demagnetizing cycles; B: the measurement standard).

**Result (2026-09-30):**

Same initial state (init_block 6, seed 1) on three meshes:

| Run | w 9 mT | w 20 mT | w 50 mT | β (9 → 50 mT) | virgin state E/K_d, n60 |
|---|---|---|---|---|---|
| C_dx2 | 6.4·10⁻⁸ | 3.3·10⁻⁷ | 5.0·10⁻⁶ | 2.54 | 0.00207, 3 |
| C_dx15 | 3.4·10⁻⁷ | 9.3·10⁻⁷ | 7.0·10⁻⁶ | 1.77 | 0.00186, 0 |
| C_dx1 | 4.4·10⁻⁷ | 1.3·10⁻⁶ | 1.0·10⁻⁵ | 1.83 | 0.00166, 0 |

- dx = 1.5 vs dx = 1: Δβ = −0.05, level −0.33 (w ≈ 30 % low). One pair only.
- dx = 2 vs dx = 1: Δβ = +0.74, level ×0.26. One pair only.
- **Correction:** the errors of `--compare` contain only the cycle
  statistics. Two independent realisations differ by σ_β ≈ 0.6 (use
  `--real`). Pairs with the same start and cubes differ much less: over the
  3 pairs coarse (dx 2 or 3) vs dx 1.5 the scatter of Δβ is 0.18, and
  `--paired` gives Δβ = +0.60 ± 0.10 (95 %: +0.16 … +1.03). Thus the coarse
  mesh makes β steeper with some confidence; dx 1.5 vs dx 1 is not decided
  (1 pair). The level trend at the smallest B_peak (coarser mesh → less
  loss) has the same sign in all ≈ 5 pairs. The planning pilot (7.8)
  measures this per dx with 3 pairs each.
- 50 mT directly from the virgin state vs after 9 and 20 mT: −5 %. Each
  amplitude can run as its own job (parallel production).

Realisation scatter at dx = 1.5 (8 runs: other initial states and/or other
cubes, protocols A and B):

| | 9 mT | 50 mT | β (9 → 50 mT) |
|---|---|---|---|
| mean | w = 1.6·10⁻⁷ | w = 8.3·10⁻⁶ | 2.29 |
| scatter per realisation | σ(ln w) = 0.41 | σ(ln w) = 0.57 | σ = 0.45 |

- The difference A vs B (level +0.21, Δβ +0.11) is well inside this
  scatter: A and B are equivalent. A is used (no demagnetizing cycles).
- **One realisation of 4 particles is not representative.** For an error
  of ±0.1 in β a parameter point needs ≈ 20 realisations; for ±0.15 ≈ 9.
  This scatter is measured at d/l_ex = 96. At d/l_ex = 300 each particle has
  more domains, and the scatter may be smaller (not measured).

### 7.8 Planning pilot (≈ 105 V100-h, ≈ 13 $, 10 GPUs, ≈ 12 h)

Protocol of all runs: virgin state (random m on 6 l_ex blocks), α = 0.1,
φ = 0.65, fp32, atol 10⁻⁵, 7 cycles, B_peak = 9 and 50 mT, β = β(9 → 50 mT).
The same seed gives the same cubes and the same start on every mesh (pairs).

| Block | Runs | Purpose |
|---|---|---|
| 1 | P1_d150 (dx 1.5, seeds 1–4); P1_d300 (dx 1.5, seeds 1–3, one job per amplitude); P1_d300_dx3 (dx 3, seeds 1–3, one job per amplitude) | realisation scatter vs particle size (d/l_ex = 96 / 150 / 300); dx 3 vs 1.5 at 1 µm (paired) |
| 2 | P2_L12 (seeds 6–8); P2_L30 (seeds 3–8); P2_rp3 (seeds 3–8), d/l_ex = 96, dx 1.5 | effect sizes: Q_eff (L_eff 12 vs 30, not paired) and r_p = 3 (paired with L12 of the same seed) |
| 3 | P3_dx1, P3_dx2, P3_dx3 (seeds 3–5), d/l_ex = 96 | mesh, paired with R_dx15_s3…s5 |

Evaluation: σ_β per d from the realisations; `--paired` for r_p, the mesh
pairs and dx 3 at d/l_ex = 300; the Q_eff effect with `--real` or from the
ensemble means. Result: the number of realisations, d/l_ex and dx of the
production (section 9).

**Result (2026-10-02, 40 runs; f = 30 MHz; P = loss density of the unit cell
including the voids, P = W · f; mean ± σ per realisation (SE of the mean)):**

| Group | n | P at 9 mT [W/cm³] | P at 50 mT [W/cm³] | β (9 → 50 mT) |
|---|---|---|---|---|
| d/l_ex 300, dx 1.5 | 3 | 5.2 ± 4.2 (2.4) | 138 ± 14 (8) | 2.06 ± 0.39 (0.23) |
| d/l_ex 300, dx 3 | 3 | 7.3 ± 1.5 (0.9) | 174 ± 9 (5) | 1.86 ± 0.11 (0.06) |
| d/l_ex 150, dx 1.5 | 4 | 5.5 ± 2.6 (1.3) | 303 ± 50 (25) | 2.37 ± 0.31 (0.15) |
| d/l_ex 96, dx 1.5 (base) | 8 | 4.7 ± 2.6 (0.9) | 333 ± 163 (57) | 2.47 ± 0.57 (0.20) |
| d/l_ex 96, dx 1 | 3 | 7.5 ± 6.4 (3.7) | 405 ± 105 (61) | 2.44 ± 0.37 (0.22) |
| d/l_ex 96, dx 2 | 3 | 2.7 ± 0.7 (0.4) | 326 ± 81 (47) | 2.82 ± 0.30 (0.17) |
| d/l_ex 96, dx 3 | 3 | 1.8 ± 0.6 (0.4) | 85 ± 35 (20) | 2.25 ± 0.35 (0.20) |
| d/l_ex 96, r_p = 3 | 6 | 1.9 ± 0.8 (0.3) | 79 ± 23 (9) | 2.16 ± 0.20 (0.08) |
| d/l_ex 96, L_eff = 30 | 6 | 5.6 ± 3.6 (1.5) | 104 ± 6 (2) | 1.84 ± 0.45 (0.18) |

- σ is the scatter of one realisation (4 particles). A real core averages
  over very many particles; the mean ± SE is the estimate for a real sample.
  σ sets the number of realisations: n ≈ (σ/SE_target)².
- At 50 mT the scatter falls with the particle size (σ(ln w): 0.56 → 0.18 →
  0.10); at 9 mT it does not (0.56 / 0.46 / 0.77): few single events carry
  the small loss. The scatter of β comes from the 9 mT point.
- Pairs with the same cubes and the same start do not reduce the scatter
  (the relaxation ends in other states on another mesh or with r_p).
  Evaluation therefore uses group means with Welch errors.
- Hints for β (d/l_ex = 96, vs base): L_eff 30: Δβ −0.63 ± 0.27, loss ×0.31;
  r_p = 3: Δβ −0.31 ± 0.22, loss ×0.24. L_eff = 30 at d/l_ex = 96 has
  d/L_eff = 3.2 (< 5, outside the valid window): measure this effect at 1 µm.
- β ≈ 2 between 9 and 50 mT. The dynamic (damping) share at 30 MHz and
  α = 0.1 is not measured; it also gives β = 2. Check with f vs 2 f first.

### 7.9 Frequency test (PLAN step 1, 2026-10-02/03): NEGATIVE

d/l_ex = 300, dx = 3, α = 0.1, virgin state, 7 amplitudes; f (30 MHz) and
2 f, 4 seeds each; base point (L_eff 12, r_p 0) and corner (L_eff 30, r_p 3).
w_h = 2 w(f) − w(2 f), dynamic share = 1 − w_h/w(f), mean ± SE over the seeds.

| B_peak | dynamic share, base | dynamic share, corner | w(2 f)/w(f) base / corner |
|---|---|---|---|
| 9 mT | 0.11 ± 0.11 | 0.31 ± 0.32 | 1.06 / 1.21 |
| 20 mT | 0.03 ± 0.09 | −0.16 ± 0.10 | 1.02 / 0.82 |
| 35 mT | 0.43 ± 0.05 | 0.33 ± 0.14 | 1.44 / 1.31 |
| 50 mT | 0.83 ± 0.17 | 0.64 ± 0.13 | 1.83 / 1.62 |
| 100 mT | 1.06 ± 0.07 | 1.19 ± 0.06 | 2.06 / 2.19 |
| 150 mT | 1.09 ± 0.07 | 1.24 ± 0.04 | 2.09 / 2.24 |

- At 9 … 20 mT the loss per cycle does not depend on f: hysteresis.
- From 50 mT the loss per cycle is ∝ f: almost all of it is damping loss
  (∝ α f). w_h is a small difference of large numbers there and is not
  measurable at 100 and 150 mT (negative values).
- Rule of PLAN step 1: FAIL at both points.
- **Correction of earlier statements:** β ≈ 2 between 9 and 50 mT (planning
  pilot 7.8) and the effects of Q_eff, r_p and d/l_ex on the loss at 50 mT
  describe mainly the damping loss at 30 MHz, not the hysteresis. Only the
  values at 9 … 20 mT are hysteresis. The LLG protocol at 30 MHz cannot give
  the hysteresis β in the study range above ≈ 30 mT.

### 7.10 Mesh error at d/l_ex = 212 (PLAN step 1b, 2026-10-03): rule not met

d/l_ex = 212, dx 3 vs dx 1.5, 4 seeds each, 9 and 50 mT (Welch):

| | Δ ln w 9 mT | Δ ln w 50 mT | Δβ (9 → 50 mT) |
|---|---|---|---|
| d/l_ex 212: dx 3 − dx 1.5 | −0.35 ± 0.22 | −0.06 ± 0.14 | +0.14 ± 0.18 |
| d/l_ex 300: dx 3 − dx 1.5 | +0.54 ± 0.46 | +0.23 ± 0.07 | −0.20 ± 0.23 |
| difference 212 − 300 | −0.89 ± 0.51 | −0.29 ± 0.15 | +0.34 ± 0.30 |

Rule (|difference at 50 mT| + 1.645 SE < 0.10): 0.54, not met. The mesh
error is not shown to be the same at both ends of the d/l_ex range. Note:
by 7.9 the 50 mT values are mainly damping loss; the mesh question must be
asked again for the hysteresis once its method is fixed.

### 7.11 Surface roughness (PLAN step 1c, 2026-10-03): no valid result (7.13)

**Withdrawn after the audit (7.13):** all 9 mT values and all values at
d/l_ex 300 are disturbed (not converged relaxation; at d/l_ex 300 each
amplitude ran as its own job from the virgin state); the smooth variant at
d/l_ex 212 and 50 mT is doubtful (w_dis/w_loop up to 2.15, 2 of 8 amplitudes
missed). Only the staircase reference at d/l_ex 212 and 50 mT is clean
(−6 % ± 14 %). Thus there is no valid comparison of the surfaces. The
conclusion "the staircase is not the main cause" is withdrawn. All runs use
the staircase surface (`--surface voxel`, default). Next: if mesh300 (PLAN
3e) shows a mesh error at dx 3, repeat this test at d/l_ex 300 with the
converged relaxation (≈ 10 V100-h); else the surface is not relevant.
The original evaluation follows.

dx = 3 with a volume-fraction surface (`--surface fraction`, only `surface`
differs from the references), mesh error vs dx = 1.5 (Welch, mean ± SE):

| | Δ ln w 9 mT | Δ ln w 50 mT | Δβ (9 → 50 mT) |
|---|---|---|---|
| d/l_ex 212, dx 3 staircase | −0.35 ± 0.22 | −0.06 ± 0.14 | +0.14 ± 0.18 |
| d/l_ex 212, dx 3 smooth | +0.01 ± 0.26 | −0.20 ± 0.09 | −0.06 ± 0.14 |
| d/l_ex 300, dx 3 staircase | +0.54 ± 0.46 | +0.23 ± 0.07 | −0.20 ± 0.23 |
| d/l_ex 300, dx 3 smooth | +0.31 ± 0.45 | +0.31 ± 0.11 | −0.01 ± 0.23 |

- The smooth surface does not remove the level error at 50 mT; the
  difference between the two ends grows (−0.29 → −0.51). The staircase is
  not the main cause of the level error. The β error is slightly smaller
  with the smooth surface, all within ≈ 1 SE. No model change.
- **The d/l_ex = 300 references are disturbed at 50 mT:** w_dis/w_loop =
  1.44 … 2.54 (dx 1.5) and 1.35 … 1.52 (smooth) because the 50 mT job starts
  directly from the not converged virgin state (same cause as at 9 mT,
  PLAN step 1d). At d/l_ex = 212 (9 → 50 mT in sequence) the balance is
  ≈ 1.0. Thus all mesh comparisons at d/l_ex = 300 and 50 mT (7.8: +26 %,
  7.10, this test) are not reliable; only d/l_ex = 212 is clean: mesh error
  at 50 mT −6 % ± 14 %, β +0.14 ± 0.18.
- Next: the relaxation test (PLAN step 1d); then a repeat of the d/l_ex =
  300 comparison with the converged relaxation (decision by Chris, ≈ 10 $).

### 7.12 Damping α (2026-10-03, base point, 1 µm, dx 3, 3 … 4 seeds each)

**Partly withdrawn (7.17):** with the converged virgin state and n = 6 vs 11, α
changes β (−0.27 ± 0.08 for α 0.01 vs 0.1). The statement "α acts on the
level, not on β" below is not valid.

Loss ∝ α^k with k = 0.12 ± 0.05 (20 mT), 0.17 ± 0.02 (50 mT), 0.29 ± 0.04
(100 mT), 0.25 ± 0.09 (150 mT): a factor 1.3 … 1.9 per decade of α, far from
∝ α. Fit per realisation over 20 … 150 mT (9 mT excluded): a power law fits
(χ²/dof ≤ 0.7, no significant curvature); β = 1.52 ± 0.18 (α 0.01),
1.53 ± 0.09 (0.03), 1.69 ± 0.04 (0.1); dβ/d log₁₀α = +0.17 ± 0.14 (1.3 SE).
**α acts on the loss level, not on β** (within the errors). A two-point
secant 50 → 100 mT had shown +0.38 ± 0.08 per decade; the fit over all
amplitudes does not confirm it (two-point values depend on the noise of
two points). Small α: the same cost per cycle, but 3 … 7× more scatter from
cycle to cycle. Pending: 3 more seeds at α = 0.01 (904 … 906).
Decomposition: w ∝ B̂² sin δ / μ; μ = 8.0 at all B̂ and α (geometry limit of
the sphere packing); β < 2 because the loss angle sin δ falls with B̂
(∝ B̂^−0.3 … −0.5). No Rayleigh rise in 20 … 150 mT. Physical or a model
effect: open (domain snapshots, PLAN step "snaptest"; reference data needed).

### 7.13 Audit of all runs (2026-10-03): the virgin state was never converged

The virgin-state protocol requires a fully relaxed state. In **no** run did
the relaxation converge (`relax_converged = False`, 5000 iterations; at
d/l_ex = 300 and dx = 3 about 25 000 … 30 000 are necessary). The work went
on although the code reported it; this was an error of the procedure.
Table: results/audit_2026-10-03.txt (per group: convergence, energy balance
of the first and the later stages, closure, offset, amplitude hit).

| Runs | First measured stage | Later stages |
|---|---|---|
| d/l_ex = 96 (sections 7.3 … 7.6, 7.8 block 2/3) | balance 1.0 … 1.5 (dx = 1: 2.4) | ≈ 1.0 |
| d/l_ex 150 … 300, amplitudes in sequence (7.8 d150, 7.9, 7.10, 7.12) | **9 mT disturbed** (balance up to 64) | ≈ 1.0 (exception: seed 902 at 50 mT for α ≤ 0.03, an event) |
| d/l_ex = 300, one amplitude per job (7.8 P1_d300*, 7.11 S300*) | **every amplitude disturbed** (up to 76) | — |

**Not reliable (withdrawn until repeated):**
- all 9 mT values and β(9 → 50 mT) at d/l_ex ≥ 150;
- the scatter vs particle size at 9 mT (σ = 0.77 at 1 µm), which entered
  the planning of the number of realisations (PLAN D1);
- all mesh comparisons at d/l_ex = 300 (7.8: +26 %; 7.10; 7.11);
- the dynamic share at 9 mT (7.9); α at 9 mT (7.12);
- dx = 1 at d/l_ex = 96 at 9 mT (balance 2.4).

**Still valid:** results at 20 … 150 mT from runs with amplitudes in
sequence (7.9, 7.10 at 50 mT for d/l_ex = 212, 7.12: β ≈ 1.5 … 1.7, μ = 8,
loss ∝ α^0.1 … 0.3) and the results at d/l_ex = 96 (balance ≈ 1).
Offsets above 0.1 at small B_peak are frequent in the virgin runs at
d/l_ex = 96 (residual net M of the unconverged state); by 7.5 they change
the loss little, but they are checked again after the fix.

Fix: `--relax_mode converge` (PLAN step 1d) and a hard gate (a run stops if
the relaxation does not converge; analysis excludes such runs). Repeats:
PLAN section "repeats".

### 7.15 Relaxation to convergence (PLAN step 1d, 2026-10-03): fix confirmed

d/l_ex 300, dx 3, base point, `--relax_mode converge`, seeds 901 … 903,
cycles from 2 on (mean ± SE over the seeds, n = 3 per α):

| | α = 0.1 | α = 0.01 |
|---|---|---|
| relaxation converged | 3/3 (26 000 … 45 000 iterations, 0.5 … 0.8 h) | 3/3 |
| w_dis/w_loop at 9 mT (before: up to 23 … 64) | 1.08 / 1.02 / 1.00 | 1.11 / 1.07 / 1.01 |
| w_dis/w_loop at 20 mT | 0.99 / 1.01 / 1.00 | 1.01 / 0.90 / 0.96 |
| β(9 → 20 mT) | 1.81 ± 0.09 | 1.63 ± 0.12 |

- Rule (balance 0.9 … 1.1 at 9 mT): met for α 0.1 (3/3); α 0.01 1 of 3
  just above (1.11).
- α 0.01 vs 0.1 (same seeds and relaxed states): Δβ(9 → 20) = −0.17 ± 0.15
  (Welch), Δ ln w = −0.07 ± 0.07 at 9 mT and −0.21 ± 0.12 at 20 mT (paired).
  No effect on β shown; α acts on the level (as 7.12).
- At α 0.01 the 9 mT loop still jumps between cycles (up to 85 % from cycle to
  cycle); at α 0.1 about 6 %.
- Seed 902 needed 45 000 iterations: `--relax_maxiter` is 80 000 for all
  repeat jobs (also at dx 1.5).

### 7.16 β(B_peak) with a valid 9 mT value (2026-10-03)

`eval_relaxfix.py` → results/relaxfix_beta.txt. 9 mT from the relaxation test
(converged), 20 … 150 mT from the runs with the same seeds and settings
(F_base_f_*, A010_*), whose later stages are valid (7.13). Splice check at
20 mT (ln w new − old, same seed): +0.02 ± 0.03 (α 0.1), +0.06 ± 0.04
(α 0.01): the two sources agree. d/l_ex 300, dx 3, L_eff 12 l_ex, r_p 0,
f = 30 MHz; mean ± SE over n = 3 realisations.

| | α = 0.1 | α = 0.01 |
|---|---|---|
| β power law 9 … 150 mT | 1.74 ± 0.03 | 1.55 ± 0.07 |
| curvature dβ/d ln B | −0.07 ± 0.09 | −0.03 ± 0.17 |
| β at 10 / 50 / 100 mT (quadratic fit) | 1.84 ± 0.13 / 1.72 ± 0.03 / 1.67 ± 0.07 | 1.57 ± 0.13 / 1.53 ± 0.15 / 1.51 ± 0.26 |
| P at 10 / 50 / 100 mT [W/cm³] | 11.5 / 201 / 650 | 10.4 / 127 / 364 |
| β per segment 9→20 / 20→35 / 35→50 / 50→70 / 70→100 / 100→150 mT | 1.78 / 1.82 / 1.64 / 1.96 / 1.43 / 1.74 (SE 0.02 … 0.25) | 9→20 1.57, 20→50 1.60, 50→100 1.28, 100→150 1.95 |

- β is constant within the errors over 10 … 150 mT (no curvature); the
  segments scatter more than the fits (short lever arms, SE up to 0.25).
- The earlier values with 9 mT (β(9 → 50) = 1.86 at dx 3, 2.06 at dx 1.5, 1 µm,
  planning pilot) were too high by
  the unconverged state.
- α 0.01 − 0.1: Δβ (power law) = −0.19 ± 0.08 (2.4 SE, Welch, n = 3 + 3);
  with n = 3 (t₉₅ ≈ 2.8) this is a hint, not a result. alphatest2 (3 more
  realisations at α 0.01, converged) decides it.

### 7.17 Baseline at 1 µm with the converged virgin state (2026-10-04)

`eval_baseline.py` → results/baseline_eval.txt. d/l_ex 300, dx 3, base point
(L_eff 12, r_p 0, φ 0.65), α 0.1, 7 amplitudes 9 … 150 mT, seeds 921 … 928,
`--relax_mode converge` (all 8 converged, 16 000 … 43 000 iterations).
Amplitudes hit within 0.7 %, closure ≤ 0.5 %, offset ≤ 0.12.

| | mean | sd (one run) | SE (n = 8) | without s921 (n = 7): mean / sd |
|---|---|---|---|---|
| β power law 9 … 150 mT | 1.76 | 0.15 | 0.05 | 1.71 / 0.04 |
| β at 10 mT | 1.93 | 0.26 | 0.09 | 1.86 / 0.17 |
| β at 50 mT | 1.72 | 0.14 | 0.05 | 1.67 / 0.06 |
| β at 100 mT | 1.62 | 0.16 | 0.06 | 1.59 / 0.13 |
| P at 10 / 50 / 100 mT [W/cm³] | 10 / 189 / 603 | scatter factor ×/ 1.50 / 1.19 / 1.11 | | |

- Energy balance w_dis/w_loop outside 0.9 … 1.1 in 5 of 56 stages:
  s921 at 9 mT (1.40) and 100 mT (1.83), s922 at 20 mT (2.07), s925 at
  35 mT (1.21), s926 at 35 mT (1.27). Cause: single events in the kept
  cycles. s921, 100 mT, cycle 4: w_dis = 6.3·10⁻⁵ vs w_loop = 1.6·10⁻⁵ and
  the number of large-angle pairs falls from 28 to 23 (an irreversible
  rearrangement of the domain structure). s921 at 9 mT and s922 at 20 mT:
  w_dis ≫ w_loop in cycles 0 … 3 (the converged virgin state still releases
  energy in the first cycles of the AC field; 2 skipped cycles are not
  always enough).
- By PLAN 4.3a the flagged runs stay in; sensitivity: without s921 the
  scatter of β (power law) falls from 0.15 to 0.04. The realisation scatter
  is dominated by rare events, not by the cube landscape.
- The scatter is much smaller than assumed in PLAN 4.3 (σ_β ≈ 0.45 at
  d/l_ex 96): σ = 0.14 … 0.26 (with s921). Scaling of the review Monte Carlo
  (SE = 0.31 σ √(136/N)): N = 26 … 86 points for SE ≤ 0.1; at N = 64:
  SE = 0.06 … 0.12 (β at 10 mT is the largest).
- **Evaluation changed (2026-10-04):** β is now the power law over 3
  neighbouring amplitudes (`analyze.features`, PLAN 4.3); the quadratic fit
  in ln-ln is dropped. Baseline (n = 8): β_pow 1.76 ± 0.05 (sd 0.15),
  β 9-20-35 mT 1.90 ± 0.10 (sd 0.27), 35-50-70 mT 1.58 ± 0.05 (sd 0.15),
  70-100-150 mT 1.77 ± 0.07 (sd 0.19). All valid runs: results/beta_windows.md,
  results/beta_all_valid.md. The numbers at 10 / 50 / 100 mT in this section
  and in 7.16 come from the quadratic fit (kept for the record).
- α 0.01 − 0.1 (Welch): clean runs (A010c n = 3 vs baseline n = 8):
  β power law −0.34 ± 0.14 (2.5 SE), at 10 / 50 / 100 mT −0.58 ± 0.22 /
  −0.28 ± 0.12 / −0.15 ± 0.10; pooled with the spliced runs (6 vs 11):
  −0.27 ± 0.08 (3.3 SE), −0.44 ± 0.14 / −0.23 ± 0.09 / −0.14 ± 0.13.
  **α changes β** (smaller damping → smaller β, most at 10 mT). The earlier
  statement "α acts on the level, not on β" (7.12, from 20 … 150 mT and
  n = 3) is withdrawn. α is a model parameter of the loss at 30 MHz, not
  only a numerical choice.

### 7.18 Wall measures at 1 µm: snapshots at 20 and 150 mT (2026-10-04)

`domains_batch.py` on V_s901, V_s902 (base point, α 0.1, converged virgin
state; 4 phases of the last cycle; analysed on the instance). Totals over the
4 particles (V = 5.69·10⁷ l_ex³; one 180° disc through a particle centre
πR² = 7.07·10⁴ l_ex²):

| | virgin state | 20 mT | 150 mT |
|---|---|---|---|
| A_w = T(2 L_eff)/0.85 [10⁵ l_ex²] | 3.35 / 3.55 | 3.34 / 3.55 | 3.35 / 3.36 |
| A_F90 [10⁴ l_ex²] | 4.5 / 4.8 | 4.5 … 4.6 / 4.7 … 4.9 | 4.7 … 5.4 / 5.1 … 5.4 |
| F90 | 2.3 / 2.4 % | 2.2 … 2.4 % | 2.3 … 2.7 % |
| l_C [l_ex] | 111 / 105 | 111 / 105 | 108 |
| switched volume per half cycle (H max ↔ H min) | | 2.0 / 2.1 % of V | 12.8 / 12.5 % of V |

(values for s901 / s902)

- **No clear 180° walls.** F90 = 2.2 … 2.7 % is below the value of one
  synthetic vortex at this size (4.3 %) and far below one 180° wall per
  particle (12.3 %). A_w/A_F90 ≈ 6 … 7: most of the texture T is a
  continuous rotation of m (curling, vortex-like), not walls. The texture
  corresponds to ≈ 1.2 disc areas per particle.
- **The texture does not change with B_peak** (20 → 150 mT: ≤ 1 %; s902
  loses 6 % after 150 mT, an irreversible rearrangement as in 7.17).
- **The switched volume grows almost linearly with B_peak** (half cycle
  2 % → 12.6 %, B ratio 7.5, exponent ≈ 0.9). With μ ≈ 8 constant this is
  the reversible rotation of m that carries ΔM ∝ B, not the motion of walls.
- Consequence for the hypothesis "β depends on the total wall area": in
  this model at 1 µm there are no distinct walls, and the amount of texture
  is constant over the amplitude range. β < 2 does not come from walls
  that appear or vanish with B_peak; it comes from the loss per unit of
  reversible rotation (damping dynamics, consistent with the α effect in
  7.17). Across realisations (baseline, n = 8) and the Sobol points the
  texture can still correlate with β; this is tested next.
- Short-scale ripple (against 2 L_eff): 6 … 8°.

### 7.19 Wall measures of the baseline vs B_peak and vs β (2026-10-04)

`eval_walls.py` → results/walls_eval.md (domains.json of B300_s921 … s928,
4 phases of the last cycle per amplitude; mean ± SE over the 8 seeds).

| B_peak | texture A_w/V [10⁻³/l_ex] | 180° discs per particle | F90 | l_C [l_ex] | switched volume per half cycle | wall path s [l_ex] |
|---|---|---|---|---|---|---|
| virgin | 5.93 ± 0.03 | 1.19 | 2.7 % | 110 | – | – |
| 9 mT | 5.92 ± 0.03 | 1.19 | 2.7 % | 110 | 1.06 ± 0.05 % | 1.8 |
| 20 mT | 5.92 ± 0.04 | 1.19 | 2.6 % | 111 | 2.10 ± 0.08 % | 3.6 |
| 50 mT | 5.91 ± 0.03 | 1.18 | 2.5 % | 111 | 4.44 ± 0.09 % | 7.5 |
| 100 mT | 5.88 ± 0.00 | 1.18 | 2.5 % | 110 | 8.58 ± 0.13 % | 14.6 |
| 150 mT | 5.89 ± 0.00 | 1.18 | 2.7 % | 109 | 12.46 ± 0.06 % | 21.2 |

- The structure does not change over 9 … 150 mT (texture, F90, l_C within
  0.7 %); F90 = 2.5 … 2.7 % stays below the value of one vortex (4.3 %):
  no distinct 180° walls at any amplitude (as 7.18, now n = 8).
- The switched volume grows as B^0.87 ± 0.02 (fit per seed): almost linear,
  the reversible rotation that carries ΔM at μ ≈ 8.
- The texture is almost the same in all seeds (5.87 … 5.92·10⁻³/l_ex,
  ± 0.4 %), except s921 (6.17·10⁻³, +5 %; it loses 0.3·10⁻³ until 150 mT:
  the irreversible rearrangement of 7.17).
- Correlation with β over the seeds: with s921, texture vs β_pow r = +0.95,
  but Spearman ρ = +0.24: one point carries it. Without s921 (n = 7):
  texture vs β r = −0.37 … +0.06 (no correlation). The switched-volume
  exponent correlates with β at 35-70 and 70-150 mT (r = +0.82 / +0.79,
  ρ = +0.86 / +0.89, p ≈ 0.02 … 0.03; 35 tests per table, ≈ 2 false hits
  expected); this is partly expected by construction (loss grows with the
  rotated volume). Hint: texture vs loss level ln w at 100 mT r = +0.82
  (p 0.03, n = 7).
- **Result for the hypothesis "β depends on the total wall area":** at the
  base point the realisations differ too little in texture to test it
  (± 0.4 %), and the texture does not change with B_peak. β follows how
  fast the rotated volume grows with B, not the amount of texture. The
  Sobol design varies L_eff, r_p, φ and d and thus the texture much more;
  the per-run texture from the snapshots is a candidate output there.

### 7.20 Loss decomposition w = c1 b + c2 b² (hypothesis test, 2026-10-04)

Hypothesis (Chris, 2026-10-04): pinned structures give only reversible
rotation (damping loss, linear viscous: w ∝ b², ∝ α); β < 2 needs a part
that grows slower, e.g. irreversible jumps with a fixed energy per jump and
a number ∝ b (w ∝ b, α-independent). Predictions: (i) the fit describes the
data within the cycle errors; (ii) c2 changes strongly with α; (iii) c1 does
not. `eval_decompose.py` → results/decompose.md (all valid runs).

| | α 0.1 (n = 11) | α 0.01 (n = 6) | ratio 0.01 / 0.1 |
|---|---|---|---|
| c1 (jump part) | 3.2 ± 0.5 ·10⁻⁵ | 5.3 ± 1.2 ·10⁻⁵ | 1.64 ± 0.45 |
| c2 (damping part) | 4.81 ± 0.14 ·10⁻³ | 1.92 ± 0.28 ·10⁻³ | 0.40 ± 0.06 |
| median χ²/dof: decomposition / power law | 9.6 / 3.3 | 2.3 / 2.5 | |
| jump share of w at 10 / 50 / 100 mT | 0.44 / 0.16 / 0.09 | 0.76 / 0.44 / 0.30 | |

- (ii) holds partly: the b² part falls to 0.40 ± 0.06 when α falls 10×
  (not ∝ α). (iii) holds within the errors: c1 changes by 1.64 ± 0.45
  (1.4 SE).
- **(i) fails at α 0.1:** the decomposition describes the data worse than a
  single power law with the same number of parameters (χ²/dof 9.6 vs 3.3);
  one run gives c1 < 0 (s921). The reason is the shape: c1 b + c2 b² makes
  β rise with B (from ≈ 1.5 to ≈ 1.9), but at α 0.1 β is highest at small
  B (1.87 at ≈ 18 mT) and lowest at ≈ 50 … 70 mT (7.17). At α 0.01 both
  forms fit equally (2.3 vs 2.5), with a large jump share at small B.
- Result: the hypothesis in this simple form is not supported at α 0.1.
  The α dependence fits it (α acts mainly on the b² part). The high β at
  small B needs a third process that grows faster than b² there (for
  example Rayleigh-like hysteresis, w ∝ b³). A direct test of the jump part
  is the count of jumps in the dissipation time series (samples_*.csv,
  p_dis per sample), not a fit of the form.

### 7.14 Wall measure: method test on synthetic states (2026-10-03, local, no GPU)

Tool `domains.py` (PLAN 3f). Reference: synthetic spheres with known walls
(profile θ = 2 atan(exp(s/L_eff)), width π L_eff) plus a ripple of 20° rms on
cubes of edge L_eff. A real state (d/l_ex 96, dx 3, converged relaxation,
20 mT) has a ripple of 13° rms.

| Measure | Result |
|---|---|
| T(2 L_eff) / wall area, d/l_ex 300 (wall spacing ≈ 5 … 8 L_eff) | 0.89 … 0.94 |
| same, d/l_ex 212 / 150 (spacing ≈ 4 … 6 / 3 … 4 L_eff) | 0.68 … 0.88 / 0.44 … 0.76 |
| same, d/l_ex 300, wall width 0.5 / 2 × π L_eff | 0.92 … 0.95 / 0.75 … 0.88 |
| T(2 L_eff), ripple only, relative to the wall area | 3 … 13 % |
| T(2 L_eff), 4 uniform particles in the FCC cell (no wall) | < 1 l_ex² (smoothing per particle; before the fix 4.2·10⁴) |
| dx 3 vs dx 1.5 (analysed on the same grid) | 0.3 % |
| switched volume / true value (wall moved by 6 and 24 l_ex) | 0.98 … 1.03 |
| switched volume, ripple pattern changed completely | 6 % of V (noise floor) |
| switched volume, uniform rotation by 10° / 30° (no wall motion) | 5.5 % / 16.6 % of V |

- Wall area A_w = T(2 L_eff)/0.85 (± 12 % over the wall widths 0.5 … 2 × π L_eff),
  valid only if the walls are ≥ 5 L_eff apart and not more strongly curved. Else it is the rotation of m on scales above
  2 L_eff. No width solves both: at 1 L_eff the ripple adds 13 … 50 %.
- T counts a 90° wall half. It does not separate walls from a continuous
  rotation (vortex). The localization (∫g)²/(V∫g²) does not separate them
  either (walls 0.61 … 0.89, helix 0.98, vortex 0.48).
- A segmentation into domains (vMF mixture + hidden Markov random field)
  over-counts the walls 2 … 8× (the BIC selects too many direction clusters;
  broad clusters take whole domains). Archived: tag
  archive/domains-segmentation-2026-10-03.
- The switched volume counts every change of m. It does not separate wall
  motion from rotation inside the domains.
- Real state d/l_ex 96 (n = 1): T(2 L_eff) = 2.8·10⁴ l_ex², switched volume
  per quarter cycle at 20 mT 1.0 … 1.9 % of V, l_C = 36 l_ex, mean wall path
  per half cycle 0.6 … 1.6 l_ex.

Literature check and three extensions (2026-10-03):
- Ripple at small fields is long-wave (magnetic SANS, Michels: l_C = L + l_H,
  l_H = √(2A/(μ0 Ms H_i)); for our fields l_H ≈ 1 … 4 µm, estimate with the
  applied field). Test with synthetic ripple of correlation 1 … 8 L_eff
  (13° and 20°): with walls T(2 L_eff)/A = 0.91 … 0.94 (robust); without
  walls the ripple alone gives 0.04 … 0.17 of the area of 3 walls (floor).
- Porod law (S ∝ A/q⁴) as a second wall measure: no q⁻⁴ range for these
  walls (width π L_eff, spacing ≈ 5 L_eff); walls, vortex and helix give the
  same form of the structure factor. Not usable (negative result, as Michels
  reports for nanocrystalline magnets).
- Rotation test F90 (share of windows of 1.5 π L_eff along random lines with
  a net rotation of m above 90°). No peak at the wall angle (oblique cuts
  smear it), but the tail separates 180° walls from ripple and continuous
  rotation if the particle is large: F90 [%] vortex / helix / one 180° wall
  / three 180° walls = 4.3 / 0 / 12.3 / 37 at d/l_ex 300, 9.4 / 6.9 / 18.4
  / 54 at 212; no separation at d/l_ex ≤ 150 (vortex 22 … 65). 90° walls are
  not detected (0.9 %). dx 3 and 1.5 agree. Real state d/l_ex 96: F90 = 66 %,
  the same as a synthetic vortex at this size (65 %): not interpretable.
  Stereology gives a second 180° wall area: a random segment of length W
  hits a wall with probability W S_V/2, thus A_F90 = 2 F90 V/W (d/l_ex 300:
  three walls 1.85 vs 1.86·10⁵ l_ex², one wall −13 %, curved bubble +32 %,
  vortex: false area ≈ 0.3 wall discs). A_w ≫ A_F90 means that part of T is
  continuous rotation.
- Added: correlation length l_C (SANS analogue; synthetic ripple 15 … 78,
  walls 48 … 63, vortex/helix ≈ 100 l_ex) and the mean wall path
  s = V_sw/A_w per half cycle (Bertotti's active magnetic objects; synthetic:
  5.4 for 6 l_ex, −9 % from the calibration).

### 7.7 Pilot (not run; replaced by the planning pilot 7.8 and PLAN.md) (d/l_ex = 150, B_peak = 150 / 100 / 50 / 25 / 9 mT, 7 cycles, ≈ 11 V100-h)

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
| all | scatter of the mean | error of ln w per amplitude ≤ 0.1 (mean over 5 cycles) | more cycles cost budget: ask Chris |

## 8. Error budget

| Quantity | Error (10 % scatter) | How `analyze.py` gets it |
|---|---|---|
| ln w at one amplitude | max(s_cyc, 0.10)/√n (n = 5 kept cycles) | s_cyc pooled over all amplitudes of one run |
| β at 10 / 50 / 100 mT | from a local weighted fit of ln w over ln b (all amplitudes within a factor 2.2); at 10 mT one-sided | `features()` |
| realisation | from the seeds of the base point | variance over the seeds minus the statistical part |
| change of a feature | √(stat.² + realisation² + base²) | `ranking()`; significance with Student-t of the seed scatter |
| dynamic share (`--pair`, 2f) | ≈ ±0.10 | two runs with 10 % scatter each |

Consequences:
- A change of β smaller than ≈ 0.2 is not a result.
- With 3 features and many runs, some |d|/σ above the threshold are false
  hits. Report them as candidates.
- The per-cycle scatter of w_loop can be much larger than 10 % (steady
  test). Then the error of the mean decides, and more cycles cost budget.

## 9. Production

The production study is planned in PLAN.md (frequency test, then a
randomized Sobol design over L_eff, r_p, φ and d/l_ex at dx = 3). The
earlier one-factor and factorial options are replaced by this plan.

## 10. Known limitations

- 4 particles per box, identical spheres, FCC order, no contacts.
- T = 0 (no thermal activation). No eddy currents (add them analytically).
- Herzer length only up to a factor of order 1.
- K_eff is set, not computed: the Herzer averaging from the grains to the
  cubes is assumed, not simulated (section 2.1).
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

Open 2026-10-03: the method for the hysteresis loss (7.9). Options: (A)
quasi-static loops by energy minimization (no damping loss by
construction; validation against LLG at 9 … 20 mT); (B) LLG at a much lower
f (cost ∝ 1/f, ≈ 100× for 100 mT); (C) smaller α and f. Step 2 of PLAN.md
waits for this decision.

Decisions 2026-10-02 for the production (PLAN.md): D1 128 + 8 Sobol points
first, extension to 256 if necessary; D2 7 amplitudes (9 … 150 mT); D3
d/l_ex 212 … 300 with d/L_eff ≥ 7 and a residual check vs d/L_eff; D4
frequency test at the base point and at the corner L_eff 30, r_p 3.


1. After the benchmark: GPU type and budget.
2. After the pilot: the r_p level and the production matrix (section 9).
3. Damping: decided, α = 0.1 (section 7.2).
4. Mesh (decided by the planning pilot, section 4): dx = 3 at d/l_ex = 300.
   Open: the control of one parameter change at dx = 3 and dx = 1.5.
5. Audit points that wait for the mesh test (then Chris decides):
   - separate the hysteresis loss w_h from the damping loss for every point
     (w_h = 2 w(f) − w(2f), +50 % cost), because the damping loss can move β by 0.2 … 0.5;
   - pilot controls at L18 instead of L30 (d/L_eff = 5 at L30 and d/l_ex = 150
     is the no-wall limit), and 2 more seeds at L18 (+0.75 V100-h);
   - pilot pass criteria as fixed tolerances with `--compare` (matched b);
   - ranking by a weighted regression over all runs (the present ranking
     compares only runs that differ in one parameter);
   - production option D: Q_eff (4) × r_p (2) at φ = 0.65 + φ = 0.55 / 0.69 at
     L18 for both r_p + 2 seeds = 14 runs, ≈ 89 V100-h (contains Q × r_p).
