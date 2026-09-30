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
| Amplitudes | B_peak = 150 / 100 / 70 / 50 / 35 / 20 / 9 mT (`--b_list`; 9 mT keeps 10 mT inside the range). The drive of each stage comes from b/h of the stage before and is corrected once after cycle 0 with b/h of cycle 0 (factor limited to 0.67 … 1.5) | study range 10 … 150 mT; no time on saturation. b of cycle 0 is within ±4 % of the later cycles (V100 data), thus each stage hits its target to ≈ 4 % |
| Drive | sinusoidal H (not controlled B) | below 150 mT the core is almost linear: B(t) has 0.2 … 3 % harmonics, b varies by 0.1 … 3 % from cycle to cycle (mesh test data). Thus sinusoidal H ≈ sinusoidal B (Steinmetz condition) |
| Initial state | open: decided by the initproto test (section 7.4). Candidates: A = virgin state (random m, full relaxation at H = 0 and α = 1: the T = 0 analogue of the anneal and the cooling without field; a real powder is never magnetized before use), B = AC demagnetization (decaying saturating cycles, as in IEC 60404-6), C = 1 saturating cycle (present) | C gave minor loops around a remanent state (offset up to 0.55 at 23 mT) |
| Protocol | A and B: ascending amplitudes (the classic Rayleigh procedure: every larger loop erases the smaller ones, the loops stay centred). C: descending. 7 cycles each (fixed), cycles 0 and 1 discarded (drive correction) | AC demagnetization. At T = 0 the soft powder does not lock into a periodic cycle (steady test): the mean over cycles is the measurement, not a single steady cycle |
| Loss per cycle | w = mean of w_loop (∮H dB) over the kept cycles; w_dis (LLG dissipation) as a check (`--w dis`) | w_loop of one cycle contains the change of the stored energy when the cycle is not closed; this part averages out over cycles. w_dis scatters less |
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

### 7.7 Pilot (d/l_ex = 150, B_peak = 150 / 100 / 50 / 25 / 9 mT, 7 cycles, ≈ 11 V100-h)

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

## 9. Production matrix (open: Chris decides after the pilot)

Base: d/l_ex = 300, φ = 0.65, L_eff = 18 l_ex, r_p = 0. One d/l_ex = 300 run
(7 amplitudes × 7 cycles) costs ≈ 6.3 V100-h (≈ 0.8 $ on a V100). A 2f partner
run for w_h costs ≈ 3.2 V100-h.

| Option | Runs | Cost | Pros | Cons |
|---|---|---|---|---|
| A: full factorial Q_eff (4) × r_p (2) × φ (2) + 2 seeds | 18 | ≈ 115 V100-h | all main effects and interactions | budget |
| B: Q_eff scan (4) at r_p = 0, φ = 0.65 + 2×2 (r_p, φ = 0.55 / 0.69) at L_eff = 18 + 2 seeds | 10 | ≈ 63 V100-h | Q_eff curve and both other factors with their interaction | interactions with Q_eff not visible |
| C: B + floor (Q_eff = 0) at d/l_ex = 300 | 11 | ≈ 70 V100-h | B + floor at production size | cost |

A dx = 2 control at d/l_ex = 300 costs ≈ 26 V100-h. Thus the mesh check stays
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
4. Mesh in the study range (after 7.5): dx = 3 fails at 9 mT (Δβ ≈ +0.5).
   Options: (i) mesh convergence test dx = 2 / 1.5 / 1 at small B_peak, then
   production at the converged dx (d/l_ex = 300: 26 … 72 V100-h per run);
   (ii) smaller d/l_ex (e.g. 150: 8 V100-h per run at dx = 1.5; this is 1 µm
   only for a material with a larger l_ex, i.e. a smaller Js or a larger A);
   (iii) restrict the study to B_peak where dx = 3 passes (≥ 90 mT, steady test).
   Initial state: A (virgin) at the fine mesh; α = 0.1.
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
