# PLAN: staircase correction of the stray field — phase 1, static test

Status 2026-10-19: stage A approved by Chris (compute both forms, the full pairwise correction and the tensor).
Main variant: the original method of Donahue & McMichael 2007 (V1, cut cells). V2 (cells with constant
volume) is dropped: it cannot be built locally (PROTOKOLL 3.4). Only the Newell solver with full PBC; only states of
this solver without Bloch points. Basis: PROTOKOLL.md of this folder (3.1 … 3.4) and the β(B_peak) record (7.x).

## 1. Question

Does the near-field correction of Donahue & McMichael let the coarse mesh (Δx 3 l_ex) reproduce the stray field of a
much finer description of the sphere surface? Which form is sufficient: the cell-local tensor (their eqs. 8 and 9)
or the full pairwise correction (their eq. 5)? Phase 1 is static: it compares fields and energies for given
magnetisation states. It does not compute losses.

## 2. Terms and method

- Coarse mesh: the production mesh, N³ cubic cells of edge Δx = 3 l_ex. Voxel cell: a coarse cell whose centre is in
  a particle (the current model). Fine mesh: the same box with q sub-cells per axis and coarse cell (q N)³.
- Fine body: the fine cells whose centre is inside a sphere of diameter d (d from box_edge, as in the run). This is a
  fine staircase of the true spheres.
- V1 cells: each fine cell of the fine body belongs to the coarse cell that contains it. A coarse cell with at least
  one fine cell of the fine body is magnetic. v_i = (number of its fine cells) / q³ is its volume fraction; its moment
  is M_s v_i Δx³ m_i. A magnetic cell that is not a voxel cell ("new cell") takes m of the nearest voxel cell of the
  same particle (centre distance; ties by the cell index). Modified cell: a magnetic cell with v_i < 1.
- Sign convention of magnum.np: H = K M, with the kernel K of demag_f and demag_g (self term about −1/3).
- Couplings, all from demag_f and demag_g with the default Newell radius (20 cell diagonals on each mesh):
  - K^c(i − j): the free-space coupling of the coarse cubes i and j, at the minimum-image offset.
  - L_ij: the mean over the n_i fine cells a of cell i that lie in the fine body, of the sum over the n_j body fine
    cells b of cell j of the free-space fine kernel K^f(a − b). For two cells that are not modified, L_ij equals
    K^c(i − j) inside the window. Because K^f is even and each 3 × 3 block is symmetric, the weighted couplings are
    symmetric:

        v_i L_ij = v_j L_ji
- Window W(i): the coarse cells j with max(|Δi_x|, |Δi_y|, |Δi_z|) ≤ p, minimum image (Donahue's "7 × 7" is p = 3).
- Models, all with the production stray-field setting (DemagField(true_pbc=True, extrapolate=True), 4 shells):
  - voxel model (production): moments M_s m_i on the voxel cells; energy E_vox, field h_vox;
  - volume fraction alone (information; the stray-field part of the model of 7.29): moments M_s v_j m_j, no
    correction; energy E_frac;
  - pairwise (P), energy E_P(p):

        h_i = h_i(volume fraction alone) + Σ_{j ∈ W(i)} (L_ij − v_j K^c(i − j)) m_j

  - cell-local (L), energy E_L(p): the same sum with m_j replaced by m_i, one symmetric 3 × 3 tensor per cell. This
    tensor is v_i times the tensor K_i of Donahue's eq. 9 with the sign of magnum.np (3.3).
  Both corrected operators are symmetric after the weight v_i. P is exact only for a window over the whole box.
- Reference: the fine mesh with M = M_s m_owner on the fine body, the same stray-field setting; energy E_ref. The
  field of a coarse cell is the mean over its n_i body fine cells. The reference is exact for piecewise-uniform M on
  the V1 cells.
- Energies per box volume V_box in units of K_d = μ0 M_s² / 2 (on the fine mesh the same sum over the fine cells):

        E / (K_d V_box) = −(1 / N³) Σ_i v_i m_i · h_i

  Staircase effect: S = E_vox − E_ref. E_K: the sum of the energies of the anisotropy terms (cube and particle) of the
  state per box volume in units of K_d, as run_loops computes them (term_energies); the criteria use |E_K|.
- Field h = H / M_s. Layer: the voxel cells within 3 cells (same window measure) of a modified cell; the same layer
  for all p and all forms. e_h = RMS over the layer of |h_form − h_ref|, divided by the RMS over the layer of
  |h_vox − h_ref|. The field error at the new cells (no h_vox there) is information.
- Double precision in all computations of this test (it compares small differences; no dynamics). The void cells get
  M_s = 0 in all static stray-field computations; the check in step 2 uses the material of the run (void M_s 10⁻⁸ M_s).

## 3. Stage A: production geometry, existing states (local CPU)

**Data.** The three clean prepared states of the AC-300 test runs in studies/beta_bpeak/runs/states_pd300/ (7.114 …
7.118): K2_0_nx_pd300 and K2_1_nx_pd300 (d/l_ex 300, φ 0.62, L_eff 12 l_ex, r_p 3; 152³ cells), K4_1_nx_pd300
(d/l_ex 300, φ 0.55, L_eff 18 l_ex, r_p 1.5). Newell solver with full PBC (newell_ext, 4 shells), single precision.
md5 e4e3aacb…, 22e2c693…, 81e54ade…. The run records give "Bloch points +0 -0" after the last relaxation.
Synthetic fields (uniform m along (100), (110), (111)) only for the check R1 of the computation.

**Settings.** q = 2 (fine mesh about 304³ … 312³ cells; q = 4 does not fit this computer). Window p = 1, 2, 3 for both
forms; p = 4 and 5 if the memory allows (pair matrix about 0.5 … 0.8 GB at p = 3, estimate of the review).

**Steps.**
1. Code (section 5) and the unit tests T0.
2. State checks against the run records: array shape and particle map from the config of the run; |m| = 1 on the voxel
   cells; Bloch points counted again with the detector of run_loops (cube_charge_t, bloch_mask_t), expected 0; total
   energy of the single-precision state, computed in double precision with the production setting and the material
   of the run, against the last relaxation record (E/K_d −0.00749368, −0.00749917, −0.00036765), tolerance
   10⁻⁶ K_d. A state that fails a check is not used; the report says
   why. E_K comes from the same verified setup.
3. V1 cells and the fine body for both geometries (R0).
4. Uniform fields (R1).
5. Smallest eigenvalue of each corrected operator (R3).
6. The three states: S, E_ref, E_vox, E_frac, E_P(p), E_L(p), e_h (R2).
7. Report results/static_A.md and .json; PROTOKOLL entry; context-free review.

**Readings (fixed before the computation).**
- R0, cells: each fine cell of the fine body belongs to a magnetic coarse cell of its own particle; no coarse cell
  holds fine cells of two particles; no voxel cell has v_i = 0 (expected, because the sphere centres lie on cell
  vertices). Reported: number of modified and new cells, distribution of v_i. If a criterion fails: stop.
- R1, computation: for uniform m along (100), (110) and (111), E_ref agrees with φ_f (1 − φ_f) / 3 and E_vox with
  φ_vox (1 − φ_vox) / 3 within 5·10⁻⁴ (relative); φ_f and φ_vox are the volume fractions of the fine body and of the
  voxel model; exact because both bodies have cubic symmetry (7.105). Else: "computation not valid", stop and report.
  E_P(p) − E_ref for uniform m is reported as the error of the cut window (information).
- R3, sign: the operator m → −v ⊙ h(m) of each corrected model is positive semi-definite in the exact case (energy
  not negative). Its smallest eigenvalue (Lanczos, about 100 … 300 field evaluations) must not be below −10⁻⁴ times
  its largest eigenvalue. Else the form is "not usable" at this p, whatever R2 says (negative modes make the
  dynamics unstable). If Lanczos does not converge, the report gives the residual, the reading is "not determined",
  and the form cannot be "sufficient" at this p.
- R2, form, for each form and p, per state:
  - e_E = |E_form − E_ref| / |S| ≤ 0.10, and |E_form − E_ref| ≤ 0.10 |E_K| (the 10 % rule of test A, 3.1);
  - e_h ≤ 0.10 (the same 10 %, chosen for this test).
  - If |S| < 0.01 |E_K|, e_E is not used for this state. If the denominator of e_h is below 0.001, e_h is not used,
    and RMS|h_form − h_ref| ≤ 10⁻⁴ counts instead. The report says which criterion was used. If |S| < 0.1 |E_K| in
    all states, the report also says that the staircase effect of the stray-field energy is already below the 10 %
    rule.
  - A form is "sufficient" at p if all criteria hold in all three states and R3 holds at this p.
  - Baseline: the same criteria are applied to the voxel model (E_vox, h_vox). K2_0 and K2_1 have the same geometry;
    the three states are two geometries.
- Outcome:
  - "go on, cell-local" if L is sufficient at some p ≤ 5;
  - "go on, pairwise" if only P is sufficient at some p ≤ 3 (or 4, 5 if computed);
  - "window not converged" if no form is sufficient, but the largest error of P over the states (e_E and e_h, only
    the criteria in use) at the largest p is more than 0.7 times the value at the p before (a fall of less than 30 %
    per step; a threshold chosen for this test);
  - "not sufficient" else.
  - "Go on" refers to the fine staircase at q = 2, not to the smooth sphere; the target of the correction is still
    open (PROTOKOLL section 2) and stage B gives q = 4.
- Information: E_frac against E_ref (the volume fraction alone); S and S / |E_K|; the errors as functions of p; time
  and memory of the precomputation and of one field evaluation.

**Cost (estimates).** 0 $, no GPU. Code and unit tests about one working day. CPU (8 threads): fine kernel with the
production setting about 1.6 … 1.7 h per geometry (two geometries; measured kernel speed of the review), coarse kernel
about 12 … 13 min per geometry, L_ij minutes to about 1 h, Lanczos about 1 … 2 h, fields minutes per state; total
about 5 … 8 h. Memory up to about 7 GB of the 15 GB. Tokens about 2 … 4 million. The runs need this computer to stay
on.

## 4. Stage B (only after stage A, separate go)

q = 2 and 4 on the B0b geometry (d/l_ex 96, φ 0.65, L_eff 12 l_ex, r_p 0, seeds 941 … 944), fine meshes 96³ and
192³; states prepared on the local CPU with the production protocol (option C of 7.114 with the Bloch-point check, at
most 3 attempts; a seed without a clean state is not used). The same readings. Purpose: the dependence on q and the
link to the fine dynamic references N0b (7.107 … 7.110).

## 5. Code (branch feature/staircase-correction, folder studies/staircase_correction/)

Extend the existing mechanisms: geometry.fcc_box and run_loops (box_edge, locked_geometry, material setup, cube axes,
cube_charge_t, bloch_mask_t), DemagField(true_pbc=True, extrapolate=True), demag_f and demag_g.
- nearfield.py: fine body and V1 cells; L_ij for the modified source cells (the fine kernel convolved with the fine
  occupancy of the source cell, averaged over the fine cells of each target cell; the pairs with an unmodified source
  from the symmetry); the forms P and L; corrected field and energy.
- check_static.py: driver for stage A; state checks; readings; results/static_A.md and .json.
- test_staircase.py, unit tests T0 (code tests on small boxes, not physics checks):
  1. no modified cells: the correction is exactly zero;
  2. symmetry v_i L_ij = v_j L_ji;
  3. open boundaries, a window over the whole small box, and a Newell radius that covers the whole box on both meshes:
     the corrected coarse field equals the fine reference to rounding;
  4. a shift of the box by a whole number of cells does not change the corrected energy (window across the box edge);
  5. a cubic-symmetric fine body with uniform m and full PBC gives φ_f (1 − φ_f) / 3 within 5·10⁻⁴.
- Raw outputs stay local with an md5 manifest; reports and scripts are committed.

## 6. Phases after phase 1

Phase 2: the chosen form as a field term next to DemagField(true_pbc=True) and an option of run_loops (forces the
Newell solver), with tests and a smoke run. Phase 3: a dynamic test of the loss with this term. Each needs its own
PLAN and the go of Chris.

## 7. Limits and risks

- Only prepared states at zero field, not states inside a drive cycle.
- Stage A has only q = 2 (fine staircase of 1.5 l_ex).
- The m of the new cells is a choice (nearest voxel cell); in a dynamic model these cells have their own m.
- The exchange is not corrected (first-order error at the surface, 3.2).
- The Bloch-point detector of the runs does not test the outermost cell layer of a particle.
- Lanczos gives the smallest eigenvalue only approximately; R3 is no proof that the operator is positive.

## 8. Phase 4 (preparation): static reference of the smooth spheres with FE (approved by Chris 2026-10-20, night run)

**Question.** The goal of 3.19 is the stray field of the coarse m (constant in each FD cell, as the FD mesh gives it)
in the smooth spheres. FD gives this only as the limit q → ∞ of X3. An independent reference is FE with the smooth
sphere surface. Before an acceptance criterion is set (decision of Chris), this phase measures what FE can resolve
and at which cost.

**Facts from the first local tests (PROTOKOLL 3.24).**
- magnum.pi 3.13.3 (private clone, docker image magnum.pi.dev-3.13.3) solves the potential on a periodic box with the
  mean field zero, the same convention as the production solver (reading R1).
- The coarse m jumps on every inner cell face. FE resolves these charge sheets badly: the absolute FE energy of a real
  state on the 10-nm mesh is 79 % (P1) or 17 % (P2) too low. Thus only differences of two bodies on the same mesh are
  used (difference method):
      E_s = E_vox(FD, exact) + [E_smooth(FE) − E_vox(FE)].
- The method is checked on a pair with an exact FD value: the 5-nm staircase (X3 at q 2) against the voxel body.

**Night run (local CPU, no cost; fe/run_night.sh).**
1. States: the coarse m of 6 real states (T3Pc_s941 … s944, X3_s951, X3_s952; init.pt). The cells that touch a sphere
   but have no m in the state take m of the nearest cell of the same particle.
2. FD (magnum.np, Newell, full PBC): E_vox, X3 at q 2 and 4, T3Pc (pairwise, p 3, centroid term) at q 2 and 4.
3. FE (P2 potential, M per tetrahedron as its L2 projection onto P1, sampled at sub-tetrahedra; cells at the sphere
   surface with 8⁵ points): mesh 48³ cubes × 6 tetrahedra (conforms to the FD cells) for the 6 states: voxel body,
   5-nm staircase (check ε), smooth spheres.
4. FE with a continuous m (trilinear in the cell centres, unit length) in the smooth spheres: the second question
   (does the whole model, roughness of m included, approach the continuum?). Raw FE value without a correction,
   information only; 2 states, meshes 48³ and 64³.
5. Mesh 64³ (not conforming) for 2 states: does E_s change with the mesh? A sampling check (512 against 4096 points
   per tetrahedron) for the voxel body on this mesh.
5a. Checks added after the review (PROTOKOLL 3.24): ε also on the 2.5-nm staircase (X3 at q 4, exact FD value);
   ε_u on the smooth spheres themselves with uniform m (exact target φ (1 − φ) / 3, a test of the method on the target
   body; uniform m is not a test of the correction); a body without a surface (the whole box magnetic, m of the nearest
   touching cell) to measure the FE error of the inner charge sheets alone.
6. fe/fe_eval.py writes results/fe_static.md: ε, X3(q) − E_s, T3Pc(q) − E_s, continuous m − E_s, all in units of
   S = E_vox − E_X3,q2 of the state, and the cost of each FE run.

**Not in this phase.** The field in each surface cell (the stricter test, as e_h of stage A); a mesh of 96³ (about
30 GB RAM for P2, estimate; vast.ai containers cannot run apptainer or docker, so it needs a VM); an acceptance
criterion.

**Risks.** The difference method assumes that the FE error of the inner charge sheets is the same for both bodies;
ε measures this only for the 5-nm staircase. The sampling of the smooth surface converges slowly (3.24: 64 → 512
points changed E_s by 4 % of S; 4096 → 32768 by 0.35 %). The 64³ mesh may exceed the 12 GB limit of the container.

## 9. Phase 5 (proposal 2026-10-21): dynamic series in q (q 4, q 8) and the way to the smooth spheres

**Question.** The production uses q 2. Statically q 2 covers only about 49 % of the way from the 10-nm roughness to
the smooth spheres, q 4 about 75 %, q 8 about 87 % (3.27, estimate with convergence of first order in 1/q). Does the
loss P change on this way, and does the correction stay exact dynamically at q 4?

**History (PROTOKOLL).** 2026-10-20: the design-A arms D3Pc4 were cancelled before their start (focus on q 2). Then
3.22 decided a dynamic check at q 4 (X3q4 and T3Pc4, seeds 941 … 944, budget 4 … 6 $); it waited for the end of
q 2, which is now done. An early proposal (3.11) already gave "P(q) dynamically at dx 3 and an extrapolation to
q → ∞"; this phase does that. q 16 is not possible with the current code: the coupling table has 343 × q³ × q³ × 6
values, about 280 GB at q 16 (3.25; a convolution form would remove the table: a patch, not built).

**Names** (3.28): Coarse (T3), Coarse m corrected (T3Pc; at q 4 and 8: T3Pc4, T3Pc8), Coarse m exact (X3; at q 4:
X3q4).

**Facts used.**
- Run time of the q-2 runs (stdout.log, median and range of 12 runs; 3.28): Coarse 0.43 h, Coarse m exact 0.85 h
  (0.57 … 1.42), Coarse m corrected 1.29 h (1.12 … 2.17), Fine m 1.31 h, Fine 3.46 h. The GPUs ran 3 … 4.6 runs at
  the same time, thus these are not controlled benchmarks, and the time per step barely grew with the number of cells.
- Time step: 6.6 ps on the 10-nm cells, 1.9 ps on the 5-nm cells (3.28).
- Coarse m corrected at q 4: field of the correction per call 1.7 s against 1.6 s at q 2 (CPU, 3.22). At q 8: set-up
  20 min (8 cores, local) … 43 … 52 min (6 threads, box), peak host memory 10.8 GB (3.25, 3.27).
- Coarse m exact: the Newell kernel of 384³ cells (q 8) took 12426 s (3.27); for 192³ (q 4) about 1/8 of it, about
  0.4 h (estimate). Run time at q 4: 6 … 10 h on a shared GPU (estimate of 3.22).
- Correction against exact field, static: within 1.3 % of S(q 2) for q 2 … 8 (3.27). N3 at q 2, dynamic: scatter
  of the pairs about 0.003 (3.26).

### 9.1 Phase A0: measure the cost on the rented box (first 30 min, before the queue starts)
- For T3Pc4, X3q4 and T3Pc8 and, as reference, T3Pc and X3 at q 2 (production geometry: d/l_ex 96, φ 0.55,
  dx 3 l_ex, N 48; the field terms as in fe/fe_states.py): the time per field evaluation (mean of 20 calls) alone on
  the GPU and with 4 X3q4 at the same time; the peak GPU and host memory; the time step of 200 RKF45 steps from a
  real state.
- Projection per run: run time of the q-2 run of the same arm and seed × ratio of the times per evaluation × ratio of
  the time steps (q 2 / q) + set-up. The number of cycles is adaptive (7 … 15); the projection uses the number of the
  q-2 run (an estimate).
- Rule (fixed now): if the projection of phase A is above its budget (9.7), stop and report before the queue starts.
  The same rule applies to phase B before B starts. A time step below 6.6 ps enters the projection through its ratio.

### 9.2 Phase A: q 4 (dynamic)
- X3q4 on the seeds 941 … 944, as decided in 3.22 (no new choice of seeds after the q-2 data).
- Rule for N3 (fixed now): if fewer than 3 seeds have an event-free pair of T3Pc4 and X3q4, X3q4 runs also for the
  seed 945, then 946 (this order), until 3 pairs exist; at most these 2 more runs.
- T3Pc4 on the seeds 941 … 952 (extension of 3.22).
- Job lines: jobs/95_q4.jobs (checked 2026-10-21: they differ from the q-2 production lines only in --staircase_q),
  plus the same lines for T3Pc4_s945 … s952 and X3q4_s945 … s946. Start state: the prepared state of D15 of each seed
  (--init_transfer; run_loops reads only init.pt, which holds its config; 12 files of 11 MB, local, md5 before the
  upload).
- Code: the current branch; magnumnp/ and run_loops.py are unchanged since 50c36c9, the clean commit of box 3
  (checked with git diff). magnum.pi is not used and not uploaded.
- Evaluation code (local, during the runs; no dependency on the results): check_events.py and eval_dynamic.py get
  the new arms, with the magnetic cells of the q-staircase of each arm (staircase_cells(fine_ids(N, q, …), q)).

**Readings, fixed before the data.** On the seed pairs without an event (rule of 3.28); the event numbers of each arm
are always reported; mean ± SE, 95 % t-interval; band ±0.15 in ln P (K4). Three outcomes for each reading: "within
±0.15" (the interval lies inside the band), "outside ±0.15" (it lies outside), else "not determined"; with fewer than
3 pairs always "not determined".
- N3(q 4) = ln P(T3Pc4) − ln P(X3q4).
- Δ_24 = ln P(T3Pc4) − ln P(T3Pc): the step from q 2 to q 4 at the same coarse m.
- The step from q 2 to the smooth spheres, Δ_2∞. With convergence of first order in 1/q (an assumption: the static
  increments halve, ratio 0.485, 3.27) the steps halve, thus Δ_2∞ ≈ 2 Δ_24 per seed. Reading on 2 Δ_24: "within ±0.15"
  means "q 2 is sufficient for the smooth spheres". Information: if the loss followed the static staircase effect,
  Δ_24 would be about +0.14 and 2 Δ_24 about +0.28 (estimate: Δ_c at q 2 is +0.27 for 1 S, q 2 → 4 adds 0.52 S).
- Δ_c(q 4) = ln P(T3Pc4) − ln P(T3): information.

### 9.3 Phase B: q 8 (dynamic, only the correction)
- T3Pc8 on the seeds 941 … 952; job lines as T3Pc4 with --staircase_q 8, the same start states. N3 at q 8 stays
  static (3.27): X3 at q 8 needs an FFT on 384³ cells in each step and is not practicable dynamically.
- Start rule (fixed now): B starts only if N3(q 4) is "within ±0.15" and the reading of 2 Δ_24 is not "within ±0.15".
  B needs a go of Chris. If he gives it in advance, the rule decides on the box of phase A without a new question.
- The set-up of the coupling table is made again in each run (run_loops has no cache, although the geometry is the
  same for all seeds). Thus at most 2 set-ups at the same time; box with at least 16 cores and 64 GB host memory.
  A cache of the table would be a patch (not proposed now).
- Readings (same outcomes and pairs as 9.2): Δ_48 = ln P(T3Pc8) − ln P(T3Pc4); Δ_2∞ ≈ Δ_24 + 2 Δ_48 (first order),
  reading as in 9.2. Check of first order per seed, Δ_48 − Δ_24 / 2: "consistent with first order" if its interval
  includes 0 (a weak test; information).

### 9.4 Phase C: the smooth spheres (evaluation only, no q 16)
- With B: ln P(∞) ≈ ln P(q 8) + Δ_48 per seed (the rest of a series that halves), with SE; an estimate.
- Without B: ln P(∞) ≈ ln P(q 4) + Δ_24 per seed; an estimate.
- If the check of first order in 9.3 gives "not consistent": no extrapolation; q 16 would need the patch above
  (decision of Chris). The static FE check of the smooth-sphere energy (section 8) stays parked.

### 9.5 Cost (estimates; one box; RTX 3090 on vast.ai, about 0.17 … 0.20 $/h; credit 20.70 $)
The runs of a phase share one box, so the cost is the time of the box.

| Phase | Runs | Time of the box | Cost |
|---|---|---|---|
| A0 | timing | 0.5 h | 0.1 $ |
| A | X3q4 × 4 at the same time (they set the time), T3Pc4 × 12 beside them | 1 … 2 days | 4 … 10 $ |
| A, rule of 9.2 | X3q4 × 0 … 2 more | 0 … 0.5 day | 0 … 2 $ |
| B | T3Pc8 × 12, 3 … 4 at the same time, set-ups 2 at the same time | 8 … 12 h | 1.4 … 2.4 $ |
| C | evaluation, local | 0 | 0 |
| Total | | about 1.5 … 3 days | about 6 … 15 $ |

The upper limit of A assumes that 4 X3q4 at the same time take up to twice the 6 … 10 h of 3.22 (GPU-bound).

**Infrastructure.** One box: RTX 3090 (24 GB), at least 16 cores and 64 GB host memory. Dispatcher run_queue.py
--spool (a job is done only by its DONE file). After the last run: download, md5, report; destroy only with the go
of Chris.

### 9.6 Risks
- The main reading (2 Δ_24) may stay "not determined": the scatter of Δ_24 is not known. N3, which changes only the
  method of the field, scatters by about 0.003; Δ_c, which also changes the surface material, by about 0.34. At most
  8 pairs are event-free (Coarse m corrected had events in 4 of 12 seeds at q 2; events at q 4 reduce this). With the
  scatter of Δ_c, n 8 gives a half-width of about ±0.28 for Δ_24 and ±0.57 for 2 Δ_24. Phase A measures this
  scatter. More pairs need new seeds, and new seeds need new runs of Fine (3.5 h each): the parked option 2 of 3.28.
- The run time of X3q4 is uncertain by a factor of about 2 (A0 and its rule handle it).
- The exchange and the anisotropy keep the coarse cells with the weights v of the q-staircase; the series in q
  changes the stray field and these weights, not the resolution of m (Limits, section 7).
- First order in 1/q is an assumption from 3 static points; the reading of 2 Δ_24 depends on it.

### 9.7 Decisions of Chris
1. Phase A0 and A: q 4 as in 3.22, plus T3Pc4 on the seeds 945 … 952, plus the rule for N3 (at most 2 more X3q4).
   Budget of A: proposal at most 10 $ (3.22 had 4 … 6 $ for the 8 runs without the extension).
2. The readings of 9.2 and 9.3 (main reading on 2 Δ_24 with the assumption of first order) and the rules of 9.1 …
   9.3, fixed before the data.
3. Phase B (q 8, budget at most 3 $): go in advance (the start rule of 9.3 decides on the box) or a separate go after A.
4. Phase C: extrapolation instead of q 16 (no patch).

## 10. Speed patch of the correction (go of Chris 2026-10-21: "Mach den Patch"; PROTOKOLL 3.35, 3.36)

**Aim.** The same correction field, computed faster; no change of the model (same pairs, coefficients, window,
centroid term). Scope: the pair sum and the centroid term. Not in scope (later patches): the set-up memory at q 8
and the host build of the fine geometry.

**Design.**
- Pair sum: a sparse matrix (CSR, int32 indices) with rows 3·target + a and columns 3·source + b, built once on the
  CPU in the set-up, then moved to the device; one sparse product per evaluation instead of gather, 15 temporaries
  and atomic index_add.
- Centroid term under PBC (FFT shape N³): i sin w_k in Fourier space is the periodic central difference
  (f[x+1] − f[x−1]) / 2, and the circular convolution commutes with it. Thus h = Σ_k d_k CD_k(K * v m) −
  K * Σ_k CD_k(v d_k m): 12 FFTs instead of 24. The open (padded) shape keeps the old code.
- Branch perf/staircase-speed from feature/staircase-correction, own worktree; merge into the feature branch only
  after the criteria below pass (not into main).

**Acceptance criteria (fixed before the measurement).**
- C1: the existing tests of test_staircase.py pass.
- C2: new tests, double precision, CPU, periodic and open, random m: max |Δh| / max |h| below 1e-12 for the pair sum
  and for the centroid term (new against old code).
- C3: on the state T3Pc_s941 (init.pt) in single precision: max |Δh| / max |h| below 1e-5, and the term energy
  within 1e-5 relative (the old code itself varies at about 1e-7 from the order of the atomic additions).
- C4 (GPU, V100, geometry d/l_ex 211, q 4): pair sum at most 3 ms (now 20.9 ms), centroid term at most 3.5 ms (now
  6.2 ms), correction at most 6 ms (now 25.4 ms). GPU memory after the set-up is reported (estimate 3.4 GB for the
  matrix at d/l_ex 301, q 4, against 1.9 GB now).
- C5 (GPU): one drive cycle of T3Pc_s941, old against new code: w_loop of the cycles agrees at the round-off level
  of the single-precision runs (reported as a number).
- Cost of C4 and C5: one 1 × V100 box for about 1 h, about 0.2 … 0.3 $ (estimate); renting needs a go of Chris.

## 11. Acceptance of speed patches (decision of Chris 2026-10-22: option 2; criteria to be approved by Chris)

**Aim.** A speed patch must not change the results beyond the round-off of single precision. Round-off cannot be
avoided: each fp32 change of the order of operations changes each field evaluation at the level of 1e-7, and the
trajectories then differ. The yardstick is the spread that round-off alone gives (null runs). Scope: each speed patch
on the production path: the speed patch of PLAN 10, A1 … A6 (time step) and B1 … B5 (correction) of PROTOKOLL 3.36.
A patch with a known change of the model (B1: about 1e-8 relative in the field, through the void moments) is a model
change: Chris approves it as such, and its predicted change replaces the round-off bound of S1.

**Terms.**
- Old: the same model and settings as the new code, without the speed patches. Correction off: the code state of the
  main study. Correction on: feature/staircase-correction at the commit before the bundle. The arms differ only by
  the patches or by the null perturbation.
- Null run, two kinds: (i) repeat: the old code again, where two old runs are not bitwise equal (the old correction
  uses atomic additions); (ii) perturbed start: the old code from a perturbed init.pt. A script of this plan writes it:
  in fp64, m + ε r in the cells with |m| above 0, ε 1e-7, r a fixed random field with |r| 1 per cell (seeds 1, 2, 3),
  renormalized, cast to fp32; the void keeps m = 0. The script reports the RMS and the maximum of |Δm| after the cast
  and the share of changed cells. The old code is not changed.
- State change: cells inside the particles with m · m0 below 0.9 between the start and the end of the cycles (the
  criterion of check_events.py, PROTOKOLL 3.26; the threshold was checked only at d/l_ex 92, 3.28, and is checked
  again on the null runs of stage 4 with the smallest m · m0 and the cells below 0.97 and 0.99).
- Event cycle: a cycle with w_loop at most 0 or |ε| above 0.20 (β(B_peak) record 7.54). Both are reported.
- ln P: the logarithm of the loss density of the measured cycle (P = f K_d w_loop). D: the difference of ln P against
  the old arm on the same state (D_new for the new arm, D_null for a null run).

**Stages 1 … 3: for each patch, before it goes into the bundle.**
- S1, operator level. Each changed function, new against old, on random m and on named states:
  remote_20261021_box2/runs/T3Pc_s941/init.pt (N 48), one phase-1 state of K1 (N 108) and one of K4 (N 156) of the
  main study (open point 2); periodic and open where the code has both.
  - Field terms: max |Δh| / max |h| below 1e-12 in double precision and below 1e-5 in single precision; the term
    energy within 1e-5 of |E_total|.
  - Integrator (A1, A3 … A6): from a fixed state, one RKF45 try: Δm and the error norm below 1e-12 (fp64) and 1e-5
    (fp32) relative; over 1000 steps the same sequence of accepted and rejected tries (fp64).
  - Device placement: one test with the meta device as default device (lesson of 3.36).
  - B4, B5 (set-up at q 8): at N 48 with q 8 (the old code runs there) and at the main geometries with q 2 and q 4.
    At q 8 the main geometries have no old reference (3.35).
- S2, bit equality. Precondition: two old runs are bitwise equal. Then, for a patch that keeps the order of all
  operations (for example A4): m after one sample interval and after one cycle is bitwise equal (torch.equal). If the
  precondition or the test fails, the patch goes to S3.
- S3, one cycle without a state change. One phase-1 state of K1 (N 108), the 9 mT stage; arms: old, 3 null runs (kind
  (i) where it applies, else (ii) with the seeds 1, 2, 3), new. Criteria: |D_new| at most 1e-4 (deciding), and |D_new|
  at most 10 × max(RMS(D_null), 1e-6). Calibration: before the first S3 decision the 3 null runs must have RMS(D_null)
  at most 1e-5; else Chris sets a new bound. If an arm has a state change, S3 runs again on another state, and the
  state change is reported.
- The speed patch of PLAN 10: C5 (test geometry, no null run) is its PLAN 10 criterion; it also goes through S3.

**Stage 4: once for the bundle, before the new code state of the main study.**
- States: 8 prepared states of phase 1 of the main study, 2 for each corner (K<c>_<k>_sat, init.pt; Newell solver,
  static saturation, no Bloch point; β(B_peak) record 7.132 … 7.135). N: K1 108, K2 152, K3 104, K4 156.
- Command line: the production line of the state with --init_from <state>/init.pt --b_list 0.006 --n_amp 1
  --max_cycles 4, without --acc, --acc_harm, --keep_stage_ckpt and --snap_* (as in 7.122); the same line in all arms.
  ln P comes from cycle 3, the last cycle, in all arms.
- Arms: old, one null run (kind (i) where it applies, else (ii) with seed 1), new.
- The analysis uses only the states without a state change and without an event cycle in all arms (rule of Chris in
  3.28). At least 6 such states are necessary; with fewer, the result is "not determined: more states".
- Deciding criterion (equivalence): the 90 % t-interval of the mean of D_new (t(0.95, n − 1), 1.89 at n 8) is inside
  ± δ. Proposal δ 0.01 (1 % of P): the difference N3 of two correction methods on the same start states was
  −0.0009 ± 0.0011 (n 8, without events, 3.28). Option δ 0.03: the cycle scatter of ln P is 2.3 … 2.7 % (7.54).
- Result: pass if the interval is inside ± δ. Fail if |mean D_new| is above δ; then the patches are tested one by one
  with S3 on the states with the largest |D_new|. Not determined if |mean D_new| is at most δ but the interval is not
  inside ± δ: Chris decides between more states and acceptance.
- Information, no fail by itself: the mean of D_new ± SE and the t-test against 0; RMS(D_new) / RMS(D_null); the
  counts of state changes and event cycles in each arm. A count difference cannot be shown at n 8 (3.28: Fisher p 0.22
  for 3 against 0 of 12). A state change in the new arm only goes to Chris as information.

**Cost (estimates).**
- S1: minutes. S3: 5 one-cycle runs at N 108 plus the kernel set-up, about 0.5 … 0.7 V100-h per patch; for 10
  patches about 5 … 7 V100-h.
- Stage 4, correction off: 24 runs of 4 cycles; one 9 mT cycle takes 107 s at N 108 and 253 s at N 156 in production
  (7.122), plus about 2 min per start. The Newell kernels (about 7 min at N 104 … 108 and about 18 … 20 min at N 152 …
  156 on the CPU of box 8) are built once for each geometry before the runs, in an empty cache folder (the cache name
  has no true_pbc and no shell count, 3.36). About 6 … 8 V100-h, about 1 $.
- Stage 4, correction on: the old and null arms run the unpatched correction, which costs about 6 (q 2) … 8 (q 4)
  times a time step without it (3.35): about 35 … 45 V100-h, about 5 … 7 $. At q 4 the old set-up may not fit a 16-GB
  GPU at N 152 and 156 (11.1 GB reserved at N 148, 3.35). At q 8 the old code does not run at the main geometries:
  stage 4 has no old arm there.

**Decisions and open points for Chris.**
1. Approval of this section and the bound δ (0.01 recommended, or 0.03).
2. The init.pt of the phase-1 states (8 for stage 4, 2 for S1 and S3; about 15 … 45 MB each) are on the production box
   of the other session (instance 53995473), which this session does not touch. Download by that session, or by this
   session with a go of Chris.
3. Which patches go into the bundle; each patch needs its own go ("patch").
4. Timing: a bundle between phase 1 and phase 2 of PLAN 4.15 breaks its condition 3 (one commit for both phases); the
   phase-1 runs would then be a test only. A merge into main needs a go of Chris.
5. Stage 4 with the correction off, or on with which q (cost above).
6. The evaluation scripts of S3 and stage 4 (perturbed start, state-change check with free run names) belong to this
   plan; they do not change the simulation code.

**Decisions of Chris (2026-10-22, about 06:00 UTC).** Section 11 approved, with δ 0.01 (the recommendation). Download
of the start states approved ("nur Download"): 10 init.pt of phase 1 (K1_0, K1_1, K1_2, K2_0, K2_1, K3_0, K3_1, K4_0,
K4_1, K4_2 _sat) were copied from instance 53995473 with scp, md5 checked (remote_20261022_prod/, gitignored,
manifest prod_states.md5); nothing else was done on that instance. Stage 4 uses K1_0, K1_1, K2_0, K2_1, K3_0, K3_1,
K4_0, K4_1; S1 and S3 use K1_2 (N 108) and K4_2 (N 156).

## 12. Build of the speed patches (go of Chris 2026-10-22: "alles freigegeben", as proposed in the chat)

- Patches: all items of the audit (PROTOKOLL 3.36) except A5 (CUDA graphs; Chris) and A6 (one batched FFT; tested in
  the fp32 campaign of 2026-07 and rejected: slower than the loop over the components at both precisions on the V100;
  bench/fp32_validation.tex). B1 is approved as a model change (about 1e-8 relative in the field).
- Branches: A1 … A4 on perf/timestep from main (own worktree ~/Documents/code/magnum.np-timestep; the worktree of main
  is used by another session). B1 … B5 on perf/staircase-speed after the patch of PLAN 10. No merge into main or into
  the feature branch without a go of Chris.
- Order: A4, then A2, A1, A3 (order of these three after the profile of the time step), then B5, B4, B2, B3, B1 (the
  model change last).
- For each patch: one commit; S1 (local fp64 on the CPU, fp32 on box 8), S2 or S3 on box 8 (section 11); the result
  in the PROTOKOLL; push of the branch.
- Stage 4 at the end, once, with the correction off (tests A1 … A4). Stage 4 of the B patches only after the decision
  of Chris on the correction in the main study and on q. The use of the new code in the main study (phase 1 and 2 of
  PLAN 4.15) is decided after stage 4.
- Box 8 (0.135 $/h) stays for the tests between the patches (go of Chris).
- Cost (estimate): S3 about 0.5 … 0.7 V100-h for each of 8 patches, stage 4 about 6 … 8 V100-h; with the waiting time
  of box 8 during the build about 3 … 7 $ in total.

**Decision of Chris (2026-10-22, about 08:00 UTC): B2 and B3 are not built.** Measured on box 8: B2 (3 × 3 block
matrix, BSR) has 45 % less matrix memory but is not faster (0 … 16 % slower on the V100); B3 (FCC translation
symmetry) saves only memory and depends on the geometry. With B4 and B5, q 8 fits a 16-GB GPU at both main geometries.
Both are a later option for q 16 or larger boxes.

## 13. Precision test of the speed patches (decision of Chris 2026-10-22, about 10:30 UTC; replaces stage 4)

**Reason.** S3 in single precision cannot decide (PROTOKOLL 3.39). The start-to-end criterion of a state change marks
about 130 000 cells in every arm, because the first cycle of a prepared state changes m in many cells. RMS(D_null) is
4.5e-5 … 1.2e-4, above the calibration bound 1e-5: the adaptive step control takes a different number of tries in each
arm (5311 … 5373 accepted steps in the cycle). Stage 4 decides on a statistic over few states, with events excluded.
Chris wants a verification that the scatter of the dynamics cannot hide.

**Principle.** A patch that only changes the order of operations changes each field evaluation by round-off. The unit
round-off is 2^-24 (6.0e-8) in single precision and 2^-53 (1.1e-16) in double precision. A real change of the computed
function (a fault or a model change) keeps its size in both precisions. In double precision, two arms with round-off
differences only take the same sequence of accepted and rejected tries with a high probability, and their trajectories
stay together. The test measures this sequence; it does not assume it.

**Set-up.**
- Each patch: its S3 run pair (old, new) again, with --precision_cycles double. This is the existing mechanism of
  2026-10-17: the run loads the init.pt and does all of its own work in double precision; no arm relaxes; the Newell
  kernel is computed again (run_loops.py gives no cache folder); the coefficients of the correction take the dtype of
  the state. The other arguments (also atol_cycles 1e-7, rtol_cycles 1e-6), the start state and the code trees are
  those of S3; only the precision changes.
- Correction off: K1_2_sat (N 108). Correction on (q 4): the start state of the B S3 (K1_2_sat transferred to the
  fraction surface, bs3/prep/init.pt).
- step_dump.py --trace saves m after 1, 2, 4, …, 256 sample intervals (one drive cycle), and after each sample
  interval the counts of RKF45 tries and of accepted steps, t, and m in 2048 probe cells (a fixed choice of magnetic
  cells, the same in all arms).
- Positive control, one for each code branch: the old code (9b15ea9 for A, 0125db3 for B) with the drive amplitude
  changed by the relative amount 1e-6 (--b_list 0.006000006 instead of 0.006). It shows that the test sees a known
  small change, and it gives the detection limit: the relative drive change that gives ρ = 1e-3 ρ_null (linear
  response).
- One geometry (Chris): the code paths do not depend on the geometry, and S1 tests the operators on several sizes and
  boundary set-ups.

**Pairs (code trees of S3).**

| patch | old | new |
|---|---|---|
| A2 | 9b15ea9 | 31aa455 |
| A1 | 31aa455 | c0da045 + the llg.py change of 3d9de71 (tree A_a1f on box 8, md5 of git diff db467c1c…) |
| A3 | c0da045 + the llg.py change of 3d9de71 | 3d9de71 |
| PLAN 10 | c63e862 | 01ed90c |
| B5 | 01ed90c | c56c0a5 |
| B4 | c56c0a5 | d538bde |
| B1 | 0125db3 | df0f811 |

A4 passed S2 (bitwise equal). B5 and B4 change only the set-up; if their S3 end states in single precision are
bitwise equal (m and t after 256 sample intervals), they pass as in S2 and do not run again. 0125db3 changes only
energy_budget.py against d538bde (not on the run path). Runs: A 4 + 1 control; B 5 (7 if B5 or B4 is not bitwise
equal) including the control.

**Criteria (fixed before the data).**
- ρ: the RMS over the magnetic cells of |m_new − m_old| after the cycle. ρ_null: the median ρ of the 3 null runs of
  S3 in single precision of the same code branch (old code, start m + 1e-7 r), that is the scatter of production runs
  that round-off alone causes. The yardstick is the same for all patches of a branch.
- Patches without a model change (A1, A2, A3, PLAN 10, B4, B5): pass if ρ in double precision is at most 1e-3 ρ_null;
  fail otherwise. The bound has a wide margin: when the tries agree, a round-off difference in double precision is
  of the order of the unit round-off times the gain that the positive control shows (estimate).
- Not determined: the counts of tries or accepted steps of the two arms differ after a sample interval (one changed
  accept or reject decision). The result then gives the probe ρ up to the interval before; Chris decides.
- B1 (approved model change; ms_void, the ratio Ms_void / Ms without a unit, is 1e-8 in production; S1 measured a field
  change of 3.2e-9 relative with it, and 1.5e-16 in double precision without void moments): S1 has verified the code.
  The double-precision run measures only the size of the model change. Criterion (approved by Chris): pass if ρ in
  double precision is at most 0.1 ρ_null, that is the model change is at least 10 times smaller than the scatter that
  round-off causes in production.
- Information, no criterion: ρ and the largest angle in degrees between m_new and m_old after each saved sample
  interval; the probe ρ after each sample interval; the counts of tries and accepted steps in both arms; D_new in
  double precision; the ρ of the pair in single precision (S3).

**Limits of the test.**
- torch.compile makes its kernels for each dtype: the double-precision runs do not run the single-precision kernels of
  production. Those are checked only by S1 in single precision on the GPU and by S3 in single precision (the A1 fault
  of 3.39 was found there).
- Open finding: in single precision, ρ of the A2 pair is 1.0e-5, 5 times ρ_null (1.9e-6); A1 and A3 give ρ_null. The
  cause is not known (assumption: a systematic difference of the rounding of the fp32 kernel, 1/n and the block term
  are now in the kernel). The double-precision test shows if A2 changes the computed function.
- Only q 4 with the correction. Production q 8 (if Chris decides for it) uses the same code paths with more pairs; B
  at q 8 is open until the q decision (section 12).

**Cost (estimate).** A: 5 runs at N 108 without the correction, up to 3 at a time. B: 5 … 7 runs with the correction
q 4, one at a time (about 2 × the memory of single precision). About 3 … 4 V100-h after the B S3, about 0.6 $ for
box 8.

**Changes to sections 11 and 12.** Section 13 replaces stage 4 (section 11) and the stage-4 steps of section 12; the
use of the new code in the main study is decided after section 13. The S3 bounds of section 11 (|D_new| at most 1e-4,
calibration 1e-5) and the S3 proposal of 3.39 are not used; S3 in single precision gives ρ_null and stays as
information. The field change of B1 is 3.2e-9 relative (measured, S1); "about 1e-8" in sections 11 and 12 was the
estimate. The scripts s4_a.sh and s4_eval.py stay as a record; they did not run.

**Decision of Chris (2026-10-22, about 11:00 UTC).** Section 13 approved as written, with the B1 criterion 0.1 ρ_null
("Ok").

## 14. Verification over the full cycle and final speed measurement (decision of Chris 2026-10-23: option 2; criteria to be approved by Chris)

**Reason.** Section 13 gave "not determined" for all 7 pairs (PROTOKOLL 3.39). In double precision, the RKF45 step
control also makes a different decision after 8 … 43 of 256 sample intervals. After that, ρ grows to a level that does
not scale with the unit round-off. Before it, the data separate clearly, but only for the sample intervals 1 … 4,
where the drive is below 10 % of its peak. Also, no test compares the single-precision kernels of production over the
cycle, and the A2 pair has 5 times ρ_null in single precision. Chris wants all patches tested in a conservative way,
and then the final speed gain.

**Principle.** A difference that grows from round-off needs time to change a step decision. Thus, if all arms start a
short window from the same state with the same time step, they make the same decisions for some steps. In this window,
the difference of the arms shows the difference of the computed functions without the scatter of the dynamics. A
positive control in the same window gives the sensitivity of the test at that phase of the cycle. A null run, which
changes only the order of the sums, shows the level of a pure round-off difference.

**Terms.**
- Arm: one run of one code tree. Pair: the old and the new arm of one patch.
- Control: the old code of the comparison with the drive amplitude changed by the relative amount 1e-6 (as in section
  13). The drive amplitude is H_amp = (0.006 / 8) Ms = 7.5e-4 Ms (run_loops.py, h_amp = b_target / mu_guess), so the
  control changes the applied field by 7.5e-10 Ms (one control unit).
- Null run: the old code of the comparison with the field terms summed in the reverse order. It changes each field
  evaluation by round-off, as a patch that changes the order of operations does.
- m0: m of the start state (init.pt) of the geometry. Particle cells: the cells with Ms above Ms_void, including the cut
  cells of the surface model. Void cells: Ms_void = 1e-8 Ms, exchange constant 0, a free m.
- Difference measures of two runs at the same point: ρ_all, the RMS over the cells with |m0| above 0 of |Δm| (the ρ of
  section 13); ρ_part, the same over the particle cells; Δ_max, the largest |Δm| over the particle cells.
- Checkpoint k0: the end of sample interval k0 of the reference run, k0 = 0, 16, 32, …, 240 (16 checkpoints; 64 is the
  drive peak, 128 the zero crossing, 192 the negative peak). Window: the sample intervals k0 + 1 … k0 + 5.
- RKF45 stage: one evaluation of the right-hand side; a try has 6 stages.
- Try record, for each RKF45 try: dt, the error norm as the tree itself computes it (old trees in _optimal_stepsize, A3
  trees in _error_norm), the decision (accept or reject), and the branch: a tuple of the cases of rkf45.py (A4 clamp of
  dt to the rest of the interval before the try; norm above 1.1, from 0.5 to 1.1, below 0.5; clamp at minscale or
  maxscale; clamp at maxstep; the retry with min(dt_opt, remaining)). The wrapper wraps only _try_step and the norm
  function of the tree (as step_dump.py --trace does now); it does not replace step().

**Code trees.**
- A geometry (K1_2_sat, N 108, correction off): 9b15ea9 (A4; old of A2), 31aa455 (A2), A_a1f (A1, as in section 13),
  3d9de71 (A3), AB.
- B geometry (bs3/prep/init.pt, q 4): c63e862, 01ed90c (PLAN 10), c56c0a5 (B5), d538bde (B4), df0f811 (B1), BA4, AB.
  0125db3 differs from d538bde only in energy_budget.py, which is not on the run path; thus the arm of d538bde is also
  the old arm of B1.
- AB: a new branch perf/speed-all from perf/timestep (3d9de71) with a merge of perf/staircase-speed (df0f811). There are
  two conflicts (git merge-tree): demag.py (B1 stores the last stray field in DemagField.h for the correction; in AB
  this goes into the A2 version) and run_loops.py (NEW_KEYS and INIT_KEYS: sat_static and h_sat of main against the
  keys of the correction). Checks: AB writes the same configuration (config.json) as the old tree for the A line
  (--h_sat 2 --sat_static) and for the B line; all tests of both branches pass locally with torch 2.7.1; one short run
  on the production path on the GPU at N 48 with the correction off and on (lesson of 3.39).
- BA4: df0f811 with 9b15ea9 (A4) cherry-picked. The B trees do not have A4: in each sample interval they make one try
  past the end, which is always rejected; AB and BA4 do not make it. Thus their tries line up only with each other.
- No merge into main or into the feature branch.

**A4 is not exact in one case.** In the old code, the try past the end gives dt_opt from its own error norm, and the
retry uses min(dt_opt, remaining). If that dt_opt is below the rest of the interval (a norm above 1.1 in the try past
the end), the old code takes a shorter step than A4. Then the step sequence differs at the level of the truncation
error, not of round-off. In S2 (A geometry, single precision, one cycle) A4 was bitwise equal, so this case did not
occur there. K1 counts the case on the B geometry (below).

**K1: windows over the cycle (double precision).**
- Reference run, one for each geometry: the old code (A: 9b15ea9; B: d538bde) with the arguments of section 13
  (--precision_cycles double, one cycle). It integrates the full cycle and saves m, t and the RKF45 time step at each
  checkpoint. The B reference also counts, over the full cycle, the tries past the end with dt_opt below the rest of
  the interval (the A4 case).
- Arm run (step_dump.py --resync): before interval k0 + 1 of each checkpoint, the wrapper sets m (with copy_ into
  state.m, so that the cache of the stray field of B1, keyed by the tensor and its version, becomes invalid), t (a
  float64 tensor on the device) and the RKF45 time step (_dt) to the reference values. The solver keeps no other state
  between two calls (rkf45.py: only _dt), and the drive is a function of state.t only (run_loops.py). Thus all arms
  start each window from the same state. Outside the windows the wrapper does not integrate; it only moves t to the
  end of the interval (this saves about 2/3 of the time; sample() and the diagnostics of run_loops.py do not change
  the state, and with --max_cycles 1 no drive correction occurs).
- Window record: the try record; m in the 2048 probe cells of section 13 after each accepted step; the full m after
  k0 + 1, k0 + 2 and k0 + 4.
- Comparisons, each with its own old arm, control and null run:
  - A geometry, old 9b15ea9: A2 (→ 31aa455), A1 (31aa455 → A_a1f), A3 (A_a1f → 3d9de71), AB off (3d9de71 → AB).
  - B geometry, old d538bde: PLAN 10 (c63e862 → 01ed90c), B5 (01ed90c → c56c0a5), B4 (c56c0a5 → d538bde), B1
    (d538bde → df0f811).
  - B geometry, old BA4: AB on (BA4 → AB). AB on is the only test of A1 … A3 with the correction on; if it fails, the A
    patches are tested one by one on the B geometry (separate go).
  - A4 with the correction on (df0f811 → BA4): the tries are aligned by accepted step (the try past the end of
    df0f811 is left out). If the B reference has no A4 case, criteria 1 and 2 apply; else the number of cases goes to
    Chris (A4 then changes the step sequence within the tolerance of the integrator), and Chris decides.
- Arms: A: reference (= old arm 9b15ea9), 31aa455, A_a1f, 3d9de71, AB, control, null. B: reference, d538bde, c63e862,
  01ed90c, c56c0a5, df0f811, control and null of d538bde, BA4, AB, control and null of BA4. B has a separate old arm,
  because the B code is not reproducible in single precision (3.39).
- Determined try: the old arm, the new arm and the control have the same decision and the same branch in all tries
  from the checkpoint up to this try. Determined point: the end of interval k0 + 1, k0 + 2 or k0 + 4, if all tries up
  to it are determined. The small continuous differences of dt are part of the measured difference.
- Criterion 1 (m): at each determined point, ρ_all, ρ_part and Δ_max of the pair are each at most 1e-3 × the same
  measure of the control. The bound is the response to a relative drive change of 1e-9 (linear in the window). In
  section 13 at k0 0, the largest ratio of a round-off pair was 1.5e-4 (A2).
- Criterion 2 (step control): at each determined point, the RMS over the tries up to it of the relative difference of
  the error norms of the pair is at most 1e-3 × the same value of the control. This tests A3 and A4, which change the
  computation of the norm and of dt, and not only m.
- Pass: each of the 16 windows has at least one determined point, and criteria 1 and 2 hold at all determined points.
  Fail: one value above a bound. Not determined: a window without a determined point. Thus a fault that changes the
  decisions at once gives "not determined" in many windows, never "pass". Each window is reported; none is dropped.
  The number of windows that are not determined is compared with that of the null run.
- Pilot, before the other arms: the reference runs, the controls and the null runs. Stop rule: if a null run leaves a
  window without a determined point, or its ratio to the control is above 1e-4 in criterion 1 or 2 at a determined
  point (a margin of 10 below the bound), I stop and report to Chris before the other arms run.
- B1 (model change): K1 gives its ratios to the control in each window (information for K3).

**K2: single-precision kernels on the checkpoint states.**
- For each tree, on the production path in single precision, on the 16 checkpoint states cast to single precision
  (m32): the field of each term and the effective field, evaluated in the order of production (for the trees with B1,
  DemagField.h before the correction, and a check that the correction used the stored stray field), and the LLG
  right-hand side; one RKF45 try (_try_step with _dt set to the reference value, not step(), which A4 changes): the six
  stage values, the increment x1 − x (formed in double precision), the error estimate (RMS) and the error norm.
- Reference: the same quantities in double precision from m32.double() (the same input as in single precision, so
  that the rounding of the input is not in the error). Reference trees: 9b15ea9 (A trees without B1), d538bde (B trees
  without B1), df0f811 (the trees with B1: df0f811, BA4, AB).
- For each quantity q and state: the error e = RMS of |q_single − q_double| / RMS of |q_double|, and the bias
  b = Σ (q_single − q_double) · q_double / Σ q_double · q_double. A coherent scale error (for example in the kernel copy
  of A2, which holds 1/n in single precision) acts like a change of the field; the bias shows it. A quantity with a
  reference of 0 (the drive at k0 0) gets its absolute difference and no ratio.
- Each single-precision evaluation runs 3 times (the pair sum of c63e862 uses atomic additions); e and b use the median.
- Criterion, for each quantity: the RMS over the 16 states of e of the new tree is at most 2 × that of its old tree, and
  at most 2 × that of the base tree of its branch (A: 9b15ea9; B: c63e862), so that the factors do not multiply along
  the chain. The same for |b|, with a floor of 2^-24 (6.0e-8, the unit round-off of single precision) for the old
  value.
- A2: if its ratio is from 1 to 2, A2 passes, and the cause of its 5 × ρ_null in S3 stays open in the report.
- The kernel cache: one folder for each geometry and precision; the script checks the dtype of a kernel when it loads
  it (the cache name has no dtype, and demag.py only casts a single-precision file up).

**K3: B1 (approved model change).**
- B1 meets the bound of section 13 already: ρ_all in double precision is 3.7e-8, 0.013 ρ_null (bound 0.1), with the step
  differences included (they add scatter). Only the rule on the step sequence made it "not determined".
- Explanation: from the double-precision values of df0f811 and d538bde on the 16 states (K2, in the order of production,
  with the check of the stored stray field): Δh in the cut cells (B1 changes the field only there; in full and void
  cells it is 0 by design), RMS and maximum, in control units.
- Hypothesis, fixed before the data: the RMS of Δh in the cut cells is 0.1 … 10 control units (estimate: about 1 … 2,
  at most about 1.5e-9 Ms). The 3.2e-9 of S1 was relative to the maximum of the correction term in a unit test (N 16,
  random m, q 2), not to the field of a production state.
- Result: B1 passes K3 if the bound of section 13 holds (it does) and the RMS of Δh is inside the band. Outside the band,
  the hypothesis is false; I report, and Chris decides.

**Speed measurement (after K1 … K3 pass).**
- Tool: prof_step.py on the production path, single precision, one run on the GPU at a time, old and new runs in turn.
  Extension: prof_step.py measures the set-up time (start to the first sample interval) and the peak GPU memory, and
  gives DemagField a cache folder for the Newell kernel (run_loops.py gives none). One folder for each geometry; if old
  and new write different cache files (md5), each code gets its own folder. The kernel build is measured once for
  each geometry and reported separately.
- Old: main 9e0ea7c with the correction off (the code state of the main study; physics code as 8640e55); c63e862 (the
  feature branch before all speed patches) with the correction on. New: AB.
- Configurations: K1_2_sat (N 108) and K4_2_sat (N 156); correction off; correction on with q 4; q 8 with the new code
  only (no gain factor: the old code does not fit the 16-GB GPU there, 3.35). If an arm does not fit the GPU (old q 4
  at N 156, or q 8 at N 156), that is the result for it, and the gain factor of that configuration is not given.
- Start states with the correction (staircase_q is an INIT_KEY): bs3/prep for q 4 at N 108, and 3 transfer runs as for
  bs3/prep: q 4 at N 156, q 8 at N 108 and at N 156.
- Each run: 30 sample intervals of warm-up, then 96 measured sample intervals (the measured part includes the drive
  peak). 3 runs for each arm and configuration (n 3, standard error). Check: the ms per RKF45 stage of the two halves of
  the measured part differ by at most 10 %; else the run repeats with 100 intervals of warm-up.
- Result: ms per sample interval (the time of production; it includes the try that A4 saves), ms per RKF45 stage,
  tries per sample interval, set-up time, peak GPU memory; the gain factor old / new ± SE for each configuration with
  an old arm; the time for one drive cycle at each N and configuration. The projection to the main study uses its
  number of cycles for each N and is marked as an estimate.

**Limits.** K1 and K2 use one state of each geometry (decision of Chris for section 13) and only q 4; q 8 uses the same
code paths with more pairs (S1 tested B4 at N 48 with q 8: 6e-16). The double-precision kernels of K1 are not those of
production; K2 covers the single-precision kernels on 16 states.

**Order and stop rules.**
1. Approval of this section and of the decisions below.
2. Local work: step_dump.py --resync with the try record; s1_states.py (an extension of s1_demag.py); k_eval.py; the
   extension of prof_step.py; the trees AB and BA4. Local tests on the CPU (torch 2.7.1, T3_s941, N 48): the resync
   restores m, t and dt bitwise; two equal arms give 0 in all measures and equal try records; the null run keeps
   windows determined and gives a small ratio to the control; an injected change of a field shows in the window where
   it starts; an injected change of the error norm shows in criterion 2; the try record of an A4 tree and an old tree
   align by accepted step.
3. GPU: start box 8 again (needs credit; a new box, if its GPU is gone, needs a go of Chris). Check the free disk (about
   60 GB are necessary). Short checks on the production path (N 48; one window at N 108).
4. Kernel builds (one folder for each geometry and precision) and the 3 transfer runs; they run beside steps 5 and 6.
5. K1 pilot; stop rule above; report.
6. K1 arms, K2, K3. Report. A fail or a window that is not determined stops the work for that patch until Chris
   decides.
7. Speed measurement only after K1 … K3 pass. Report.
- Disk on the box: about 1.5 GB per K1 arm (full m at 3 points in 16 windows; about 22 arms, about 33 GB) and about 7 GB
  per K2 reference tree (3 trees, about 21 GB). The raw data stay on the box; the aggregate results go to git.

**Cost and time (estimates).**
- Local work: about 8 … 10 h (scripts, AB and BA4, tests); no GPU.
- Box time: checks about 1 h; kernels and 3 transfer runs about 1 h (beside other runs); K1 pilot about 1 h; K1 arms
  (windows only, about 1/3 of a cycle each; A 3 at a time, B one at a time) about 1.5 … 2 h; K2 (about 15 set-ups,
  kernels from the cache) about 1.5 … 2 h; speed measurement (30 runs) about 4 h. Sum about 10 … 12 h, with waits
  about 14 h: about 1.9 $ at 0.135 $/h; with a reserve for repeats or a new box, 2 … 3 $.
- Tokens: about 3 … 5 million.
- Wall clock: about 2 days, with the long GPU phases overnight.

**Decisions of Chris for this section.**
1. K1: the bound 1e-3 × the control for criteria 1 and 2 (a relative drive change of 1e-9; a new yardstick, not 1e-3
   ρ_null).
2. K2: the factor 2, and the floor 2^-24 of the bias.
3. B1: the bound of section 13 (met) plus the band 0.1 … 10 control units.
4. A4: the rule for the case of a shorter step in the old code.
5. The stop rule of the pilot (margin 10).
6. A new box, if the GPU of box 8 is gone.

**Changes to section 13.** Section 14 replaces the criteria of section 13 for the patches without a model change. For
B1 the bound of section 13 stays, and K3 replaces its rule on the step sequence. The data of section 13 stay as
information.
