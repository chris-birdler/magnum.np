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
