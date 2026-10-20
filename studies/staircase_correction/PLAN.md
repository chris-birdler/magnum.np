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
