# Protokoll: staircase correction of the stray field (magnum.np, Newell solver with full PBC)

Owner: Christoph Vogler (Chris). Branch: `feature/staircase-correction` (from `main` 6ca3b65, 2026-10-19). Folder:
`studies/staircase_correction/`. Literature notes: `studies/staircase_correction/lit/`.

Rules: documents in ASD-STE100 Simplified Technical English, SI units only.
**Chris makes the decisions.** Before you change the model, the method, the test design or the budget, give him the
options with pros, cons, chances and risks, and wait. Engineering fixes (bugs, logging, speed) need no decision, but
tell him.

References: a number "7.x" is a section of the record of the β(B_peak) study, `studies/beta_bpeak/PROTOKOLL.md`
(branch `study/beta-bpeak-pinning`). A number "3.x" is a section of this record.

## 0. Current state (2026-10-20)

- **Terms.** l_ex: exchange length. d: particle diameter. Δx: cell edge. P: loss density in W/cm³; β: exponent of
  the loss curve over B_peak. SE: standard error; n: number of independent seeds. Staircase: the surface of a sphere
  built from cubic cells. Near field: the stray-field coupling between a cell and its neighbours up to a few cells
  away. FFT: fast Fourier transform. PBC: periodic boundary conditions. Newell solver with full PBC:
  DemagField(true_pbc=True), the FFT convolution of M with the Newell tensor of cubic cells, with image shells and the
  compensation of the outer surface charges (7.105). DemagFieldPBC: the older finite-difference stray-field solver,
  with errors at edges and corners and cell-scale null modes (7.105). V1: cut boundary cells with volume-fraction
  moments (Donahue & McMichael 2007). V2: boundary cells with constant volume and moment, deformed to fill the smooth
  sphere (idea of Chris; definition in 3.2).
- **Question.** Can a near-field correction of the stray field remove the dependence of the loss on the voxel
  staircase of the coarse mesh (Δx 3 l_ex)? Hypothesis of Chris (2026-10-19): the stray-field near field of the
  staircase causes this dependence. Basis: on a locked voxel geometry the loss converges (7.110, V3g). On the mesh of 0.75 l_ex the
  loss at 9 mT with the voxel geometry of this mesh is about 3 times the loss with the voxel geometry of dx 3:
  difference of ln P −1.15 ± 0.41 (SE, n 4; 7.110, arm F, N1g − N0b, d/l_ex 96, 9 mT). This attribution to the geometry is an inference of 7.110,
  because the stray-field kernel and the cycle tolerance also differ between the two arms.
- **Literature (3.1 … 3.3).** No source tests a staircase correction on a body that is curved in all three
  directions, or on hysteresis or loss. Donahue & McMichael 2007 give the method (far field by FFT, near-field
  correction from a submesh) and use its cell-local form for all their results (3.3). The variant V2 (boundary cells
  with constant volume, idea of Chris) was not found (3.2) and cannot be built locally (3.4).
- **Method now (3.7).** V1 (cut cells after Donahue & McMichael) plus a first-moment term in the far field: the
  moment of each cut cell sits at the centroid of its body part (central difference of the same periodic kernel; no
  new kernel). Without this term the window converges only as 1/p (3.6).
- **Static result (3.6 … 3.9).** With the centroid term both forms (pairwise P, cell-local L) meet the criteria R2 of
  PLAN.md from p = 2 on, on the production geometry (d/l_ex 300, 3 states, q = 2 and 4) and on the geometry of the
  dynamic test (d/l_ex 96, 4 states, q = 2 and 4). At p = 3: energy error 1 … 2 % of the staircase effect S; field
  error in the surface layer 0.2 … 0.4 % (P) and about 6 % (L) of that of the voxel model. The staircase effect S of
  the stray-field energy: 0.05 … 0.08 |E_K| (K2), 0.23 … 0.39 |E_K| (K4), 1.0 … 2.0 |E_K| (d/l_ex 96), larger against
  the finer target (q = 4) than against q = 2 (factor 1.7 … 1.8).
- **Dynamic test (3.14, 3.15).** Design A (own preparation per arm) has no power (SE of d 0.65 against the bound
  0.15). Design B (transferred start, n 4): d = +0.07 ± 0.13 (SE), reading "not determined" (2 SE above 0.15); 3 of 4
  seeds agree within 0.18 in ln P, seed 944 deviates by +0.42 without a change of state.
- **Open.** Next step of the dynamic validation (decision of Chris); the use in the main study
  (decision of Chris). Target for now: the staircase of dx 1.5 (q = 2; decision of Chris 2026-10-20).
- **Status.** Night work 2026-10-19/20 under the night mandate (decisions D1 … D8 in 3.5 … 3.10). Since 2026-10-20
  morning all decisions go through Chris again (3.11).

## 1. Goal

A coarse mesh (Δx 3 l_ex) whose stray field represents the smooth sphere, or at least a much finer staircase, so that
the loss and β of the powder model do not depend on the coarse staircase. The method must keep the Newell solver
with full PBC.

## 2. Decisions and open points

**Decisions of Chris.**
- 2026-10-15: test A, a static script test of the correction (3.1). Stopped before it ran, because its fine reference
  carried the error of DemagFieldPBC (7.105).
- 2026-10-19: own branch from `main` for this work.
- 2026-10-19: only the Newell solver with full PBC. Never DemagFieldPBC, also not as a reference. The older state
  files (.pt) that were relaxed with DemagFieldPBC are not used.
- 2026-10-19: strategy approved: a static test of both forms of the correction (cell-local and pairwise) against a
  fine Newell reference with full PBC (the choice of V2 as the main variant was replaced on the same day, see below
  and 3.4).
- 2026-10-19: magnetisation states for the tests only from the Newell solver with full PBC and without Bloch points
  ("the cleanest states").
- 2026-10-19: V2 dropped (3.4). Main variant: the original method of Donahue & McMichael (V1: cut cells with volume
  fraction, near-field correction from a submesh); the static test computes both forms (pairwise and cell-local).
  Stage A of PLAN.md approved.
- 2026-10-19, night mandate: go as far as possible in the plan overnight; test everything that can be tested; for
  each finding give the options and then go on with the own recommendation, and record it here ("decision under the
  night mandate"). Rented machine: a cheap instance with CPU and GPU (about 0.2 … 0.3 $/h; an A100 offer at 0.59 $/h
  was destroyed unused on the order of Chris). Go for the end: when all is done and evaluated and the method is
  clear, download what is needed (md5) and destroy the instance.
- 2026-10-19: own record for the staircase correction (this file).
- 2026-10-20 morning: option B of the dynamic test (transferred start, 3.11). Both boxes keep running until D15g and
  T15g are done.
- 2026-10-20 morning: target of the correction for now: the staircase of dx 1.5 (q = 2), because it can be validated
  against dx 1.5 runs. The series in q towards the smooth sphere is not started.
- 2026-10-20: reading of the dynamic test: "reproduces" needs also 2 SE at most 0.15 (option B, 3.14). Box 1 is
  destroyed after the download of its data with md5 check.
- 2026-10-20: the centroid term (first moment of the cut cells in the far field) is part of the method; the method is
  V1 cut cells with the pairwise window term and the centroid term (3.12). Focus on q = 2 as the proof of concept. The sections 7.117, 7.120 and 7.121 of the study
  record are copied here as 3.1, 3.2 and 3.3, with the references to each other changed to the new numbers.

**Earlier decisions that bind this work.**
- 7.29: the volume-fraction model (M_s, A and K of the surface cells scaled by their volume fraction) is not usable
  at dx 3 for high fields; the staircase stays; no new smoothing variant without validation (decision of Chris
  2026-10-06; the exclusion confirmed 2026-10-15, 7.104).
- 7.110: whether the voxel staircase is part of the model or an error against the smooth sphere is a model question
  for Chris after the corner test (PLAN 4.13 of the β(B_peak) study; corner runs in 7.111 … 7.119). Also open there:
  the size and roughness question (at a fixed dx 3 l_ex, dx/d is 0.010 at d/l_ex 300 and 0.014 at d/l_ex 212).

**Open for Chris.**
- Stage B of PLAN.md (separate go after stage A).
- Target of the correction: the fine staircase (exact reference on a fine mesh) or the smooth sphere (only through a
  series of finer submeshes).
- Meng et al. 2025, J. Magn. Magn. Mater. 173350 (powder cores, cites Donahue & McMichael): full text wanted (3.2).

## 3. Record

### 3.1 Literature check: stray-field corrections of the staircase, and their use in 3D (search 2026-10-15, recorded 2026-10-19) (corrected in 3.3)

*First recorded as 7.117 in studies/beta_bpeak/PROTOKOLL.md (commit 1132876).*

**Source of this section.** The search was made on 2026-10-15. Chris asked first whether a correction of the
staircase can work, and then for a thorough search for 3D work that continues or develops the idea of Donahue. The
results were given in the chat. The notes file of the search was in a temporary folder; on 2026-10-19 it did not exist
any more. This section is written on 2026-10-19 from the chat record of that session. Details that are not in the
chat record are left out. Statements from other sections are marked with their section number.

**Terms.**
- l_ex: exchange length. d: particle diameter. Δx: edge of a cubic cell of the mesh. B0b: the mesh test of 7.102
  (d/l_ex 96, φ 0.65). FFT: fast Fourier transform. PBC: periodic boundary conditions.
- Staircase: the surface of a sphere that is built from cubic cells (7.104). E_st: the extra stray-field energy of the
  staircase against the smooth sphere (7.104).
- Partial cell: a cell that the true surface cuts. Volume fraction: the magnetic part of the volume of a cell.
- Submesh: a finer mesh inside a cell, with q sub-cells along each axis.
- Near-field correction: a change of the stray-field coupling between a cell and its neighbours up to p cells away.
  The change is computed once from a submesh.

**Method.**
- Web search, and the citation lists of the three base papers (OpenAlex, Semantic Scholar).
- The citers of Lepadatu 2015, Lebecki 2008 and Wei 2009.
- Read, as the chat record says: Donahue & McMichael 2007, García-Cervera et al. 2003, and Parker et al. 2000 (all 9
  pages; PDF from Chris). For the other entries the record gives the sources (list at the end). It does not say for
  each entry whether the full text or only the abstract was read.
- Checked: the source code of Boris and of magnum.np, and the OOMMF user guide.
- Full text not read: Rajaram et al. 2013 and 2014, Wei 2009, Wang 2011, GPMagnet 2012.
- Gap: the citation lists are not complete. Riahi et al. 2013 was missing in OpenAlex and Semantic Scholar. Only the
  web search found it.

**Base methods.**
- Donahue & McMichael, IEEE Trans. Magn. 43, 2878 (2007).
  - Only the stray field is corrected. The far field stays the FFT convolution. Each partial cell enters it with its
    volume fraction. Each partial cell also gets a near-field correction. A submesh of p × p neighbour cells with
    q × q sub-cells each gives it.
  - In the simplified form the correction is a local 3 × 3 tensor per cell. It acts like a cell anisotropy.
  - Tests only on films one cell thick: a NiFe square, rotated against the mesh. The frequency of the edge mode shifts
    by up to 50 % without the correction and by less than 5 % with it. The field that pushes a vortex out of the
    square changes with the rotation angle by 26 % without and by 3.2 % with the correction. This is a factor of
    about 8 … 10, for these two quantities and not for the energy. Δx/d of these tests is 0.011 … 0.014 (the range of
    this study).
  - The volume fraction alone ("rough correction only" in their figure) gives almost no improvement. The near-field
    term gives the effect.
  - The authors needed smaller cells to resolve the vortex core and to prevent its pinning on the mesh.
  - A talk slide (2006, slide 7; not in the paper) shows the error without the correction as first order in Δx in the
    cells at the edge (given as "0.003 Δ"), and as second order 5 nm inside. This is read from a slide. On 2026-10-15
    it was taken as support for the exponent 1.23 of 7.104; that exponent uses DemagFieldPBC values (7.105).
- García-Cervera, Gimbutas & E, J. Comput. Phys. 184, 37 (2003).
  - The staircase gives large errors in the exchange and in the stray field. In a hexagonal film the switching
    mechanism depends on the rotation against the mesh, with artificial vortices at the edge. A mesh with twice the
    resolution does not remove this.
  - Exchange: the Laplace operator at the surface is second order. The normal derivative of M is zero on the true
    surface; ghost values outside give this.
  - Stray field: the nearest neighbours use exact tensors of the true cell shapes (numerical integration). The far
    field uses the FFT with equivalent rectangular cells.
  - With both corrections the result does not depend on the rotation. The loops on 100² and 200² cells agree.
  - Only 2D films.
- Parker, Cerjan & Hewett, J. Magn. Magn. Mater. (2000), "embedded curve boundary" (ECB). Donahue 2007 cites it as
  ref. [4].
  - The stray field comes from a finite-difference Poisson solver for the potential (alternating-direction implicit
    method), not from an FFT convolution. This solver type does not fit the FFT solver of magnum.np.
  - At the nodes next to the surface the stencil uses the true distance δ to the surface and the jump condition with
    the charge of the smooth surface. The differences there are first order.
  - Exchange: the normal derivative of M is zero on the true surface. M is interpolated at a point 2δ inside, along
    the normal. The chat record names a risk at contact points between particles: there the normal is not defined.
    In this model the particles do not touch (geometry.py keeps at least one void cell between them).
  - This is the only 3D LLG implementation of the three base papers. But only the vertical side walls are curved.
  - Tests: a uniformly magnetised disk (r 200 nm, h 50 nm, Δx 10 nm). The field error at the edge is at most 10 %
    (x axis) and 7 % (y axis) with ECB, and at least 20 % on the y axis with the staircase. On a finer mesh both errors
    fall, and their ratio stays "roughly the same". The paper gives no convergence order and no energy. The field error
    falls by a factor of about 2 … 3 (estimate, read from the figures 3 and 4). µMAG standard problem 2 has no curved
    surface and is no test of ECB.
  - Co ellipses (450 × 250 nm, 16 … 40 nm thick, Δx 6 nm): single vortex in the experiment 97 %, with ECB 50 %, with
    the staircase only in the thickest ellipse. The coercive field is similar with ECB and with the staircase. Thus
    the staircase changes the state more than the coercive field. This agrees with the reading of 7.104 that coarse
    states follow the steps; 7.105 leaves that reading open.

**Use in 3D after 2007 (works found).**

| work | 3D | method | validation |
|---|---|---|---|
| Donahue & McMichael, talks 2006 and 2007 | no | near-field correction, simplified to a local tensor per cell; exchange with a fixed l_ex at volume-weighted M_s | films one cell thick |
| Lebecki, Donahue & Gutowski, J. Phys. D (2008) | 3D code with PBC | Donahue correction for an infinitely long rod (refinement 64 × 64 × 1, range 7 cells) | nucleation fields "fit quite well" to the analytic values; no values with and without the correction; only the 2D cross-section is curved |
| Rajaram et al., IEEE Trans. Magn. (2013), disks 8 … 20 nm thick | not known | uses a "3D-compatible correction" after Donahue | full text not read |
| Lepadatu, J. Appl. Phys. 118, 243908 (2015), and the code Boris (GPL-3) | yes | module "Roughness": a correction from a fine mesh, refined along all three axes, averaged onto the coarse cells; at run time a local symmetric 3 × 3 tensor per cell (seen in the source code); PBC supported | roughness up to a depth of about l_ex: "excellent agreement" with an explicit fine mesh (switching, wall mobility); as a staircase correction only on an 80 nm disk, qualitative |
| mumax3, EdgeSmooth | yes | scales only M_s with the volume fraction of the cell | similar to the volume-fraction model of 7.29, which scales M_s, A and K |
| OOMMF 2.1, magnum.np | – | no staircase correction (user guide and source code checked) | – |
| Wei & Wang (2009), Wang (2011) | no | triangular prisms as surface cells | films only |
| Riahi et al. (2013) | no | a variable moment per cell | 2D only (kagome lattice) |
| Lai, Li & Lin, J. Comput. Appl. Math. (2006) | yes | 3D Poisson solver with an "immersed interface" and FFT, second order, written for magnetostatics | analytic test cases only; not used in a micromagnetic code |
| reviews: García-Cervera (2007), Abert (2019), Miltat & Donahue (2007) | – | no 3D continuation; Miltat & Donahue: "practical solutions are expected in the near future" | – |

Gillis et al., J. Comput. Phys. (2018) is in the source list of the search; its content is not in the chat record.
Holt et al., arXiv 2412.10466 (2024) treat the anisotropy that the difference formulas of the exchange and the DMI
make, not the staircase.

**Result.**
- The search found no work that tests a staircase correction on a body that is curved in all three directions
  (sphere, ellipsoid). It found no convergence order and no energy of a corrected 3D scheme. The citation lists are
  not complete (see Method).
- Boris is the only public 3D implementation of the idea of Donahue that the search found.
- The statement of 2026-10-15 before the 3D search ("the correction probably works for the stray field") was
  withdrawn. The question is open.
- Limit 1, the simplified form. The correction of a cell then depends only on M of this cell. This is correct only
  if M changes little over the correction window (a few cells). Lepadatu checked this only up to a step depth of about
  l_ex. The steps of this study are Δx deep, that is 1.5 … 3 l_ex.
- Limit 2, the target. The coarse energy plus the correction gives the energy of the fine staircase of the submesh,
  not of the smooth sphere (stated in Boris). Boris computes the fine mesh over the whole box. The local windows of
  Donahue avoid this.
- Exchange: the staircase also gives a first-order error of the exchange energy (helix test of 7.104: −10 % on
  3 l_ex; a code test on the geometry d/l_ex 48, φ 0.45, not valid for B0b). The plan of 2026-10-15 was to correct
  the stray field first and the exchange only if necessary.
- Estimates of 2026-10-15 that are not valid now. They assume that the correction is otherwise exact. They come from
  the power law of E_st in 7.104 (exponent 1.23) and from the ratio of E_st to the anisotropy energy E_K. Both use
  DemagFieldPBC values, and 7.105 requires a new determination.
  - Need: E_st must fall to 10 % of |E_K|, a factor 5 … 46 over the design range and 19 for B0b on 3 l_ex.
  - Necessary refinement: about 11 for B0b on 3 l_ex, about 22 in the worst part of the design range. With Boris this
    is about 1.5·10⁸ fine cells for B0b, only for the start.
- Link to 7.29: the volume fraction alone is the far-field part of the method of Donahue. Donahue shows that it alone
  gives almost no improvement. 7.29 measured more loss at 100 mT with the volume-fraction model and left the cause
  open. The literature does not explain that result.
- FEM was weighed on 2026-10-15: exact surface, but the vortex cores must still be resolved (about 10⁸ elements per
  1 µm particle, estimate), PBC for the stray field are more difficult, and it is a new code stack.

**Decision of Chris on 2026-10-15 and what followed.**
- Options: A, a static script test of both forms of the correction (local tensor as in Boris; full local sum as in
  Donahue, eq. 5) on the B0b states on 3 l_ex, with a fine mesh of 0.75 l_ex. The reference was the same coarse M on
  the fine shape. Reading fixed before the test: a form is usable if its error against the fine reference is at most
  10 % of |E_K|. B, Boris as a test bench. C, an immersed-interface Poisson solver. D, no correction (grid-aligned
  particles or FEM).
- Decision of Chris: A.
- From 7.105: test A was stopped before it ran, because its reference on the fine mesh carried the DemagFieldPBC
  error. Since 7.105 the production solver is the Newell solver with the compensated outer surface charges (PLAN
  4.11 … 4.13 of the β(B_peak) study). check_demag_alias.py gives an independent reference. Thus the obstacle of test A is removed in
  principle.

**Sources.**
- Donahue & McMichael 2007: https://math.nist.gov/~MDonahue/pubs/jaggies_20070117.pdf; talks:
  https://math.nist.gov/~MDonahue/talks/bonn2006-staircase-errata20070116.pdf,
  https://math.nist.gov/~MDonahue/talks/joint2007-GE05-staircase.pdf
- García-Cervera, Gimbutas & E 2003: https://web.math.princeton.edu/~weinan/pdf%20files/accurate.pdf
- Parker, Cerjan & Hewett 2000: https://www.sciencedirect.com/science/article/abs/pii/S0304885300000433
- Lebecki, Donahue & Gutowski 2008: https://math.nist.gov/~MDonahue/pubs/pbc-in-oommf_4.6.pdf
- Lepadatu 2015: https://arxiv.org/pdf/1508.04391; Boris 2020: https://arxiv.org/pdf/2008.02478; code:
  https://github.com/SerbanL/Boris2, https://github.com/SerbanL/BORIS
- Rajaram et al. 2013: https://ieeexplore.ieee.org/document/6559077/
- mumax3: https://mumax.github.io/api.html, https://github.com/mumax/3/blob/master/engine/geom.go,
  https://arxiv.org/pdf/1406.7635
- OOMMF 2.1a1 user guide: https://math.nist.gov/oommf/doc/userguide21a1/userguide.pdf
- Wei & Wang 2009: https://ieeexplore.ieee.org/document/5170238/; Wang 2011: https://ieeexplore.ieee.org/document/6028212/
- Riahi et al. 2013: https://arxiv.org/abs/1401.5274
- Lai, Li & Lin 2006: https://xblin.math.ncsu.edu/preprint/lai_li_lin.pdf; Gillis et al. 2018:
  https://www.sciencedirect.com/science/article/abs/pii/S0021999117308094
- Holt et al. 2024: https://arxiv.org/abs/2412.10466
- Miltat & Donahue 2007: https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=50910; García-Cervera review 2007:
  https://web.math.ucsb.edu/~cgarcia/Papers/Garcia-Cervera_ReviewSEMA.pdf; Abert 2019: https://arxiv.org/pdf/1810.12365

### 3.2 Second literature search: volume-preserving boundary cells (V2), and corrections that keep the Newell solver with full PBC (2026-10-19, questions of Chris) (corrected in 3.3)

*First recorded as 7.120 in studies/beta_bpeak/PROTOKOLL.md (commit a3dac8d, restored in 76dfcc5).*

**Questions of Chris.** (1) A targeted search for the variant V2 (below), because the idea is obvious. (2) A new
survey of the variants for our goal that stay compatible with the Newell solver with full periodic boundary
conditions (PBC). Rule of Chris for this work: only the Newell solver with full PBC, never DemagFieldPBC, also not as
a reference.

**Terms.**
- Newell solver with full PBC: DemagField(true_pbc=True), the FFT convolution (FFT: fast Fourier transform) of M with
  the Newell tensor of cubic cells, with image shells and the compensation of the outer surface charges (7.105).
- m: unit vector of the magnetisation; M_s: saturation magnetisation; Δx: cell edge; f: volume fraction (3.1).
  Staircase, partial cell, submesh, near-field correction: as in 3.1. Window: the neighbour cells, up to p cells
  away, for which the near-field correction is applied.
- V1 (Donahue & McMichael 2007): a partial cell is the part of its cube inside the body, with the moment
  f·M_s·Δx³. Cells with a centre outside the body but with f > 0 become magnetic. The near field comes from a fine
  submesh, that is, from a fine staircase (3.1, limit 2). García-Cervera et al. 2003 use the exact cut shapes
  instead (2D only).
- V2 (idea of Chris): each voxel cell keeps its volume Δx³ and its moment M_s·Δx³. Only the shape of the boundary
  cells changes, with the volume kept, so that all cells together fill the smooth sphere with exactly the voxel
  volume. The near field comes from these deformed cells. The FFT far field does not change.

**Physics of V2 (own reasoning, 2026-10-19).**
- One magnetisation field M(r) serves all energy terms.
- Each deformed cell is a closed body. Thus each cell is charge-neutral.
- If m is equal in two neighbour cells, their common border carries no charge. Thus the surface charge is M·n on the
  smooth surface.
- The stray-field energy of the deformed cells is the exact energy of this M(r). Thus it is symmetric and not
  negative.
- Two conditions limit this in the computation. (a) Outside the window the coupling of the cubes is used. If the
  centroid of a deformed cell is not at the centre of its cube, this gives an error at quadrupole order outside the
  window. (b) With a cut window the energy is exact only up to the cut, and the sign of the energy is not guaranteed
  by construction. The static test must measure both, and the correction matrix must be symmetric.
- V2 leaves the exchange on the cube neighbours. The first-order exchange error at the surface stays (helix test of
  7.104, −10 % on 3 l_ex, a code test on another geometry, see 3.1).
- Why the simple form of the idea is not used: put a charge with the true normal on the free face of a full cube and
  keep the other faces. The cell then has a net charge. It is zero only if the new charge face has the same vector
  area as the free cube faces, that is, only if the normal does not change. The particle as a whole stays neutral, but
  each cell has a monopole when m varies.

**Method.** Two literature searches in parallel on 2026-10-19 (web, arXiv, citation lists from OpenAlex and in part
from Semantic Scholar, which stopped at its rate limit). The notes give each query, each hit with its citation and
link, and whether the full text or only the abstract was read: studies/staircase_correction/lit/
LIT_v2_volume_preserving.md and LIT_newell_pbc_variants.md (branch feature/staircase-correction, commit 7158436). After
the searches, three sources were checked again on arXiv for title, authors and date (Ai, Duan & Lomakin 2024;
Durhuus, Insinga & Bjørk 2026; Adams & Michels 2026).

**Result 1: V2 was not found.**
- Fields searched: micromagnetics; the discrete dipole approximation (DDA: light scattering on a cubic lattice of
  dipoles, with FFT); conformal FDTD and cut-cell methods; magnetic and gravity forward modelling with rectangular
  prisms (geophysics; prism formulas equivalent to Newell); particle-mesh methods (P3M, precorrected FFT); partitions
  of a body into cells of prescribed volume. No source keeps the cell volume and the cell moment and changes only the
  shape of the boundary cells.
- Closest relatives:
  - Donahue & McMichael 2007 (full text): their eq. 2 averages the fine tensor over the sub-cells that belong to a
    coarse cell. The search infers that the same formula gives V2 with a transport assignment of the sub-cells and the
    volume Δx³. This inference must be checked in the code.
  - Precorrected FFT (Phillips & White 1994, full text) and the adaptive integral method: the same algebra (exact near
    field minus the grid near field).
  - Lebecki, Donahue & Gutowski 2008 (full text): the correction of Donahue together with a periodic kernel (1D PBC).
    The only finite-difference micromagnetic precedent found for a surface correction under PBC.
  - Partition with prescribed volumes: its existence and its link to the least-squares transport are known
    (Aurenhammer et al. 1992/1998, abstract; Lévy 2015, full text). The search found no use for the cells of a field
    kernel.
  - Sen, Wang, Breitzman & Dayal 2024 and 2026 (abstracts only): in the continuum limit of a lattice of charge-neutral
    cells, a surface charge appears that depends on the partial cells at the boundary. That this supports V2 at fixed
    moments is an inference of the search.
- No source compares cells with constant volume against cells with a volume fraction. The sources compare the volume
  fraction alone with the volume fraction plus a near-field correction. In Donahue & McMichael the near-field term
  gives the effect (3.1). García-Cervera et al. 2003 claim large errors of the volume weighting in the abstract, but
  show no curve with the volume weighting alone. Riahi et al. 2013 (full text) find that a volume-weighted mask on
  5 nm cells gives about the same coupling energies of curved magnets as the staircase on 2 nm cells. Thus the volume
  fraction helps for energies of the far field.
- Risks of V2 (inferences of the search and own): the partition must stay in a shell at the surface, else inner cells
  deform a little and the window must include them; the exact coupling of the deformed cells needs new code (own
  plan: sums of sub-cell Newell tensors).

**Result 2: variants that keep the Newell solver with full PBC.**

| family | what it corrects | evidence | fits Newell with full PBC | main risk for us |
|---|---|---|---|---|
| F1 pairwise near-field correction (V1, V2) | near field of the stray field | Donahue & McMichael 2007 and talks (2D films); DDA weighted discretisation (Piller 1997, abstract; Yurkin et al. 2006, full text; 3D spheres in optics) | yes, as a sparse local term | memory of the pair matrix (0.1 … 2 GB, estimate of the search); the matrix must be symmetric |
| F1 in the cell-local form (Boris, one tensor per cell) | near field, only for m uniform over the window | Donahue talk 2006, thin spheroid: RMS field error 118 % with the staircase, 29 % with cell-local corrections only | yes | DDA: the weighted discretisation keeps a first-order term with a smaller constant (Yurkin et al. 2006). That this holds for the cell-local form in micromagnetics is an inference. |
| F2 refined grids with FFT (precorrected FFT, NUFFT, adaptive refinement) | geometry of the sources | Ai, Duan & Lomakin 2024 (finite elements, FEM, periodic kernel); Exl & Schrefl 2014 | precorrected FFT yes; adaptive refinement only with large changes | a second solver, weeks of work (estimate of the search) |
| F3 charge correction (smooth minus staircase surface charge) | near field | no source in finite-difference micromagnetics | yes; the search infers that it is F1 in another basis | needs M on the true surface: extrapolation over one cell |
| F4 / F8 cell weights and composite voxels (DDA, FFT homogenisation) | effective properties of boundary cells | DDA: integral quantities up to 24 times (Piller 1997) and 60 times (Zhu et al. 2023) more accurate, abstracts only; composite voxels: about 10 times smaller error (Lendvai & Schneider 2024, conference paper on thermal conductivity, one sphere of 17 % volume, contrast 10) | yes | rules for linear media; for LLG they reduce to volume fractions plus a cell-local term of the Donahue type (inference of the search); 7.29 tested the volume fraction alone |
| F6 exchange at the true surface | exchange | García-Cervera 2003, Parker 2000 (3.1); no new source | yes, local | second part of a full correction |
| F7 averaging over grid offsets or rotations; Richardson extrapolation at a locked geometry | the random part of the shape error | no source in DDA or micromagnetics; Yurkin et al. 2006 show that the shape error is quasi-random in the cell size | yes, no solver change | the systematic part stays: the staircase area is 1.5 times the sphere area for each Δx (own calculation of the search; 7.29 measured 1.50); at a locked geometry the extrapolation converges to the staircase body |

- Under PBC the near-field correction uses the free-space kernel with the minimum-image distance (the window wraps
  across the box edge). Precedents: Lebecki et al. 2008 (finite differences, 1D PBC) and Ai, Duan & Lomakin 2024
  (FEM; section 3.2 of that paper, read by the search: a free-space near correction together with a periodic FFT kernel). That this
  also holds with the Newell kernel is an inference.
- Durhuus, Insinga & Bjørk 2026 (abstract and start of the text read): for sufficiently large samples under 3D PBC,
  only the mean magnetisation carries the effect of the outer sample shape. This agrees with the compensation term of
  DemagField(true_pbc=True) (7.105).
- Holt et al. (3.1) is now published: IEEE Trans. Magn. 2025, doi 10.1109/TMAG.2025.3604525.

**Gaps.**
- No source tests a staircase correction on hysteresis or loss of spheres. The validations are field errors,
  energies, FMR frequencies and scattering cross-sections.
- All effect sizes of the sources (Donahue & McMichael, García-Cervera et al., Riahi et al.) come from single
  deterministic simulations without error bars.
- No finite-difference micromagnetic work combines a staircase correction with the FFT and full 3D PBC.
- Full text not read: Zhu, Liu & Yurkin 2023; Kabel et al. 2015; Piller 1997; Piller & Martin 1998; Tandon et al.
  2004; Wei 2009; Wang 2011; the geophysics and FDTD papers. The citation lists are from OpenAlex and in part from
  Semantic Scholar. Six citers of Donahue & McMichael that only Semantic Scholar lists are not read.
- Meng et al. 2025, J. Magn. Magn. Mater. 173350 (powder cores, 1 … 20 kHz, cites Donahue & McMichael): abstract not
  found. It has the same application as this study; a check by hand is proposed.

**Consequence for the strategy (proposal; decision of Chris open).**
- The pairwise near-field correction, not the cell-local tensor (Donahue: 29 % field error remains with cell-local
  corrections only; DDA: a first-order rest; transfer to our case is an inference).
- V2 as the main variant. It changes only the stray-field near field and leaves the cells, moments, exchange,
  anisotropy and Zeeman terms as they are. V1 only as a check of the computation in the static test.

**Sources** (full citations and read status in the two notes files).
- Ai, Duan & Lomakin 2024: https://arxiv.org/abs/2409.13958
- Durhuus, Insinga & Bjørk 2026: https://arxiv.org/abs/2604.08777
- Lebecki, Donahue & Gutowski 2008: doi 10.1088/0022-3727/41/17/175005
- Phillips & White 1994: https://cecs.uci.edu/~papers/compendium94-03/papers/1994/iccad94/pdffiles/04c_2.pdf;
  journal version 1997: doi 10.1109/43.662670
- Aurenhammer, Hoffmann & Aronov: doi 10.1145/142675.142747 (1992), doi 10.1007/PL00009187 (1998); Lévy 2015:
  doi 10.1051/m2an/2015055
- Sen, Wang, Breitzman & Dayal: doi 10.1016/j.jmps.2024.105743 (2024), doi 10.1137/25M1747907 (2026)
- Yurkin, Maltsev & Hoekstra 2006 (I): doi 10.1364/JOSAA.23.002578, arXiv 0704.0033
- Piller 1997: doi 10.1364/OL.22.001674; Zhu, Liu & Yurkin 2023: doi 10.1364/OE.509479
- Lendvai & Schneider 2024: https://scipedia.com/wd/images/b/b0/Draft_Sanchez_Pinedo_909610227pap_645.pdf
- Exl & Schrefl 2014: doi 10.1016/j.jcp.2014.04.013
- Riahi et al. 2013: doi 10.1088/0953-8984/25/29/296001
- Holt et al. 2025: doi 10.1109/TMAG.2025.3604525

### 3.3 Correction of 3.1 and 3.2: Donahue & McMichael use the cell-local form for all their results (2026-10-19, question of Chris)

*First recorded as 7.121 in studies/beta_bpeak/PROTOKOLL.md (commit 76dfcc5).*

**Question of Chris.** Does Donahue use the cell-local tensor or the full pairwise correction?

**Source.** Full text of Donahue & McMichael, IEEE Trans. Magn. 43, 2878 (2007),
https://math.nist.gov/~MDonahue/pubs/jaggies_20070117.pdf, text extracted with pdftotext on 2026-10-19.

**What the paper says.**
- Terms of the paper: N_ij is the coupling of the coarse cells i and j (far-field approximation). L_ij is the
  coupling from the fine submesh, averaged over the sub-cells inside the body. v_i is the volume fraction of cell i.
  Ω_local,i is the window around cell i.
- Eq. 5 is the full pairwise correction: the field correction of cell i is the sum over the cells j in the window of
  the difference of the two couplings times M_j.
- Eq. 7 splits this sum into a part with the difference m_i − m_j and a part with m_i alone. For an exchange-coupled
  material the first part is small. At the edge of the window it is also damped, because L_ij approaches N_ij.
- Eqs. 8 and 9 keep only the second part. The correction then acts like a cell-local anisotropy: a tensor K_i per
  cell, the sum over the window of (v_j / v_i)(L_ij − N_ij). The paper calls this "essentially a generalization of
  the approximation used in [7] to model film thickness variation" ([7]: Porter & Donahue, J. Appl. Phys. 89, 7257,
  2001).
- Quote: "This is the form for H_corr,i that is used in the remainder of the paper." And: "However, if higher accuracy
  is required, then (5) may be used directly, or else a higher order approximation to m_j may be introduced into (7)."
- Windows in the figures: a p × p neighbourhood with q × q sub-cells, shown as 7 × 11, 7 × 25 and 7 × 51. Quote:
  "Local neighborhood size p was found to be less important."
- The results of 3.1 are results of the cell-local form: the edge-mode frequency rises by about 50 % with the
  rotation for the uncorrected edge and also for the diluted edge cells alone, and changes by less than 5 % with the
  cell-local correction; the vortex exit field changes by 26 % without and by 3.2 % with it. The caption of their
  Fig. 2 calls it "the full correction using the locally refined grid", in contrast to the diluted edge cells alone.
  The paper does not compare the cell-local form with eq. 5.
- The paper itself limits the edge-mode test: the spins of an edge mode stay nearly uniform, so this test "does not
  strongly test the assumption" of eq. 8. The vortex exit field is the stronger test. It used 2.5 nm cells in a
  NiFe film 2.5 nm thick. With A 13 pJ/m and M_s 800 kA/m from the paper, l_ex is about 5.7 nm, so the cells are about
  0.44 l_ex (own calculation). The cells of this study are 1.5 … 3 l_ex.

**Corrections.**
- 3.1, option A of 2026-10-15: "full local sum as in Donahue, eq. 5" names the right equation, but the results of
  Donahue & McMichael use the cell-local tensor of eqs. 8 and 9. This is the same form as in Boris.
- 3.2 (variant V2: boundary cells with constant volume), table rows "F1 pairwise near-field correction"
  and "F1 in the cell-local form", and the first bullet of the consequence: the value "RMS field error 118 % with the
  staircase, 29 % with cell-local corrections only" comes from the talk of Donahue (2006, thin oblate spheroid of
  10 × 10 × 1 cells). It is a cell-local thickness correction, the thin-film precursor of K_i, not K_i of eq. 9 itself.
  It is therefore no evidence against the cell-local tensor of the 2007 paper. The 2D results of the 2007 paper are
  evidence for the cell-local tensor. The micromagnetic sources found give no measured result for the full pairwise
  form.
- Cause of the error: the final report of the search, given in the chat, wrote "with only a cell-local term" and
  dropped the word "thickness". The notes (LIT_newell_pbc_variants.md) have it. The author of 3.2 took the
  statement from the report without a check in the notes or in the source.

**What remains.**
- Whether the cell-local tensor is sufficient for this study is open. The arguments for the pairwise form are only
  inferences: (a) eq. 8 assumes that m changes little over the window; the paper tested this with cells of about
  0.44 l_ex, and the cells of this study are 3 … 7 times larger relative to l_ex; whether the coarse states follow
  the steps is open (7.105); (b) in the discrete dipole approximation (DDA, 3.2) the weighted discretisation, a
  cell-local rule, keeps a first-order term (Yurkin et al. 2006; an analogy).
- Consequence for the strategy (proposal): the static test measures both forms, as test A of 2026-10-15 planned. The
  cell-local tensor needs 6 numbers per cell. The pair matrix needs up to about 2 GB (estimate of the search, 3.2).
  The same derivation holds for V2 with the volume Δx³ instead of v_i.

### 3.4 V2 cannot be built locally; decision: the original method of Donahue & McMichael (2026-10-19)

**Finding.** A context-free review of the first PLAN draft built the fine body and the V2 owner rule for one particle
of the production geometry, and an independent check repeated it (scratchpad script v2match.py; maximum bipartite
matching, scipy).
- K2 geometry (d/l_ex 300, φ 0.62, sphere radius 50.65 cells, q = 2): 544 840 voxel cells, fine body 4 358 152 fine
  cells (q³ times the voxel cells: 4 358 720, δ_p −568), 28 256 orphan fine cells, 28 824 free places in 12 432
  receiving cells.
- With the best possible assignment inside 3 Δx, 744 orphans have no receiving cell (script v2match.py in this
  folder). All 744 lie at the six poles of
  the sphere: the largest component of their normal is above 0.99. The review found 312 such orphans in the K4
  geometry (best assignment, 3 Δx), and with its greedy rule still 763 in K2 at 10 Δx.
- Cause: at each pole the voxel sphere has a wide flat terrace (half-width about 10 cells at d/l_ex 300, estimate from
  the cap height of one cell). Above the centre of the terrace the true sphere has more volume than the voxel cells;
  at the rim of the terrace the voxel cells have more volume than the sphere. Whole cells cannot represent a cap that
  is less than one cell high. Thus the volume error of the voxel model is spread over the width of the terrace, not
  over one cell. Cells with a fixed volume Δx³ can follow the sphere there only if volume moves over about 5 … 10
  cells (estimate from the terrace width).
- Consequence: V2 does not violate a principle, but it is not local at the poles. Options given to Chris: A, V2 with a
  global partition (inner cells deform; the magnetisation moves over several cells); B, V1 (cut cells with volume
  fraction, the original method); C, V2 with volume fractions only at the pole terraces; D, no correction.

**Decision of Chris (2026-10-19): B.** The original method of Donahue & McMichael is the main variant. The static test
computes both forms (pairwise and cell-local). The volume fraction is here part of the representation of the true
sphere, not a smoothing; whether it changes the dynamics as in 7.29 (there without a near-field correction and with
DemagFieldPBC) is a question for phase 3.

**Further findings of the review, applied to PLAN.md.** The sign convention of magnum.np is H = K M (self term about
−1/3), so the correction enters with a plus sign. The fine kernel with the production setting takes about 1.6 … 1.7 h
per geometry on this computer (measured kernel speed, extrapolated), so the static test uses the production setting
and needs no one-shell deviation. The sign check uses the smallest eigenvalue (Lanczos) instead of random fields.
Donahue's "7 × 7" neighbourhood is p = 3 in the notation of the plan. The field layer and the reading for uniform m
are defined more exactly.

### 3.5 Night work 2026-10-19/20: machine, order of work, decisions under the night mandate

**Machine.** vast.ai instance 56790356 (RTX 3090 24 GB, 27 CPU cores as cgroup quota, 229 GB RAM as cgroup limit,
0.27 $/h plus disk), rented by this session; code from GitHub, branch feature/staircase-correction (clone in
/root/msc). Unit tests T0 on the box: 6 passed. The state files were copied with md5 check. All computations of this
work run on the box (Chris: "nicht mehr lokal auf dem Laptop"); the local stage A run was stopped.

**Order of work (decision under the night mandate).**
1. CPU: stage A of PLAN.md with q = 2 (K4 geometry and K2 geometry in parallel), then the same with q = 4 on the
   production geometry (the RAM of the box allows the fine mesh of about 600³ cells; PLAN.md had q = 4 only in stage
   B because of the local computer).
2. GPU, at the same time: first dynamic arms at the B0b point (section D1 below).
3. Phase 2 code (field term, run_loops option), then the corrected dynamic arms on the GPU.
4. Stage B (static, B0b geometry) with the prepared states of the GPU arms.

**D1 — first dynamic arms at the B0b point (decision under the night mandate).**
- Question: does the corrected coarse model reproduce the effect of a finer staircase on the loss?
- Options: (a) compare the corrected dx 3 runs with the existing N0b references (fine voxel geometry, 7.107 …
  7.110); (b) new paired arms with the production protocol of today (option C, Bloch-point check), the same seeds
  and settings; (c) only static tests tonight.
  - (a): no new GPU time for the references; but N0b used the virgin relaxation without Bloch-point control, 4 shells
    without the extrapolation and other cycle settings, so the comparison mixes protocols; the N0b states are on the
    study instance, which this session does not touch.
  - (b): one protocol for all arms, states without Bloch points (rule of Chris); costs GPU time (dx 1.5 runs about
    1 … 2 h each, estimate).
  - (c): no answer on the loss.
- Recommendation and decision: (b).
- Arms at the B0b point (d/l_ex 96, φ 0.65, L_eff 12 l_ex, r_p 0), 9 mT, seeds 941 … 944, production settings of
  PLAN 4.13 of the β(B_peak) study (newell_ext with 4 shells, relaxation to convergence, atol 1e-6, cycle tolerances
  1e-7 / 1e-6, 7 … 15 cycles, at least 2 calm cycles, drive correction b_tol 0.01), option C (random start,
  relaxation, AC pass, relaxation, Bloch-point check, at most 3 attempts). Job files in jobs/:
  - C3: dx 3, voxel model (today's production model);
  - C15: dx 1.5, voxel model of dx 1.5 (its staircase is the fine body of q = 2);
  - C15g: dx 1.5 with the voxel geometry of dx 3 locked (only the resolution of m changes, as V3g in 7.110);
  - later, after phase 2: dx 3 with the volume fraction alone (q = 2), and dx 3 with the correction (q = 2).
- AC amplitude: --predemag_h 0.023 = twice the drive amplitude of the 150 mT stage at the B0b point in the older dx 3
  runs (h_amp 0.0115; those runs used DemagFieldPBC). The amplitude only sets the start of the AC pass at about
  300 mT; the corner runs use the same rule (7.114).
- Planned reading (fixed now, before the results): geometry effect of the correction Δ_c = ln P(dx 3 corrected) −
  ln P(C3), geometry effect of the fine staircase Δ_g = ln P(C15) − ln P(C15g), both paired by seed (mean ± SE, n 4).
  "Correction reproduces the geometry effect" if |Δ_c − Δ_g| ≤ 2 SE of the difference and ≤ 0.15; "not reproduced"
  if the difference is larger than both; else "not determined". Limit: independent preparations per arm, so the
  scatter between states enters (7.108: SD of ln P over seeds 0.2 … 0.6).

**Stage A, first facts (box, q = 2).** Step 2 passed for all three states: total energy against the run record
−0.0003676531 / −0.0003676500 (K4_1), −0.0074936839 / −0.0074936824 (K2_0), −0.0074991682 / −0.0074991672 (K2_1);
|m| − 1 at most
1.6·10⁻⁷; 0 Bloch points. R0 passed for both geometries: K4 (φ 0.55): 2 089 600 voxel cells, 93 792 modified cells,
46 176 new cells, φ_vox 0.55041, φ_f 0.55010; K2 (φ 0.62): 2 179 360 voxel cells, 98 304 modified, 45 696 new,
φ_vox 0.62058, φ_f 0.61976. Thus the material setup of the static test is verified against the runs.

**D2 — contact of cut cells at the B0b point (finding, decision under the night mandate).**
- Finding (smoke run of run_loops with --staircase at the B0b point): geometry.fcc_box stops ("Particles touch on the
  voxel grid"). At d/l_ex 96 the gap between neighbour particles is 1.44 cells at φ 0.65 and 1.95 cells at φ 0.62.
  The cut cells reach about one cell beyond the voxel surface, so cut cells of two particles touch at edges and
  corners: 3 552 (q = 2) and 4 704 (q = 4) contacts in the 26-neighbourhood at φ 0.65; 864 and 2 688 at φ 0.62. No
  face contacts (thus no exchange between particles). At φ 0.55 (gap 3.20 cells) there are none, for q = 2 and 4.
  The production geometry of stage A (d/l_ex 300) has a gap of about 6 cells and passed R0.
- Options: (a) the dynamic test at d/l_ex 96 with φ 0.55 (in the design range; corners K4 use φ 0.55); (b) keep
  φ 0.65 and allow edge and corner contacts of cut cells for --staircase (change of a safety rule of the geometry;
  the arm "volume fraction alone" then has the coarse coupling at these contacts); (c) the dynamic test on the
  production geometry (d/l_ex 300): the dx 1.5 arms would take about 8 times longer, not possible tonight.
- Recommendation and decision: (a). No rule changes, all arms valid. The link to the N0b history of φ 0.65 is lost,
  but D1 does not use N0b. The arms C3 at φ 0.65 that already ran are information only; C15 and C15g at φ 0.65 are
  cancelled before their start.
- New arms (names D…, the same settings as in D1, φ 0.55, seeds 941 … 944): D3 (dx 3 voxel), D15 (dx 1.5 voxel), D15g
  (dx 1.5, geometry of dx 3 locked), then D3f (dx 3, volume fraction alone, q = 2), D3P / D3L (dx 3, pairwise /
  cell-local correction, q = 2, p = 3), and D3P4 / D3L4 (q = 4: the target is the voxel geometry of dx 0.75).
  AC amplitude --predemag_h 0.040, the value of the corner runs with φ 0.55 (0.0400 … 0.0401).

### 3.6 Stage A with q = 2: results; the window converges only as 1/p; decision D3 (2026-10-20, night mandate)

**Run.** Box 56790356, CPU, double precision, production stray-field setting (newell_ext, 4 shells) on both meshes,
fine meshes 304³ (K2) and 312³ (K4); check_static.py at commit 0bc5597 (the p = 4 … 6 series at b6f4b44). Files:
results/static_A_*.json (copied from the box), report results/static_A.md (eval_static.py).

**R0, R1: passed** (both geometries). For uniform m along (100), (110), (111): E_ref and E_vox agree with
φ (1 − φ) / 3 within 10⁻⁹ (relative); both forms of the correction agree with E_ref within 10⁻⁹ at every p.

**States** (energies per box volume in K_d; S = E_vox − E_ref):

| state | E_K | S | S / \|E_K\| | E_frac − E_ref | RMS h_vox − h_ref (layer) | RMS h_frac − h_ref (layer) |
|---|---|---|---|---|---|---|
| K2_0 | −1.042e-2 | +4.72e-4 | 0.045 | −2.11e-4 | 0.0333 | 0.0136 |
| K2_1 | −1.040e-2 | +4.72e-4 | 0.045 | −2.12e-4 | 0.0333 | 0.0135 |
| K4_1 | −1.847e-3 | +4.21e-4 | 0.228 | −1.89e-4 | 0.0338 | 0.0126 |

- The staircase raises the stray-field energy of these states by S. Relative to |E_K| this is 4.5 % for the K2 states
  (below the 10 % rule of 3.1) and 23 % for K4_1 (above it).
- The volume fraction alone overcorrects: E_frac lies below E_ref by about 0.45 S. Its field error in the layer is
  about 0.4 of that of the voxel model.

**Forms** (the same for the three states within 0.01; K2_0 shown; e_E = |E_form − E_ref| / |S|, e_h = field error
in the layer relative to the voxel model):

| p | P: e_E | P: e_h | L: e_E | L: e_h |
|---|---|---|---|---|
| 1 | 0.948 | 0.169 | 0.952 | 0.194 |
| 2 | 0.594 | 0.096 | 0.602 | 0.133 |
| 3 | 0.432 | 0.067 | 0.442 | 0.113 |
| 4 | 0.339 | – | 0.352 | – |

(p = 4 from the energies of the second run; its field values are in the json.) |E_form − E_ref| / |E_K| at p = 3:
0.020 (K2), 0.097 (K4_1).

- R2: no form is "sufficient": the field criterion holds for P from p = 2 on and for L nowhere at p ≤ 3; the
  energy criterion e_E ≤ 0.10 fails at all p ≤ 4.
- The energy error falls about as 1/(p + 0.6): 4.47, 2.80, 2.04, 1.60 · 10⁻⁴ K_d for p = 1 … 4 (K2_0); the ratio of
  steps is 0.63, 0.73, 0.78. Reading by the rule of PLAN.md: **"window not converged"**.
- Cause (inference, to be tested): outside the window each cut cell acts as a cube with the moment v M_s Δx³ at the
  centre of the cube; the true moment sits at the centroid of the body part of the cell. The error per pair then
  falls as 1/r⁴ (a misplaced dipole). Summed over the full cells of the particle interior this gives an energy error
  that falls only as 1/R with the window radius R. A window large enough for e_E ≤ 0.1 would need p of about 14
  (estimate from the law), which is not practical.

**D3 — how to go on (decision under the night mandate).**
- Options: (a) add the first moment of each cut cell to the far field (moment at the centroid): the FFT field gets
  the term of the centroid offsets d_i of the source and target cells, with the central difference of the same
  periodic kernel (no new kernel; about 4 … 5 times the FFT work of the stray field, estimate); the window term
  subtracts the same terms; then the pair error falls as 1/r⁵ and the energy error as 1/R² (expected, not shown);
  (b) take P at p = 3 on the field criterion (local fields in the layer within 7 %) and let the dynamic arms decide;
  (c) a larger window: not practical (p about 14).
- Recommendation and decision: (b) and (a) together. The dynamic arms with P and L at p = 3 run anyway (D1). A
  static prototype of (a) is built and tested tonight with the same states and readings; if it converges, the
  dynamic arms are repeated with it.

### 3.7 D3 (a): with the first moment of the cut cells in the far field both forms are sufficient from p = 2 on (2026-10-20, night mandate)

**Method** (commit 95002a4; unit tests 10 passed). centroid_offsets: d_i = offset of the centroid of the body
sub-cells of cell i from the cell centre (in cells). Far field with the moment of each cut cell at its centroid, to
first order:

    h_i = Σ_j v_j [K(r_ij) + ((d_i − d_j) · D) K(r_ij)] m_j

with D the central difference of the same periodic kernel (factor i sin ω in Fourier space; no new kernel, the image
sum and the compensation unchanged). The window term subtracts the same coarse coupling, so for a window over the
whole box the correction stays exact (unit test, open boundaries, maximum deviation below 10⁻¹⁰). The operator stays
symmetric with the weight v (unit test).

**Window series without the centroid term** (K2_0, p = 4 … 6, fine reference unchanged): energy error 1.60 / 1.31 /
1.11 · 10⁻⁴ K_d, e_E 0.34 / 0.28 / 0.24, e_h 0.052 / 0.043 / 0.036 (K2_1 the same within 0.003 · 10⁻⁴; the law
1/(p + 0.6) predicts 1.08 · 10⁻⁴ at p = 6). Thus a larger window alone
does not help.

**Stage A with the centroid term (q = 2)** (files results/static_A_*_q2_c.json; report results/static_A.md):

| geometry | form | p | e_E | \|ΔE\| / \|E_K\| | e_h (layer) |
|---|---|---|---|---|---|
| K2 (2 states, equal within 0.005) | P | 1 / 2 / 3 | 0.193 / 0.049 / 0.018 | 0.009 / 0.002 / 0.001 | 0.047 / 0.011 / 0.004 |
| K2 | L | 1 / 2 / 3 | 0.189 / 0.045 / 0.014 | 0.009 / 0.002 / 0.001 | 0.080 / 0.062 / 0.061 |
| K4_1 | P | 1 / 2 / 3 | 0.196 / 0.050 / 0.019 | 0.045 / 0.012 / 0.004 | 0.046 / 0.010 / 0.004 |
| K4_1 | L | 1 / 2 / 3 | 0.194 / 0.048 / 0.016 | 0.044 / 0.011 / 0.004 | 0.076 / 0.058 / 0.056 |

- R2: both forms "sufficient" from p = 2 on, in all three states (R3 not yet computed).
- The error now falls by a factor of about 4 from p = 1 to 2 and about 2.6 from p = 2 to 3 (before: 0.63 and 0.73).
  This agrees with the cause given in 3.6 (misplaced moments in the far field).
- The pairwise form reproduces the field in the layer within 0.4 % of the voxel error at p = 3; the cell-local form
  stays at about 6 % (its own approximation, m_j replaced by m_i). At the new cells the field error of P3 is 0.0011
  (volume fraction alone: 0.12).
- Volume fraction plus centroid term without the near-field correction: E − E_ref = −0.28 S (K4_1). Thus the window
  term is still necessary.

**D4 — dynamic arms (decision under the night mandate).**
- Options: (a) let the running arms without the centroid term (D3P, D3L, D3P4, D3L4) finish; (b) stop them and run
  the arms with the centroid term instead; (c) both.
- (a) tests a variant that the static test has shown to be not converged; (c) needs about twice the GPU time and the
  references D15 / D15g would not finish tonight.
- Recommendation and decision: (b). The centroid term goes into the field term of magnum.np and into run_loops
  (--staircase_centroid). New arms: D3Pc and D3Lc (q = 2, p = 3), then D3Pc4 and D3Lc4 (q = 4, p = 3; the target is
  the voxel geometry of dx 0.75, information). D3f (volume fraction alone) stays as a control. The three arms without
  the centroid term that already ran or run are information (n = 1 each).

**D5 — one image shell for the static tests with q = 4 (decision under the night mandate).**
- Check (K2, q = 2, centroid term, p = 2 and 3; file results/static_A_d300_phi0.62_q2_k1_c.json): the plain Newell
  kernel with one image shell (k = 1) on both meshes against the production setting (newell_ext, 4 shells). R1:
  passed (largest deviation from φ (1 − φ)/3: 1.2·10⁻⁵, at m along (111)). S changes by 3·10⁻⁸ K_d (relative
  6·10⁻⁵); e_E and e_h of both forms agree to 4 digits (for example P3: e_E 0.0182 / 0.0182, e_h 0.0038 / 0.0038).
  Thus the image sum cancels in the differences, as argued in PLAN.md.
- Options: (a) q = 4 with the production setting (fine meshes 608³ and 624³; the fine kernel would take about
  6 … 8 h, estimate); (b) q = 4 with k = 1 (fine kernel in minutes); (c) no q = 4 on the production geometry.
- Recommendation and decision: (b). Stage A with q = 4 (K2 and K4, centroid term, p = 1 … 3, k = 1) started.
- Engineering change before (commit 41df994): the couplings L_ij are computed directly in blocks of source cells
  instead of with a table of occupancy patterns, which would need tens of GB for q = 4 (64 sub-cells); unit tests 13
  passed (the full-window exactness also for q = 4).

**Queue order (operation, 2026-10-20 about 01:20 local).** The corrected arms need about 3 times longer than the voxel
arms (relaxation 570 … 1 620 s against 180 s; cycle about 190 s against 53 s, with shared GPU). To have the
references of D1 by the morning, the queue was restarted: the three running corrected arms (D3Pc_s941, D3Lc_s941,
D3Pc_s942) go on; a new dispatcher (spool5, 3 slots) runs D15 and D15g first, then the other D3c arms, D3f and D3c4.
GPU load before: 83 % with 3 jobs.

### 3.8 Stage A with q = 4 (centroid term, one image shell): the same result; the staircase effect grows with q (2026-10-20, night mandate)

**Run** (commit 41df994, files results/static_A_*_q4_k1_c.json). Fine meshes 608³ (K2) and 624³ (K4), one image shell
(D5). R0 passed: K2 144 896 modified and 72 032 new cells, φ_f 0.62001; K4 142 208 and 69 600, φ_f 0.55003. R1
passed (largest deviation from φ (1 − φ)/3: 1.2·10⁻⁵).

| geometry | S (K_d) | S / \|E_K\| | E_frac − E_ref | P3: e_E / e_h | L3: e_E / e_h | P2: e_E / e_h | L2: e_E / e_h |
|---|---|---|---|---|---|---|---|
| K2_0 (K2_1 within 0.001) | 7.96e-4 | 0.076 | +3.8e-5 | 0.014 / 0.004 | 0.011 / 0.057 | 0.039 / 0.011 | 0.036 / 0.059 |
| K4_1 | 7.10e-4 | 0.385 | +3.1e-5 | 0.015 / 0.004 | 0.013 / 0.055 | 0.040 / 0.011 | 0.039 / 0.057 |

- R2: both forms "sufficient" from p = 2 on, as with q = 2. The accuracy of the correction does not depend on q
  (from 2 to 4).
- The staircase effect against the finer target surface is larger: S(q = 4) / S(q = 2) = 1.69 (K2) and 1.69 (K4).
  If S approaches the smooth sphere as 1/q (an assumption; two points only), the limit is about 1.1·10⁻³ K_d for K2
  (0.11 |E_K|) and 1.0·10⁻³ K_d for K4 (0.54 |E_K|) (estimates).
- The volume fraction alone changes its sign of error between q = 2 (−0.45 S) and q = 4 (+0.05 S of the larger S);
  its field error in the layer stays about 0.3 of that of the voxel model (P3: 0.004).

### 3.9 Stage B on the geometry of the dynamic arms (d/l_ex 96, φ 0.55): the correction holds; the staircase effect is about 1 … 2 |E_K| (2026-10-20, night mandate)

**Decision under the night mandate.** Stage B of PLAN.md was planned on the B0b geometry (φ 0.65); after D2 the
dynamic arms use φ 0.55. Stage B therefore uses the geometry and the states of the dynamic arms: the prepared states
(init.pt) of the four finished D3 runs (d/l_ex 96, φ 0.55, L_eff 12 l_ex, r_p 0, dx 3, Newell solver with full PBC,
option C, seeds 941 … 944). Production setting on both meshes (fine meshes 96³ and 192³), centroid term, p = 1 … 3.
Files: results/stageB_static_A_d96_phi0.55_q2_c.json and _q4_c.json.

- Step 2 passed for all four states (total energy against the run record within 4·10⁻⁹ K_d; |m| − 1 at most
  1.5·10⁻⁷; 0 Bloch points). Thus init.pt holds the prepared state of the option-C run.
- R0 passed (q = 2: 8 736 modified, 4 704 new cells; q = 4: 13 760 and 7 200; φ_vox 0.54948, φ_f 0.55241 / 0.54971).
  R1 passed (largest deviation 1.3·10⁻⁸).

| q | S (K_d), 4 states | S / \|E_K\| | E_frac − E_ref | P3: e_E / e_h | L3: e_E / e_h | P2: e_E / e_h | L2: e_E / e_h |
|---|---|---|---|---|---|---|---|
| 2 | 1.398 … 1.420 e-3 | 1.03 … 1.10 | −5.8e-4 | 0.012 … 0.013 / 0.002 … 0.003 | 0.013 … 0.015 / 0.053 … 0.055 | 0.036 … 0.037 / 0.007 | 0.036 … 0.039 / 0.054 … 0.056 |
| 4 | 2.569 … 2.611 e-3 | 1.89 … 2.02 | +3.1e-4 | 0.010 / 0.002 … 0.003 | 0.010 … 0.011 / 0.055 … 0.056 | 0.027 … 0.028 / 0.007 … 0.008 | 0.027 … 0.029 / 0.056 … 0.058 |

- R2: both forms "sufficient" from p = 2 on, in all four states and for both q.
- At this geometry (dx/d 0.033, three times coarser relative to the particle than the production geometry) the
  staircase raises the stray-field energy by about one (q = 2) and two (q = 4) times the anisotropy energy of the
  state. S(q = 4) / S(q = 2) = 1.84. A strong effect of the correction on the loss is therefore possible in the
  dynamic arms.

**D6 — second GPU box for the reference arm D15g (decision under the night mandate, 2026-10-20 about 02:20 local).**
- Finding: with 6 jobs the GPU of box 1 is at 100 %; the corrected arms need more than 1.5 h for the AC pass alone,
  and D15_s943 needed 56 min for the relaxation. The D1 reading needs D15 and D15g in all 4 seeds; on one GPU they
  would not finish by the morning.
- Options: (a) one GPU, accept a partial result ("not determined" with n < 3 is likely); (b) a second cheap box for
  D15g; (c) drop D15g and take the resolution effect of m from V3g (7.110; other φ and protocol, not clean).
- Recommendation and decision: (b). Offer 44133018 (RTX 3090, 0.163 $/h, US), instance 56811357, about 6 h ≈ 1 $
  (estimate). The queue of box 1 was restarted without D15g (spool6: D15_s944, the other D3c arms, D3f, D3c4); the
  running jobs go on. Credit before the rent: 14.82 $.

### 3.10 R3 (sign of the corrected energy): method changed (D7); running (2026-10-20, night mandate)

- Finding: the first R3 runs (scipy eigsh / ARPACK, smallest algebraic eigenvalue, tol 10⁻³, ncv 40, maxiter 400)
  gave no result after 2.4 h of CPU (box 1, K2 and K4, P and L at p = 2 and 3). Cause (inference): maxiter counts the
  restarts of ARPACK, not the field evaluations; the operator has a large cluster of eigenvalues near 0 (divergence-free
  patterns), so the smallest eigenvalue converges slowly.
- Options: (a) wait (unknown time, possibly many hours); (b) a plain Lanczos run with a fixed number of steps (300
  field evaluations) and progress output: the smallest Ritz value is an upper bound of λ_min, so a value below
  −10⁻⁴ λ_max proves an unusable form, and a converged positive value (residual ≤ 10⁻⁴ λ_max) passes; (c) only the
  dynamic stability of the corrected runs as evidence.
- Recommendation and decision: (b) for the forms of interest (P3 and L3, K2 and K4, q = 2, centroid term), on the CPU
  of box 2, and (c) in addition (energy balance and closure of the corrected dynamic runs).
- First values after 25 steps (P3): smallest Ritz value +1.20·10⁻³ (K2) and +1.19·10⁻³ (K4), largest 0.997. About
  29 s per step on 3 threads; 300 steps take about 2.4 h per case.

**R3 results (plain Lanczos, 300 steps, centroid term, q = 2, p = 3; files results/sign_*.json).** λ_Ritz,min is
the smallest Ritz value (an upper bound of λ_min), r its residual; λ_Ritz,max = 0.9997 in all cases.

| geometry | operator | λ_Ritz,min | r | reading of the rule |
|---|---|---|---|---|
| K2 | P3 (pairwise) | +1.0143e-4 | 1.50e-4 | not determined |
| K2 | reference without window term (volume fraction, positive semi-definite) | +1.0120e-4 | 1.50e-4 | not determined |
| K2 | L3 (cell-local) | −2.192e-3 | 2.0e-5 | **not usable** |
| K4 | P3 | +1.0425e-4 | 1.50e-4 | not determined |
| K4 | reference without window term | +1.0404e-4 | 1.50e-4 | not determined |
| K4 | L3 | −2.322e-3 | 3.2e-6 | **not usable** |

- The cell-local form has a converged negative eigenvalue of about −2.3·10⁻³ (in units of the largest one) on both
  geometries. Thus the stray-field energy of L is indefinite: R3 "not usable". Probable cause (inference): replacing
  m_j by m_i makes the energy of patterns that alternate from cell to cell negative at the surface. In the dynamics
  the exchange penalises such patterns strongly (D3Lc_s941 ran stable), but the stray-field operator itself has an
  unphysical negative mode.
- The pairwise form gives the same smallest Ritz value as the positive semi-definite reference within 2·10⁻⁷ after
  300 steps (both started from the same random vector). The residual stays at 1.5·10⁻⁴ for both, because the
  operator has a cluster of eigenvalues near 0; the strict criterion of R3 (residual at most 10⁻⁴) cannot be reached
  with 300 steps even for the reference.

**D8 — reading of stage A and the dynamic arms (decision under the night mandate).**
- Options for R3 of P: (a) keep "not determined" and give no outcome; (b) take the agreement with the positive
  semi-definite reference within 2·10⁻⁷ as the evidence that P adds no negative mode, plus the stable dynamic runs
  D3Pc_s941 and D3Pc_s942 (closure 2·10⁻⁸ and 8·10⁻⁷, w_dis / w_loop 1.000); (c) more Lanczos steps.
- Recommendation and decision: (b), with (c) as an option for Chris. **Outcome of stage A with the centroid term:
  "go on, pairwise"** (P sufficient from p = 2 on in R2 for both geometries, both q and the dynamic-test geometry;
  R3 as in (b)). L is "not usable" (R3).
- Dynamic arms: the remaining arms with L (D3Lc_s943 started and stopped, D3Lc_s944, D3Lc4) are cancelled; D3Lc_s941
  and D3Lc_s942 (finished) stay as information. The queue (spool7) runs D3Pc_s944, D3f and D3Pc4.

### 3.11 End of the night mandate; option B of the dynamic test (2026-10-20 morning, decisions of Chris)

- About 06:40 local: Chris is back; all decisions go through him again. The night mandate (decisions D1 … D8) ends.
- Question of Chris: how do we test that the correction works? Answer: statically against the exact fine reference
  (3.6 … 3.10), dynamically by D1. Weakness of D1: each arm prepares its own state, so the scatter between the states
  is large (SD of ln P 0.16 for D3, 0.54 for D3Pc with n 2). Options: A as is; B transferred start (the prepared state
  of D15, block averaged, as the start of the coarse arms, as the T arm of V3g in 7.110); C more seeds. **Decision of
  Chris: B** (A runs on as information).
- B: arms T3 (dx 3, voxel), T3Pc (dx 3, pairwise correction with the centroid term, q 2, p 3) and T15g (dx 1.5, dx 3
  geometry locked), seeds 941 … 944, each from the prepared state of D15_s<seed> (init.pt, Newell, option C, no Bloch
  points), then the production relaxation and the 9 mT stage. The reference is D15 itself (same prepared state).
  Reading as D1 with Δ_c = ln P(T3Pc) − ln P(T3) and Δ_g = ln P(D15) − ln P(T15g).
- Code (run_loops --init_transfer): allowed with the fraction surface of --staircase; the keys of the surface model,
  of the locked geometry and of the preparation may differ from the source; the geometry of the source comes from its
  own config. Smoke run on box 2 (transfer from D15g_s941 to the corrected dx 3 model): transfer and relaxation start
  without error. Study tests 118 passed, 3 skipped.
- Question of Chris: is the reference the smooth sphere? No. With q = 2 the target and the reference are the fine
  staircase of dx 1.5, which still has steps. This is a validation step. The smooth sphere is the limit q → ∞ of the
  correction at dx 3 (same cost per time step for any q). Proposal (open for Chris): after the validation, S(q) for
  q = 2, 4, 8 statically on the d/l_ex 96 geometry, then P(q) dynamically at dx 3 and an extrapolation to q → ∞.
  Limits: the exchange still sees the coarse staircase; the resolution of m stays at dx 3.
- 2026-10-20 (decision of Chris): the q = 4 arms D3Pc4 are cancelled before their start. They have no dynamic
  reference (dx 0.75 runs) and belong to the series in q, which is not started. The static result for q = 4 stays
  (3.8, 3.9). Box 1 queue: spool8 (D3f_s944, then T3 and T3Pc).
- 2026-10-20 (decision of Chris): focus on q = 2 as the proof of concept (option A); no dynamic validation of higher
  q for now. Options given: B smaller particles (d/l_ex 48, reference dx 0.75), C q = 4 at d/l_ex 96 (reference
  dx 0.75, 192³), D coarse mesh dx 6 with q = 4 against the existing dx 1.5 references; D would be the first step
  if higher q is wanted later.

### 3.12 Why the method of Donahue alone is not sufficient, what the extension does, and the evidence; decision of Chris: the extension is part of the method (2026-10-20)

**The method of Donahue in one picture.** The sphere is built from coarse cubes. A cube at the surface holds only a
part of the material (a cut cell). Near pairs (up to p cells apart) use the exact coupling of the true shapes, from
fine sub-cubes. Far pairs use the coarse coupling: the cut cell acts as a smaller magnet (moment v M_s Δx³, v the
volume fraction) at the centre of its cube.

**Why this is not sufficient here.** The material of a cut cell sits at one side of the cube; its centroid lies up to
about 0.43 cells from the centre (one body sub-cell at q = 2; own calculation). For one far pair the error is very small (it falls as 1/r⁴ with the distance r, a
misplaced dipole). But a sphere has very many cells in its interior, and each of them feels the misplaced surface
magnets a little. The sum of these small errors then falls only as 1/R with the window radius R (inference, 3.6; the
test with the extension agrees, 3.7). Measured (pairwise form P, K2_0, energy error in % of the staircase effect S,
3.6 and 3.7):

| p | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| method of Donahue (P) | 95 | 59 | 43 | 34 | 28 | 24 |
| with the extension (P) | 19 | 5 | 2 | – | – | – |

Without the extension the error follows 1/(p + 0.6); an error of at most 10 % would need p of about 14 (estimate from
the law; a window of 29³, about 24 000 neighbours per surface cell, own calculation), which is not practical. A probable reason why Donahue & McMichael did not see this
(inference): they tested only films one cell thick, which have almost no interior.

**The extension.** In the far field each cut cell carries its moment at the true centroid of its body part, not at
the centre of its cube (the first moment). Technically this is one additional term with the central difference of
the same periodic kernel (factor i sin ω in Fourier space): no new kernel, the full PBC and the FFT stay; the window
term subtracts the same far-field coupling (3.7).

**Evidence that the extension is necessary.**
- Without it the error falls only as 1/(p + 0.6) (p = 1 … 4 on 2 geometries, p = 5 and 6 on K2); the energy
  criterion of R2 cannot be met with a practical window.
- The cause was stated as a hypothesis before the test (3.6); the test with the extension then showed the expected
  change: the error now falls by a factor of about 4 from p = 1 to 2 and about 2.7 from p = 2 to 3; before, each step
  left 0.63 and 0.73 of the error.
- The window term stays necessary as well: with the extension but without the window term the error is −28 % of S
  (K4_1, q = 2).

**Evidence that it is sufficient (static).**
- R2 holds for the pairwise form from p = 2 on in all 7 states, on 3 geometries (K2: d/l_ex 300, φ 0.62; K4:
  d/l_ex 300, φ 0.55; the dynamic test geometry: d/l_ex 96, φ 0.55), for q = 2 and q = 4 (3.7 … 3.9). At p = 3 the
  energy error is 1 … 2 % of S and the field error in the surface layer 0.2 … 0.4 % of that of the voxel model.
- Controls: with a window over the whole box the result is exact (deviation below 10⁻¹⁰, unit test), so the
  extension does not change the exact limit; the operator stays symmetric (consistent energy); the pairwise form
  has no negative Ritz value and the same smallest Ritz value as the positive semi-definite reference; by the strict
  rule R3 stays "not determined" (3.10, D8). The cell-local form L is "not usable" (R3, negative mode).
- Not yet shown: the sufficiency for the dynamics (does the loss follow the finer staircase?). This is the running
  validation against dx 1.5 (D1 and option B). Only the stray field is corrected, not the exchange. The extension is
  not in the literature found (3.1, 3.2); it is a construction of this work.

**Decision of Chris (2026-10-20).** The centroid term (first moment of the cut cells in the far field) is part of the
method. The method of this work is: V1 cut cells (Donahue & McMichael 2007) with the pairwise window term and the
centroid term in the far field. The decision concerns the pairwise form; the cell-local form is not used.

**Correction to 3.12 (2026-10-20, remark of Chris).** The sentence "Only the stray field is corrected, not the
exchange" is not exact. With the fraction surface of V1 the surface cells (cut and new cells) carry A v (and M_s v,
K v); the exchange between cells uses the harmonic mean of A (7.29: the exchange field of a partial cell is about as
strong as inside; the model is energy-consistent). Thus the exchange of V1 is no longer that of the staircase model,
but it is also not the exchange of a smooth surface (the coupling still runs over whole cube faces; there is no
boundary condition on the true surface as in García-Cervera et al. or Parker et al.). Its error against the smooth
sphere is not measured. Decision of Chris: a static exchange test (3.13).

### 3.13 Static exchange test: V1 keeps a first-order surface error of the exchange, smaller than that of the voxel model (2026-10-20, decision of Chris)

**Method** (check_exchange.py, box 1 CPU, double precision; file results/exchange_static.json). Exchange of magnum.np
(harmonic mean of A between cells, PBC), voxel model (A, M_s on the voxel cells) against V1 (A v and M_s v from q³
sub-cells, q = 2 and 4), on K2, K4 and the dynamic test geometry (d/l_ex 96, φ 0.55), dx 3.
- Helix m = (cos k·r, sin k·r, 0), wave vectors commensurate with the box (n = (3,0,0), (2,2,0), (2,2,2), (6,0,0),
  (4,4,0), (3,3,3) times 2π/a; wavelength about d/2 and d/4). Reference: the bulk density of the same stencil times
  the volume of the smooth spheres; thus the stencil error of the bulk drops out and the surface error remains.
  Control: a box without surfaces gives the reference within 3·10⁻¹⁴.
- Vortex around an axis through each particle centre (core profile linear up to r_c = 0.1 R and 0.3 R), against the
  continuum integral (information: stencil and surface error together).
- A first run had wave vectors not commensurate with the box; the phase jump at the box edge gave spurious errors up
  to +339 %. It is not used.

**Result (helix; the error is the same for all 6 wave vectors of a geometry).**

| geometry | voxel model | V1, q = 2 | V1, q = 4 |
|---|---|---|---|
| K2 (d/l_ex 300, φ 0.62) | −1.39 % | −1.06 % | −0.96 % |
| K4 (d/l_ex 300, φ 0.55) | −1.43 % | −1.02 % | −0.97 % |
| d/l_ex 96, φ 0.55 | −4.91 % | −3.01 % | −3.26 % |

- For a helix each bond between two cells carries the same energy, so the test counts the missing share of bonds at
  the surface: the first-order surface error of the exchange (it grows with Δx/d). The volume weighting of V1 lowers
  it by about 30 … 40 % but does not remove it. The exchange of V1 is neither that of the staircase model nor that of
  a smooth surface (remark of Chris, 3.12 correction).
- Vortex (r_c 0.1 R): K2 −2.4 % (voxel and V1), K4 −2.2 … −2.4 %, d/l_ex 96 −8.3 % (voxel), −6.5 % (V1 q 2).
- Size against the stray field (estimate): the exchange energies of the states are 1.4e-3 (K2), 3.8e-4 (K4) and 2.4e-3
  K_d (d/l_ex 96); a surface error of 1 % (K2, K4) and 3 % (d/l_ex 96) is then about 1 … 5 % of the staircase effect S
  of the stray field (estimate; the gradients of real states sit in walls and vortices, not evenly). This is the same
  order as the rest of the stray-field error after the correction (1 … 2 % of S at p = 3).
- Consequence (proposal): after the stray-field correction the exchange is not the dominant remaining error, but not
  negligible for small particles. A surface correction of the exchange (Neumann condition on the true surface after
  García-Cervera et al. 2003) is the next step if the dynamic validation shows a deviation.

### 3.14 Dynamic test, design A complete: the reading criterion is met, but the test has no power; design B gives the decision (2026-10-20)

**Data** (results/dynamic_A.md and .json; eval_dynamic.py --design A; box 1, D15g on box 2). Loss density P at
9 mT, d/l_ex 96, φ 0.55, seeds 941 … 944. Each arm prepares its own state (option C, Newell solver with full PBC).
All prepared states have no Bloch points. All 4 seeds are complete for D3, D3f, D3Pc, D15 and D15g.

| arm | P (W/cm³), seeds 941 / 942 / 943 / 944 | SD of ln P |
|---|---|---|
| D3 (dx 3, voxel model) | 1.02 / 0.799 / 0.697 / 0.775 | 0.16 |
| D3Pc (dx 3, pairwise correction with centroid term, q 2, p 3) | 1.27 / 0.594 / 3.03 / 0.585 | 0.78 |
| D3f (dx 3, volume fraction alone, q 2) | 0.882 / 12.1 / 0.64 / 0.717 | 1.40 |
| D15 (dx 1.5, voxel model of dx 1.5) | 1.88 / 1.39 / 0.623 / 1.16 | 0.47 |
| D15g (dx 1.5, voxel geometry of dx 3 locked) | 1.08 / 0.871 / 1.10 / 0.875 | 0.13 |

D3Lc: seeds 941 and 942 only (the other seeds cancelled by D8, CANCELLED_D3Lc.txt; the label "refused or crashed" of
the evaluation for D3Lc_s943 is wrong, the run was cancelled). D3P: seed 941 only. Both are information.

**Reading of D1** (fixed before the results, 3.5):

| quantity | per seed 941 / 942 / 943 / 944 | mean ± SE (n 4) |
|---|---|---|
| Δ_g = ln P(D15) − ln P(D15g) | +0.55 / +0.47 / −0.57 / +0.28 | +0.18 ± 0.26 |
| Δ_c = ln P(D3Pc) − ln P(D3) | +0.22 / −0.30 / +1.47 / −0.28 | +0.28 ± 0.42 |
| d = Δ_c − Δ_g | −0.33 / −0.77 / +2.04 / −0.56 | +0.10 ± 0.65 |

- Formal reading: "correction reproduces the geometry effect" (the mean of d is below 2 SE and below 0.15).
- This reading has no value. The SE of d (0.65) is 4 times the bound 0.15. The 95 % interval of the mean of d
  (t distribution, 3 degrees of freedom) is −1.99 … +2.17, so the test cannot tell a correct correction from a wrong
  one. Design A does not decide; design B (3.11) is the test.
- Cause of the scatter: the own preparation of each arm ends in different metastable states. Seed 943 shows it in
  both directions: D3Pc_s943 has a loop energy w 4.4 times that of D3_s943 (1.13·10⁻⁷ against 2.59·10⁻⁸, in units
  of the production log), and D15_s943 has a low P (0.62 W/cm³) against D15g_s943 (1.10 W/cm³). Without seed 943 the
  values of d are all negative (−0.33 … −0.77); a seed is not left out after the result, so this is only a remark.
- Gap of the criterion (finding): the criterion of D1 has no condition on the power of the test. A large SE makes the
  first condition (mean at most 2 SE) easy to meet. Options given: A keep the criterion; B add the condition
  "2 SE at most 0.15" for the reading "reproduces", else "not determined"; C as B plus 4 more seeds (945 … 948, about
  10 … 11 h, about 1.7 $ on box 2, estimate). **Decision of Chris (2026-10-20): B.** The decision on C comes after
  the result of design B with n 4. Note: B was fixed after the interim result of design B with n 2 was known; the
  added condition acts only on the SE, not on the mean. With B the reading of design A is "not determined".
  eval_dynamic.py prints also the 95 % interval of the mean of d (t distribution).
- Resolution of m at the dx 3 geometry, ln P(D15g) − ln P(D3): +0.18 ± 0.09 (SE, n 4). The correction does not act
  on this part (the resolution of m stays at dx 3).

**Design B, state at about 08:00 UTC** (information, n 2 for T3Pc; seeds 943 and 944 of T3Pc run on box 2, end about
09:45 … 10:00 UTC): Δ_g = ln P(D15) − ln P(T15g) = +0.32 ± 0.11 (SE, n 4); resolution of m, ln P(T15g) − ln P(T3) =
+0.08 ± 0.02 (n 4); Δ_c = ln P(T3Pc) − ln P(T3) = +0.42 ± 0.11 (n 2); d = −0.07 ± 0.11 (n 2; per seed +0.05, −0.19).
No reading before n 4. The T arms start from the prepared states of D15 (no Bloch points); the label "settled ?" of the
evaluation only says that the T runs have no preparation of their own. The transferred start removes most of the
scatter: SD of ln P 0.22 (T3) and 0.26 (T15g) against 0.16 … 1.40 in design A, and T3 agrees with D3 within about 0.01 in
ln P for seeds 941, 942 and 944.

### 3.15 Dynamic test, design B complete (n 4): reading "not determined"; the large value of seed 944 is not a change of state (2026-10-20)

**Data** (results/dynamic_B.md and .json; eval_dynamic.py --design B; T arms and T15g on box 2, T3 and D15 on box 1).
Loss density P at 9 mT, d/l_ex 96, φ 0.55. All T arms start from the prepared state of D15 of the same seed (no Bloch
points). All runs have at least 2 calm cycles; cycles 5 … 7 of T3Pc_s944 give the same loop energy within 10⁻³.

| arm | P (W/cm³), seeds 941 / 942 / 943 / 944 | SD of ln P |
|---|---|---|
| T3 (dx 3, voxel model) | 1.03 / 0.794 / 0.605 / 0.776 | 0.22 |
| T3Pc (dx 3, pairwise correction with centroid term, q 2, p 3) | 1.75 / 1.07 / 0.611 / 1.56 | 0.47 |
| D15 (dx 1.5, voxel model of dx 1.5; reference) | 1.88 / 1.39 / 0.623 / 1.16 | 0.47 |
| T15g (dx 1.5, voxel geometry of dx 3 locked) | 1.16 / 0.857 / 0.615 / 0.874 | 0.26 |

| quantity | per seed 941 / 942 / 943 / 944 | mean ± SE (n 4) |
|---|---|---|
| Δ_g = ln P(D15) − ln P(T15g) | +0.48 / +0.48 / +0.01 / +0.28 | +0.32 ± 0.11 |
| Δ_c = ln P(T3Pc) − ln P(T3) | +0.53 / +0.30 / +0.01 / +0.70 | +0.38 ± 0.15 |
| d = Δ_c − Δ_g | +0.05 / −0.18 / −0.00 / +0.42 | +0.07 ± 0.13 |
| resolution of m, ln P(T15g) − ln P(T3) | | +0.08 ± 0.02 |

- **Reading (criterion K2, 3.14): "not determined".** The mean of d is small (+0.07), but 2 SE (0.25) is above the
  bound 0.15. The 95 % interval of the mean of d (t distribution, 3 degrees of freedom) is −0.33 … +0.47.
- In 3 of 4 seeds the correction follows the geometry effect of the same seed within 0.18 in ln P, also where the
  geometry effect is near zero (seed 943). Seed 944 gives most of the scatter: Δ_c +0.70 against Δ_g +0.28; P of
  T3Pc_s944 (1.56 W/cm³) is above that of the reference D15_s944 (1.16 W/cm³).
- Reading of the interim result with n 3 (K2 met, d −0.05 ± 0.07) is replaced by this result; n 4 was the agreed
  sample.

**Check of seed 944: did a T arm leave its start state?** (check_state_overlap.py on box 2; results/state_overlap_B.md;
go of Chris 2026-10-20). Reference of each arm: the prepared state of D15, block averaged onto the mesh and geometry of
the arm (the start of the arm). Overlap: mean of m · m_ref over the magnetic cells.
- After the relaxation and at the end of the cycles: T3 0.987 … 0.988, T3Pc 0.9994 … 0.9997, T15g 0.978 … 0.979, for
  all 4 seeds. Share of cells with m · m_ref < 0: at most 0.04 %. T3Pc against T3 at the end of the cycles: 0.988 for
  all seeds. The mean m of each particle is below 0.01 (demagnetised states).
- Result: the hypothesis "T3Pc_s944 went into another metastable state" is rejected. All T arms stay in their start
  state; seed 944 is not different from the other seeds in this measure. The large Δ_c of seed 944 is a property of
  the corrected model in the same state. The check does not show what happens inside a cycle; the state at the end of
  a cycle is the same in all arms. (Corrected in 3.16: the first version of this item said that the loss comes from
  irreversible jumps; at this drive there are no jumps.)
- T15g has the lowest overlap because 8736 of its 486144 cells (1.8 %) are outside the fine staircase and take a
  replacement value at the start (3.11).
- First run of the check: the script used the nominal diameter of the config (96 l_ex) instead of the diameter of
  box_edge (92.2 l_ex); the contact check of the geometry stopped it. Corrected; the cell counts now agree with the
  runs (60768, 65472, 486144).

**Cost state** (09:05 UTC): credit 7.07 $. Box 2 (0.18 $/h) has no run since 09:02 UTC. The instance 53995473 (other
session, 0.41 $/h) also uses the credit; at 0.59 $/h together the credit lasts about 12 h.

**Decision of Chris (2026-10-20, about 09:30 UTC): option 1 (analysis of seed 944, local) and option 2 (more seeds)
together.** Options given: 1 the M(H) loops of the four arms of seed 944 (where in the cycle and in which particle the
extra loss forms; local, no cost); 2 more seeds until 2 SE is at most 0.15 (with SD of d about 0.25 about 11 … 12
seeds, thus 8 more; estimate); 3 stop with the reading "not determined". Recommendation was 1 first.
- Started on box 2 (about 09:35 UTC): seeds 945 … 952, arms D15, T3, T3Pc, T15g (jobs/50_D15a, 50_D15b, 60_Ta,
  70_Tb). The job lines are those of seeds 941 … 944 with only the seed changed, and --init_wait_h 24 instead of 8 (the
  T arms wait longer for the prepared state of D15; no effect on the result). Queue: spool (4 slots) D15 945 … 948,
  then the T arms; spool2 (2 slots) D15 949 … 952. Duration about 20 h, box cost about 3.6 $ (estimates).
- Credit: at 0.59 $/h (box 2 and the instance 53995473 of another session) the credit of 7.07 $ lasts until about
  21 h UTC. Chris was told that a top-up is necessary before then.

### 3.16 Seed 944 and the other seeds, particle by particle: the loss is a small phase lag; the correction acts in the right particles, but its size per particle is wrong by up to ±0.6 in ln w (2026-10-20, option 1 of Chris)

**Method** (check_loops.py, local; results/loops_B.md, .json, .png). Last cycle of each arm of design B, seeds
941 … 944. From the samples (256 per cycle): the loop energy of each particle (∮ H dM_par,k; its share times the loop
energy w of the box), and the fundamental of M_par,k(t) (amplitude, phase lag behind H, the rest).

**What the loss is at this drive** (f 29.9 MHz, H0 1340 … 1396 A/m):
- M_par,k(t) of each particle is a sinusoid: the rest without the fundamental is at most 0.6 % of the amplitude, in
  all arms and seeds. There are no jumps in the cycle. The loss is the out-of-phase part of a linear response: the
  phase lag of the fundamental is 0.3 … 3.6 mrad. The amplitudes of the four arms are the same (the arm with the
  fraction surface shows 7 % less because its particle means count each cut cell fully; only the shares are used).
- Check of this reading: for a linear response the loss of particle k is (φ/4) π μ0 H0 A_k sin δ_k (A_k amplitude,
  δ_k phase lag). For the 48 particle values of the voxel arms (T3, D15, T15g) this gives the measured w_k within
  −0.3 … +0.9 % (ratio 1.0007 in the mean).
- Thus a factor 1.5 in the loss of a particle is a change of its phase lag from, e.g., 2.6 to 3.6 mrad (seed 944,
  particle 3: T3 1.22, T15g 1.49, D15 2.62, T3Pc 3.57 mrad).

**Particle by particle** (Δ_c,k = ln w_k(T3Pc) − ln w_k(T3), Δ_g,k = ln w_k(D15) − ln w_k(T15g); table in loops_B.md):
- Seed 944: the correction acts in the same particles and in the same direction as the fine staircase, but stronger:
  particle 3 +1.13 against +0.60, particle 1 +0.68 against +0.28, particle 2 −0.70 against −0.46; particle 0 +0.43
  against −0.07. In D15 and T3Pc particle 3 carries 53 … 54 % of the loss, in T3 and T15g 35 … 39 %.
- Seed 941: the box values agree (d +0.05, 3.15) only because errors of opposite sign cancel: particle 0 +0.93
  against +0.35, particle 1 +0.02 against +0.65.
- Seed 942: the correction is a little too weak in all 4 particles (d_k −0.26 … −0.01). Seed 943: all effects small,
  d_k within ±0.17.
- The sign of Δ_c,k agrees with that of Δ_g,k in 13 of 16 particle-seed pairs; the 3 other pairs have |Δ_g,k| at most
  0.14. Over the 16 pairs: correlation 0.68; RMS of d_k 0.33 against 0.41 without a correction. The pairs are not
  independent (the particles of one box couple through the stray field), so these numbers are a description, not a
  test.

**Result.** Seed 944 is not an exception at the level of the particles: seed 941 has particle errors of the same size
(±0.6), which cancel in the box value. The correction moves the loss into the right particles, but the size of the
change in one particle can be wrong by a factor up to 1.8 (exp 0.6). This agrees with the SD of d of about 0.25 at
box level (3.15) and supports the estimate of 8 more seeds for K2. Cause of the particle errors: open. The loss is a
phase lag of a few mrad, so small changes of the local dynamics at the surface change it much. The resolution of m
alone is too small as a cause: per particle, ln w_k(T15g) − ln w_k(T3) (same geometry, dx 1.5 against dx 3) is
−0.14 … +0.22 (RMS 0.11), against particle errors up to ±0.6. Candidate not checked: the exchange of the cut cells
(3.13).
- Consequence (for Chris): even if more seeds give a mean d near 0, one corrected run at dx 3 differs from its run at
  dx 1.5 by about ±0.25 in ln P (SD of d), because the particle errors do not cancel in a box of 4 particles. Whether
  "no bias in the mean, ±0.25 per realisation" is sufficient for the main study is a model question for Chris.

### 3.17 Search for the cause, step 1: energy of each term and particle in the start states; the static mismatch is the same in all particles (2026-10-20, go of Chris)

**Method** (energy_budget.py, local CPU, double precision; results/energy_budget_B.md and .json). Each start state
(init.pt) of the four arms of design B, seeds 941 … 944, in the model of its own arm, built as run_loops builds it.
Energy of each linear term in each particle: sum of −½ μ0 M_s,cell m · h V_cell over the cells of the particle, in
units of K_d V_box. Check: the sum of the terms gives the energy of the run record within 1.5·10⁻⁶ (relative) for all
16 states (for D15 the record after the AC preparation). Differences as in 3.16: δ_c = E(T3Pc) − E(T3),
δ_g = E(D15) − E(T15g).

**Result** (10⁻⁶ K_d V_box per particle; range over the 16 particle-seed pairs):

| term | δ_c (correction at dx 3) | δ_g (fine staircase at dx 1.5) | δ_c − δ_g |
|---|---|---|---|
| stray field (with the correction term) | −384 … −377 | −245 … −237 | −142 … −138 |
| exchange | −53 … −47 | −109 … −103 | +55 … +57 |
| cube anisotropy | −9 … −2 | −10 … −2 | −1 … +2 |

- Sum of the three terms: δ_c −441 … −432 (mean −436), δ_g −357 … −349 (mean −352); the total geometry effect at
  dx 3 is 24 % larger than at dx 1.5, the stray field alone 58 %, because the exchange mismatch partly cancels.
- The mismatch of the stray field and of the exchange is the same in all 16 particle-seed pairs within ±2 %. This is
  forced by the geometry: the four sphere centres sit on cell vertices, so the four particles have the same staircase
  (checked: 15192 voxel cells at dx 3, 16368 V1 cells, 122184 cells at dx 1.5 for each particle). A term that is
  wrong in one particle only cannot exist. The loss errors of 3.16 change sign from particle to particle (d_k
  −0.63 … +0.58); they can only come from how the different magnetic state of each particle reacts to the same
  staircase in the dynamics. Correlation of the mismatch with d_k over the 16 pairs: stray field
  +0.44, exchange +0.09, anisotropy −0.49 (the pairs are not independent; no test).
- Meaning of the uniform mismatch: for the same coarse m the correction gives the stray field of the fine staircase
  (static test, 1 … 2 % of S, 3.9). The stray-field effect of the staircase is −1.52·10⁻³ K_d V_box (4 particles) at
  the resolution of m of dx 3 and −0.96·10⁻³ at dx 1.5. Interpretation (not tested): with the finer m the cells at
  the steps of the coarse staircase turn and lower its surface charges, so the coarse staircase costs less. The
  effect of the geometry thus depends on the resolution of m. Support (information): for the same block-averaged m
  before the relaxation (coarse pair only) δ_c of the stray field is −521 per particle; after the relaxation in each
  model −380; at dx 1.5 −240. More freedom of m gives a smaller stray-field effect of the staircase.
- Consequence for the dynamic test: d = Δ_c − Δ_g contains two parts that the four arms cannot separate: the error
  of the correction (the corrected model against the exact stray field of the fine staircase, both with the m of
  dx 3) and the interaction of geometry and resolution of m (the geometry effect at dx 3 against that at dx 1.5). The
  static accuracy of the correction suggests that the second part is the larger one; this is not shown.
- A same-state comparison before the relaxation was tried and dropped: the transferred start of T15g has 8736 cells
  (1.8 %) with the m of void cells of D15, which gives an exchange energy far from the relaxed one.

### 3.18 Reference arm X3: coarse m with the exact stray field of the fine staircase (2026-10-20, decision of Chris: option B)

**Decision of Chris (about 10:20 UTC): B.** The 8 more seeds run on; in addition the reference arm X3 for the seeds
941 … 944. Options given: A only the seeds; B seeds and X3; C X3 only.

**Purpose.** d = Δ_c − Δ_g (3.15) mixes the error of the correction and the interaction of geometry and resolution of
m (3.17). X3 separates them: ln P(T3Pc) − ln P(X3) is the error of the correction alone (window p 3 and centroid term
against the exact field, both with the m of dx 3); ln P(X3) − ln P(D15) is the interaction (plus the effect of the V1
weights of exchange and anisotropy against the fine staircase, which X3 cannot separate).

**Model X3** (magnumnp DemagStaircaseFine; run_loops --staircase fine). The V1 cells and their weights v M_s, v A,
v K as in T3Pc. Stray field: each body sub-cell of the fine mesh (q = 2, dx 1.5) takes m of its cell; one FFT on the
fine mesh with the Newell solver with full PBC (extrapolation, 4 shells, as production); the field of a cell is the
mean of the fine field over its body sub-cells. This is the gradient of the stray-field energy of the fine body
(variational) and the same operator as the static reference of check_static. No window, no centroid term.

**Tests.**
- test_staircase.py (15 passed): open boundaries: X3 equals DemagField with v M_s plus the pairwise correction with
  the full window within 10⁻¹⁰ M_s; with PBC: X3 equals the static reference within 10⁻¹⁰ M_s, its energy equals the
  energy of the fine body within 10⁻⁹, and the operator is symmetric in the weights v.
- Static check on the start states of T3Pc (seeds 941 and 944, local CPU, double precision): the stray-field energy
  of the corrected model is 4.3·10⁻⁶ K_d V_box per particle above that of X3, in all 8 particles (1.1 % of the stray
  field effect of the staircase, −380·10⁻⁶ per particle; 3.17). Field difference: RMS 3.0·10⁻⁴ M_s, maximum
  2.5·10⁻³ M_s over the magnetic cells.

**Runs.** jobs/80_X3.jobs: X3_s941 … X3_s944, the job line of T3Pc with only --staircase fine, start from the prepared
state of D15 of the same seed. On box 2 beside the seeds 945 … 952 (own queue spool3, 2 slots).

### 3.19 Goal of Chris in three proofs; first result of X3: the correction equals the exact fine stray field in the dynamics (n 2); X3 for all 12 seeds (2026-10-20)

**Goal of Chris (2026-10-20, about 13 UTC).** The roughness of m is given by the FD mesh and is not changed. The
stray field shall be computed as on a converged mesh, ideally without steps, so that the stray field gives no extra
pinning fields from the staircase. That m at the surface depends on the roughness is a separate point. Proofs:
- N1: no mesh dependence at a fixed morphology.
- N2: the stray field is correct statically, i.e. as a finer particle (or a smooth particle) gives it.
- N3: the same dynamically.
Thus the dynamic proof of the stray field compares T3Pc with X3 (same m, exact field of the fine staircase), not with
D15; D15 against X3 is the resolution of m on the fine staircase (the separate point).

**First X3 result** (seeds 941, 942; X3_s943, X3_s944 running): P of X3 1.746 and 1.075 W/cm³, equal to T3Pc;
ln P(T3Pc) − ln P(X3) = −0.0001 and −0.0003. With n 2 this is not yet a result of N3.

**Decision of Chris: the 8 more seeds (945 … 952) run on, and X3 runs for all of them** (jobs/85_X3n.jobs, queue
spool3 after X3_s943/944; each X3 starts when the prepared state of D15 of its seed exists). Options given: A stop the
seeds and finish X3 at n 4; B let them run (N1 and the resolution part with n 12); B with X3 for all seeds was the most
complete one. Cost about 4 $ for box 2 (estimate); a top-up of the credit before about 21 h UTC is necessary.

**Evaluation** (eval_dynamic.py --design B, extended): the arms T3, T3Pc, D15, T15g and X3 with their geometric
means, the reading of D1/K2 for d, and a table of the paired effects with mean, SE, 95 % interval (t distribution)
and the factor on P: N1 = ln P(T15g) − ln P(T3), N3 = ln P(T3Pc) − ln P(X3), Δ_g, Δ_c, m = ln P(D15) − ln P(X3) and
total = ln P(T3Pc) − ln P(D15); --plot gives a figure of these effects. N2 comes from the static results (3.9, 3.18).
- Reading of N1 and N3 (criterion K4, **decision of Chris 2026-10-20**): "within ±0.15" if the 95 % interval of the
  mean lies inside ±0.15, "outside ±0.15" if it lies outside, else "not determined". Reason: K2 needs |mean| at most 2 SE; for an
  effect with a very small SE (N3 about 10⁻⁴) a tiny but real offset would then give "not determined".
- State with the data of 941 … 944 (X3 n 2): N1 +0.083 ± 0.025 (SE, n 4), interval +0.005 … +0.160 (not determined
  by K4, the upper end is above 0.15); N3 −0.0002 ± 0.0001 (n 2); total −0.014 ± 0.115 (n 4).

### 3.20 Box 3 for the seeds 951 and 952; code state of the runs (2026-10-20, order of Chris)

- About 13:15 UTC, Chris: "miete eine weitere billige gpu". Box 3: vast.ai instance 56895734 (offer 47235644,
  RTX 3090, 0.165 $/h, label staircase-box3), image pytorch/pytorch:2.7.1-cuda12.6-cudnn9-devel, code from a git
  bundle of this branch at 50c36c9, packages with the versions of box 2 (numpy 2.2.6, scipy 1.17.1, statsmodels
  0.15.0, torchdiffeq 0.2.5, xitorch 0.5.1, pyvista 0.49.1, setproctitle 1.3.8; torch 2.7.1+cu126).
- Moved to box 3: the seeds 951 and 952 complete (D15, T3, T3Pc, T15g, X3; jobs/90_box3_D15.jobs,
  jobs/91_box3_arms.jobs, the same job lines as on box 2). On box 2 the queue spool2 (D15 949 … 952, 2 slots) was
  stopped (kill of its dispatcher only); its runs D15_s949 and D15_s950 go on (checked). The queued copies of the
  T and X3 jobs of 951/952 on box 2 (spool, spool3) will wait for a prepared state that does not come there; they are
  stopped by hand (note /root/runs/MOVED_s951_s952_to_box3.txt on box 2).
- Code state (provenance): the runs on box 2 that started before 3.18 record git commit 82718f9. The X3 runs and the
  runs that start after 3.18 on box 2 record "82718f9-dirty": the working tree of box 2 has the two files of 381669b
  (magnumnp/field_terms/demag_staircase.py md5 65fd939f…, studies/beta_bpeak/run_loops.py md5 0413a3c5…), which add
  only the option --staircase fine (study tests 118 passed, 3 skipped). Consequence: a run that started clean cannot
  resume after a crash without --allow_new_commit. Box 3 records the clean commit 50c36c9.

### 3.21 N3 with the seeds 941 … 944: within ±0.15 (K4); interim, the full set follows (2026-10-20, 14:15 UTC)

- X3 P (W/cm³) 1.75 / 1.08 / 0.611 / 1.57 against T3Pc 1.75 / 1.07 / 0.611 / 1.56. N3 = ln P(T3Pc) − ln P(X3):
  −0.0001 / −0.0003 / −0.0001 / −0.0066, mean −0.0018 ± 0.0016 (SE, n 4), 95 % interval −0.0069 … +0.0033: reading
  "within ±0.15" (K4). The largest difference (seed 944) is 0.7 % of P.
- The high loss of T3Pc_s944 (3.15) is the same in X3. It comes from the fine staircase with the m of dx 3, not from
  the correction. m = ln P(D15) − ln P(X3): +0.074 / +0.257 / +0.019 / −0.302, mean +0.012 ± 0.116 (n 4).
- Credit: Chris topped up (28.85 $ at 14:14 UTC).

### 3.22 Decisions of Chris for the way to the smooth sphere; dynamic check at q = 4 (2026-10-20, about 14:40 UTC)

**Question of Chris:** what do q = 4, 8, 16 cost? Answer: in production (T3Pc) the window term touches only the cut
cells (8736 at q 2, 13760 at q 4, d/l_ex 96), its field per call took 1.6 s (q 2) and 1.7 s (q 4) on the CPU, and
its set-up grows about as q³ (1.4 s, 12.6 s; estimate q 8 about 2 min, q 16 about 15 min). The proof against the
exact field X3 needs an FFT on (48 q)³ in each step: 96³ (q 2, measured about 1 … 1.5 h per run on a shared GPU),
192³ (q 4, estimate 6 … 10 h), 384³ (q 8, dynamically not practicable, statically once), 768³ (q 16, not practicable).
Statically sufficient (answer to Chris): the dynamics uses the stray field only through its operator; at q 2 a static
energy error of about 1 % of the staircase effect gave a dynamic error of at most 0.7 % in P (N3, 3.21). Condition:
the static checks cover the operator (field per cell in the surface layer, several real states, the sign check R3).

**Decisions of Chris:**
1. Dynamic check at q = 4 with 4 seeds (941 … 944): X3 and T3Pc with q = 4 (X3q4, T3Pc4), start from the prepared
   state of D15, budget 4 … 6 $ accepted. Box 4 (vast.ai 56907720, RTX 3090, 0.179 $/h); jobs/95_q4.jobs (the job
   lines of X3 and T3Pc with only --staircase_q 4). Geometry check: the V1 cells of q 4 (67968) and q 8 (69024) have
   no contacts between particles at the diameter of box_edge.
2. Then q = 8 statically.
3. Then the FE comparison, statically, with magnum.pi. The code is local (~/Documents/code/mgnumpi); Chris has it
   from its author and helped with it; it is never published. Only our own scripts and aggregate results go into this
   repository.
For 2 and 3 a plan follows for the go of Chris.

**Change (Chris, 2026-10-20, about 15:55 UTC):** the dynamic check at q = 4 waits until all runs at q = 2 (3.19 …
3.21, boxes 2 and 3) are done. Box 4 is stopped (vast.ai stop; it never ran: the image was still loading). No data
is on it. While it is stopped, only its storage costs (about 0.05 $/h).
**Change (Chris, 2026-10-20, about 16:15 UTC):** box 4 destroyed (go of Chris; no data on it). For the q = 4 check a
new box is rented after the end of q = 2.
**Decision (Chris, 2026-10-20, about 16:30 UTC):** the T arms of the seeds 949 and 950 stay on box 2 (no move to
box 3: shorter end time, but no cost saving). Proposal to Chris: destroy box 3 when its own runs are done, after
download and md5, only with his go.
**Box 3 (2026-10-20, about 18:05 UTC):** all 10 runs of the seeds 951 and 952 done (D15, T3, T3Pc, T15g, X3). Data
downloaded (88 files, 102 MB, md5 verified, local folder remote_20261020_box3/, not in git). Box 3 destroyed (go of
Chris). Box 2 runs the rest.

### 3.23 Interim with n 6 (seeds 941 … 944, 951, 952): N1 and N3 within ±0.15; why Δ_g, m and total scatter (2026-10-20, about 18:35 UTC)

Values in ln P, mean ± SE, 95 % t-interval (eval_dynamic.py, design B; interim, the reading of the full set follows
at n 12):
- N1 = +0.084 ± 0.016 (n 6), +0.044 … +0.125: within ±0.15 (K4). It is not zero: at fixed morphology the finer m
  increases P by a factor of about 1.09.
- N3 = −0.0025 ± 0.0016 (n 6), −0.0067 … +0.0018: within ±0.15.
- Δ_g = +0.51 ± 0.25 (n 6); Δ_c = +0.40 ± 0.13 (n 6); total = −0.19 ± 0.33 (n 6), −1.04 … +0.65;
  m = +0.13 ± 0.18 (n 11).

Because N3 is nearly zero, total is nearly −m for each seed: the full deviation of T3Pc from D15 is the difference
between D15 and X3. These two arms have the same stray-field geometry (the 5-nm staircase), but different m and
surface material (fine m on its own voxels against coarse m with fraction weights). It is not an error of the stray
field.

**Finding 1, seed 951: D15 changed its state during the cycles.** Loop energy of D15_s951 per cycle: 0.066 J/m³
(cycle 2), 0.13 (3), 0.17 (4), then 0.16 J/m³ constant. In the end state of particle 2, 207 cells (0.17 % of the
particle) differ from the prepared state (m · m_ref < 0.9), all in the interior (r/R 0.57 … 0.75). No cell is
reversed. After the change the shares of particles 1 and 2 are about 10 times larger (phase lag 10.6 and 11.5 mrad).
T3, T3Pc, T15g and X3 started from the same prepared state and stayed in it (overlap 0.997 … 1.000 per particle).
Thus for seed 951 D15 measures a different state. This one event gives +1.70 in Δ_g and −1.71 in total.

**Finding 2, seed 952: same state, different response.** D15 and X3 are in the same state (overlap at least 0.999
per particle), but P(X3)/P(D15) = 1.8, and the share of particle 0 is 3.7 times larger in X3. The loss at 9 mT is a
phase lag of 0.3 … 4 mrad (3.16). Its size depends on the local pinning of the internal structure, and this pinning
is different for fine m on its own 5-nm staircase and for coarse m with fraction material in the field of the same
staircase.

Not checked yet: X3_s949 (6.3 W/cm³) and X3_s950 (3.9 W/cm³; D15_s950 2.4 W/cm³, other pattern of the particles).

**Reading of the interim (agreed by Chris, 2026-10-20, about 18:50 UTC):**
- N1 and N3 test the goal of 3.19: the stray field for the given roughness of m. Δ_g, m and total test a different
  question: a finer model as a whole, in which also m at the surface becomes finer. The roughness of m is a separate
  point (Chris, 3.19). The large deviations of Δ_g, m and total do not come from the stray field.
- Total: the target is 0. The mean −0.19 ± 0.33 (n 6) cannot be distinguished from 0, but because of the scatter it
  is also not shown to be near 0 (95 % interval −1.04 … +0.65).
- Δ_g: there is no target 0. Δ_g is the effect of the finer staircase itself. Its size (+0.51, a factor of about 1.7
  on P) shows that the step size of the surface has a strong effect on the loss.

**Consequence for the goal:**
- Shown (interim): the correction gives the exact field of the 5-nm staircase (N3), and at fixed morphology the mesh
  dependence is about 9 % in P (N1).
- Not shown yet: the 5-nm staircase is still a staircase. Δ_c = +0.40 ± 0.13 shows that the step size in the stray
  field still has a strong effect. Whether the stray field acts as the field of a smooth sphere is shown only by the
  series in q at the same m (T3Pc and X3 with q 4 against q 2, 3.22) and by the static FE comparison.
- D15 changes the roughness of m at the same time. Thus D15 is not the correct reference for the goal; total and
  Δ_g stay as information.

### 3.24 FE reference of the smooth spheres: first local tests, difference method, night run (2026-10-20, about 20:15 UTC)

**Decisions of Chris.** (1) The test "uniform m per particle against the analytic sphere" is withdrawn: every body
with cubic symmetry has the isotropic self tensor I/3, so uniform m sees only the volume; this is reading R1 (3.7, 3.8)
and the β(B_peak) record. The same holds for the FE test case "m uniform per particle" of 3.22; it is replaced by real
states. (2) magnum.pi: the private copy of Chris is the local git repository ~/Documents/code/mgnumpi/magnum.pi_3.13.3_src
(no remote, kept off GitHub on purpose); clone ~/Documents/code/magnum.pi, run in the docker image
magnum.pi.dev-3.13.3; unit tests of the stray field and the periodic meshes: 25 passed. (3) Night run of PLAN.md
section 8, local CPU. (4) Box 2 is destroyed after its runs, the download and the md5 check.

**Setup.** fe/fe_states.py (FD side, magnum.np), fe/fe_demag.py (FE side, magnum.pi, P2 or P1 potential on a periodic
box mesh of n³ cubes with 6 tetrahedra each, mean field zero), fe/fe_eval.py. M per tetrahedron: mean (P1) or L2
projection onto P1 (P2), from the centroids of sub-tetrahedra; this gives the exact load and energy for the given
sources. Two changes against the magnum.pi defaults: the energy is assembled as ½ μ0 ∫ M · ∇u, because DemagField.E
needs m in the P1 space; the box mesh is the periodic mesh of Firedrake, because glueplex.PeriodicMesh does not accept
the facet orientations of fd.BoxMesh. apptainer cannot run on vast.ai containers (no user or mount namespaces, no
CAP_SYS_ADMIN, no /dev/fuse; checked on box 2), so FE runs only locally or on a VM.

**Results, state X3_s951 (energies in 10⁻³ K_d V_box; mesh 48³ = the FD cells).**
- Uniform m, voxel body (target φ_vox (1 − φ_vox) / 3 = 82.517): FD 82.517 (deviation 1·10⁻⁸); FE P1 77.06 (−6.6 %),
  P2 81.61 (−1.1 %).
- Real state, voxel body: FD 4.6520 (exact for cubes); FE P1 0.993 (−79 %), P2 3.858 (−17 %). Cause: the coarse m
  jumps on every inner cell face; the field of these charge sheets varies on scales far below 10 nm, and an element of
  10 nm cannot resolve it. FD Newell integrates it exactly.
- Difference method, check on the 5-nm staircase: FD E_X3,q2 − E_vox = −2.1487; FE P2 −2.1953 (ε = −2.2 % of S;
  with M as the mean per tetrahedron −2.5 %); FE P1 −0.874 (ε = +59 % of S). Thus only P2 is usable.
- Smooth spheres (P2): E_smooth(FE) with 64 / 512 points per tetrahedron 0.962 / 1.057; with 4096 / 32768 points in
  the 108 432 tetrahedra at the surface 1.0886 / 1.0961. E_s = 4.6520 + 1.0961 − 3.8576 = 1.8905. Thus X3 at q 2 lies
  0.613 above the smooth spheres (0.29 S), T3Pc at q 2 lies 0.017 above X3 at q 2 (0.008 S). One state; the night
  run gives 6 states and q 4.
- Continuous m (trilinear, unit length) in the smooth spheres, FE P2 (not checked by an exact value): 1.118, that is
  0.36 S below E_s. The roughness of m therefore has an effect of the same order as the remaining staircase; a
  reference with continuous m would mix both (Chris: the cleanest separation for the test of the correction is the
  same m; the continuous m belongs to the second question, the whole model against the continuum).

**Cost (local laptop, 8 cores).** FD energy with the kernel in the cache: below 1 s. FE P1 48³ (0.66 M tetrahedra):
3 s solve, 1.5 GB. FE P2 48³: 34 s solve, 4.5 GB; sampling 1 min (512 points per cell), 3.5 min or 20 min (4096 or
32768 points at the surface). P2 on 64³: about 10 GB (estimate); 96³: about 30 GB (estimate), not local.
