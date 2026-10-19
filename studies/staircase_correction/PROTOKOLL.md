# Protokoll: staircase correction of the stray field (magnum.np, Newell solver with full PBC)

Owner: Christoph Vogler (Chris). Branch: `feature/staircase-correction` (from `main` 6ca3b65, 2026-10-19). Folder:
`studies/staircase_correction/`. Literature notes: `studies/staircase_correction/lit/`.

Rules: documents in ASD-STE100 Simplified Technical English, SI units only.
**Chris makes the decisions.** Before you change the model, the method, the test design or the budget, give him the
options with pros, cons, chances and risks, and wait. Engineering fixes (bugs, logging, speed) need no decision, but
tell him.

References: a number "7.x" is a section of the record of the β(B_peak) study, `studies/beta_bpeak/PROTOKOLL.md`
(branch `study/beta-bpeak-pinning`). A number "3.x" is a section of this record.

## 0. Current state (2026-10-19)

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
  with constant volume, idea of Chris) was not found (3.2).
- **Status.** Stage A of PLAN.md (phase 1, static test, V1 with both forms) approved; code in progress.

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
- 2026-10-19: own record for the staircase correction (this file). The sections 7.117, 7.120 and 7.121 of the study
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
