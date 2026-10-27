# Literature search: "V2" volume-preserving deformed boundary cells

Date: 2026-10-19. Status: DONE (search closed 2026-10-19).

## VERDICT (short)
1. V2 (constant-volume deformed boundary cells, moments kept, FFT far field + exact near-field
   correction, transport partition) was NOT found in any field searched. No close equivalent.
2. Closest: Donahue & McMichael 2007 [H1, H17] and García-Cervera et al. 2003 [H3] (same far/near split,
   but cut cells with volume-fraction moments); pFFT/AIM [H4, H22] (same algebraic scaffold);
   Lebecki 2008 [H8] (edge correction + periodic kernel); Yurkin 2006 [H30] (global volume match).
3. No source compares constant-volume deformed cells with volume-fraction cells. Evidence exists only
   for "volume fraction alone" vs "volume fraction + near-field correction" (the correction does the work).
4. Gaps: see "Citation-list screening summary" and the final report.

## Idea under search (V2)
Each voxel cell keeps volume dx^3 and moment Ms*dx^3. Only the shape of the
boundary cells is deformed (volume-preserving) so that the union of cells is the
smooth body. Stray field = normal FFT/Newell convolution (far field) + local
near-field correction (exact coupling of deformed cells minus coupling of cubes,
for pairs in a small window). Equivalent: staircase surface charges are moved
onto the smooth surface; each cell stays charge-neutral and keeps its dipole
moment. Partition of the body into cells of prescribed volume: least-transport
assignment (power diagram / semi-discrete optimal transport).

## Legend
- [FT] full text read; [ABS] only abstract/landing page read; [SNIP] only a
  search snippet or a secondary description (not cited as evidence).

## Queries run
(see the log below; each line = one query, engine = web search tool)

## Hits

(H-numbers below; query log Q1-Q35 interleaved)

### Query log (part 1)
Q1  finite difference micromagnetics staircase boundary cells volume-preserving deformed cells demagnetization (web, extended)
Q2  Donahue McMichael 2007 "Micromagnetics on curved geometries using rectangular cells" (web)
Q3  precorrected FFT micromagnetics magnetostatic field arbitrary geometry near-field correction (web)
Q4  García-Cervera Gimbutas E "Accurate numerical methods for micromagnetics simulations with general geometries" (web)
Q5  Lepadatu "Effective field model of roughness in magnetic nano-structures" (web)
Q6  DDA shape error boundary dipoles "weighted discretization" Piller edge meshes (web)
Q7  Zhu Liu Yurkin 2023 "Reducing shape errors in the DDA using effective media" (web)
Q8  DDA boundary dipoles shifted position centroid surface normal conformal cells shape error (web, extended)
Q9  Phillips White 1997 precorrected-FFT (web)
Q10 OpenAlex: citers of Donahue&McMichael 2007 (W2130474120, 44 citers), of García-Cervera et al. 2003
    (W2128873244, 64 citers), of Lepadatu 2015 (W2211305628, 21 citers). All titles screened (list below).

## Hits — micromagnetics (contrast set, checked for V2)

H1 [FT, NIST preprint jaggies_20070117.pdf, 3 pages]
M. J. Donahue, R. D. McMichael, "Micromagnetics on curved geometries using rectangular cells:
error correction and analysis", IEEE Trans. Magn. 43, 2878 (2007), doi:10.1109/TMAG.2007.892865.
https://math.nist.gov/~MDonahue/pubs/jaggies_20070117.pdf
- Coarse cell i gets M_i = m_i v_i Ms (v_i = volume fraction inside the body) -> volume-fraction
  weighting. Averaged coupling L_ij from a fine submesh with indicator eta_i' (subcell inside body);
  far field N_ij (same dipole moment => leading error is quadrupole); near-field correction
  H_corr = sum_{j in local window} (N_ij - L_ij) M_j  (eq. 5). Simplified to a cell-local
  "anisotropy" K_i (eqs. 8-9). 2D thin-film validation (edge-mode frequency vs rotation angle, vortex
  exit field); 3D-compatible by construction.
- Key result: uncorrected staircase AND "diluted edge cells" (Ms*v only, no near-field correction)
  both give ~50 % edge-mode frequency increase with rotation; full correction < 5 %.
  Their explanation: charges are "spread out either over the longer stair-step edge perimeter or
  spread over boundaries between cells with different magnetization".
- V2? NO. Cells are cut (volume fraction), not deformed with constant volume. But the FORMALISM
  (eq. 2: L_ij = average of fine tensor over the subcells assigned to coarse cells i and j) is
  general: [own inference] V2 = the same formalism with the assignment "subcell -> coarse cell"
  taken from a volume-preserving partition instead of geometric containment, and v_i = 1.

H2 [FT, arXiv 1508.04391, method section read]
S. Lepadatu, "Effective field model of roughness in magnetic nano-structures", J. Appl. Phys. 118,
243908 (2015), doi:10.1063/1.4939093.
- Roughness as a cell-local diagonal "anisotropy" field from (rough minus smooth) demag sums on the
  same mesh (cells removed -> cell counts N_V, N_Vr). Uniform-M approximation (eq. 6).
- V2? NO (no deformation, no moment preservation; local diagonal term only).

H3 [first pass: ABS only; SUPERSEDED -> full text read later, see 'Update H3' below]
C. J. García-Cervera, Z. Gimbutas, W. E, "Accurate numerical methods for micromagnetics simulations
with general geometries", J. Comput. Phys. 184, 37-52 (2003), doi:10.1016/S0021-9991(02)00014-1.
- Abstract: boundary-corrected algorithm for incomplete boundary cells; staircase gives large errors
  and orientation-dependent switching. (Details: see 'Update H3'.)

## Hits — fast-solver analogue (same structure: grid FFT + near-field precorrection)

H4 [FT, 4-page ICCAD'94 version]
J. R. Phillips, J. K. White, "A precorrected-FFT method for capacitance extraction of complicated 3-D
structures", ICCAD 1994; journal version: IEEE Trans. CAD 16(10), 1059-1072 (1997) [journal: ABS only].
https://cecs.uci.edu/~papers/compendium94-03/papers/1994/iccad94/pdffiles/04c_2.pdf
- Steps: project panel charges onto uniform grid (matching potential on test points; at least
  monopole matched), FFT convolution, interpolate back (W^T), precorrect near pairs:
  P_cor = P_ab - W_a H_ab W_b^T (exact minus grid-represented coupling).
- [own inference] V2's near-field correction is a pFFT precorrection in which the "projection" W is
  the identity map deformed cell -> its cube with the same dipole moment.

## Hits — DDA (light scattering)

H5 [FT, arXiv 0704.0038]
M. A. Yurkin, A. G. Hoekstra, "The discrete dipole approximation: an overview and recent
developments", JQSRT 106, 558-589 (2007), doi:10.1016/j.jqsrt.2007.01.034.
- Shape-error remedies listed (Sect. 3.1, 3.4): Lorentz-Lorenz / effective-medium mixing in boundary
  cells with volume fraction f (Evans & Stephens); Weighted Discretization (Piller, Opt. Lett. 22,
  1674 (1997)): boundary cell split by a planar surface into principal/secondary parts, effective
  susceptibility and effective self-term via a boundary-condition tensor (volume-fraction based);
  adaptive dipole sizes (Bourrely et al.); surface-corrected LDR (Collinge & Draine 2004) modifies
  polarizability only. No method that keeps the cell volume and deforms or moves boundary cells.

H6 [ABS only; Optica full text not reachable by script]
C. Zhu, C. Liu, M. A. Yurkin, "Reducing shape errors in the discrete dipole approximation using
effective media", Opt. Express 31(26), 43401-43415 (2023), doi:10.1364/OE.509479.
- EMA (volume-fraction mixing) in boundary dipoles; up to 60x better integral quantities for spheres
  when the dipole is only several times smaller than the diameter. Volume-fraction type, not V2.

## Hits — theory of boundary charges of a lattice of neutral cells

H7 [ABS] S. Sen, Y. Wang, T. Breitzman, K. Dayal, "Two-scale analysis of the electrostatics of
dielectric crystals: emergence of polarization density and boundary charges", Multiscale Model.
Simul. (2026), doi:10.1137/25M1747907. Also: same authors, J. Mech. Phys. Solids (2024),
doi:10.1016/j.jmps.2024.105743 [ABS].
- Continuum limit of a lattice of charge-neutral unit cells gives bulk polarization PLUS a surface
  charge on the boundary; the surface charge depends on the partial cells at the boundary.
- [own inference] Theoretical support that, at fixed cell moments, the macroscopic surface charge is
  set by the shape of the boundary cells -> the V2 lever (reshape boundary cells, keep moments).

### Query log (part 2)
Q11 OpenAlex/S2: open-access locations for GC-G-E 2003, Zhu 2023, Lebecki 2008, Kim&Lipnikov 2016, Holt 2025
Q12 web: García-Cervera "Numerical micromagnetics: a review" -> author page lists PDFs (Stray/stray.pdf)
Q13 OpenAlex works/doi for: Wei 2009, Horley 2011, Rajaram 2013, Kumar 2013, Lai-Li-Lin 2005, Varga 2012,
    Knittel 2012, Abert 2013, Kritsikis 2008, Exl&Schrefl 2014
Q14 OpenAlex search: Riahi curved nanomagnets finite difference; triangular prism cells micromagnetic FFT
Q15 arXiv API: abs:"stair-step" AND magnet*; staircase AND demagnetizing; staircase AND "stray field";
    "finite difference" AND micromagnetic AND curved; voxel AND demagnetization AND surface;
    "edge smoothing" AND magnet*; "cut cell" AND magnet*; "volume-preserving" AND cells AND boundary AND field
    -> only Riahi (1401.5274) and Abert (1204.4302) relevant.

## Update H3 (García-Cervera, Gimbutas, E 2003) -> [FT, author preprint
http://web.math.ucsb.edu/~cgarcia/Stray/stray.pdf, version 10 Aug 2002, read as page images (pp. 1-12)]
- Abstract: staircase by "either eliminating the incomplete cells or replacing them by complete cells.
  Sometimes the magnetizations at the boundary cells are weighted by the volume of the sample in the
  corresponding cell. We show that this leads to large errors in the computed exchange and stray fields."
- Method (2D thin film, Sect. 3.2): cut boundary cells of arbitrary shape (cell ∩ body); exact mutual
  demag tensors D_ijkl by boundary integrals + adaptive quadrature; split H = H_near + H_far; near field
  = 3x3 nearest neighbours by direct summation with exact cut-cell tensors (eq. 29); far field: boundary
  cells replaced by full rectangles, magnetization "rescaled in such a way that the net magnetization
  remains unchanged" (= volume-fraction weighting, same moment), FFT. They cite Phillips & White
  pFFT (their ref. [18]) as "a similar method". Also a second-order ghost-cell exchange correction.
- Validation: Arrott hexagon (1 um x 0.25 um x 20 nm), 0° vs 45° grid rotation; corrected switching
  path independent of orientation; corrected 100x100 grid = staircase 200x200 grid (hysteresis loop).
- Fig. 12 (hysteresis loops, checked as image): curves = no corrections 45° (100^2, 200^2 grid),
  0° 200^2, boundary-corrected 45° (100^2, 200^2). There is NO volume-weighted-only curve -> the
  abstract's statement that volume weighting "leads to large errors" is an assertion, not a shown result
  in this preprint.
- V2? NO: cut cells with reduced volume, moment of cut cell = Ms * cut volume. 2D only.
- Also summarised in C. J. García-Cervera, "Numerical micromagnetics: a review", Bol. Soc. Esp. Mat.
  Apl. (2007) [FT, author PDF Garcia-Cervera_ReviewSEMA.pdf, Sect. 3.2].

H8 [FT, kops.uni-konstanz.de repository copy]
K. M. Lebecki, M. J. Donahue, M. W. Gutowski, "Periodic boundary conditions for demagnetization
interactions in micromagnetic simulations", J. Phys. D 41, 175005 (2008),
doi:10.1088/0022-3727/41/17/175005.
- Real-space Ewald-like PBC (1D, along z) in OOMMF. For rods of larger radius they COMBINED the PBC
  with the Donahue-McMichael edge correction ("correction factors calculated by considering finer mesh
  (cell subdivision 64 x 64 x 1) interacting over the range of seven local cells"); nucleation fields
  agree with Aharoni-Shtrikman theory. -> published precedent that a local near-field edge correction
  coexists with a periodic demag kernel (1D PBC, 2D cross-section). Not V2.

H9 [FT, arXiv 1608.07309]
E. Kim, K. Lipnikov, "The mimetic finite difference method for the Landau-Lifshitz equation",
J. Comput. Phys. 328, 109-130 (2017) [volume/pages from OpenAlex; online 2016-10-16], doi:10.1016/j.jcp.2016.10.016.
- Exchange on general polyhedral (also distorted) meshes. Stray field (Sect. 3.5): per-element surface
  charges m_E·n_E on each polyhedron's boundary (eq. 3.25), but implemented with FFT ONLY on uniform
  cuboid meshes. No deformed-cell FFT/near-field scheme. Not V2.

H10 [FT, arXiv 1401.5274]
H. Riahi, F. Montaigne, N. Rougemaille, B. Canals, D. Lacour, M. Hehn, "Energy levels of interacting
curved nanomagnets in a frustrated geometry: increasing accuracy when using finite difference methods",
J. Phys.: Condens. Matter 25, 296001 (2013), doi:10.1088/0953-8984/25/29/296001.
- OOMMF, 2D kagome nanomagnets. "Grey-scale mask": moment of a boundary cell proportional to the area
  fraction inside the shape (realised by thickness modulation). Spread of first-neighbour coupling and
  ground-state splitting "drastically reduced"; 5 nm grey ~ 2 nm black/white; reference FE code
  FEELLGOOD. -> evidence FOR volume-fraction weighting on inter-magnet (far-field dominated) energies.
  Not V2.

H11 [ABS] D. Wei, S. Wang, Z. Ding, K.-Z. Gao, "Micromagnetics of ferromagnetic nano-devices using the
fast Fourier transform method", IEEE Trans. Magn. 45, 3035 (2009), doi:10.1109/TMAG.2009.2017260.
- Analytic demag matrices for cubic AND right-angle triangular prism cells, FFT, thin-film devices
  with arbitrary in-plane geometry. Boundary cells = half cells (triangular prisms) -> cut-cell type,
  not volume-preserving [inference from abstract; full text not read].

H12 [ABS] S. Rajaram, D. K. Karunaratne, S. Sarkar, S. Bhanja, "Study of dipolar neighbor interaction on
magnetization states of nano-magnetic disks", IEEE Trans. Magn. 49, 3129 (2013),
doi:10.1109/TMAG.2013.2249051. Application paper (cites D&M); no new discretisation per abstract.

H13 [ABS] P. P. Horley et al., "Numerical simulations of nano-scale magnetization dynamics", InTech
book chapter (2011), doi:10.5772/23745. Review-type chapter; no V2 per abstract.

H14 [ABS] M.-C. Lai, Z. Li, X.-B. Lin, "Fast solvers for 3D Poisson equations involving interfaces in a
finite or the infinite domain", J. Comput. Appl. Math. 191, 106 (2006), doi:10.1016/j.cam.2005.04.025.
Immersed-interface Poisson solver (no abstract in OpenAlex; title only). Not cell-moment based.

H15 [ABS] S. Holt, A. Petrocchi, M. Lang, S. A. Pathak, H. Fangohr, "Discretization anisotropy in finite
difference micromagnetic simulations", IEEE Trans. Magn. (2025), doi:10.1109/TMAG.2025.3604525.
Exchange-stencil anisotropy; remedies: stencil choice and smaller cells. Not V2.

H16 [ABS] E. Kritsikis, J.-C. Toussaint, O. Fruchart, H. Szambolics, L. D. Buda-Prejbeanu, "Fast
computation of magnetostatic fields by nonuniform fast Fourier transforms", Appl. Phys. Lett. 93 (2008) [article number not verified], doi:10.1063/1.2995850. NUFFT for FEM charges: true-surface
charges with N log N cost. Different discretisation (FEM), not V2.

### Query log (part 3)
Q16 web: "surface charge" correction staircase voxel magnetostatic FFT smooth boundary charge-neutral dipole (extended)
Q17 web: conformal FDTD volume-preserving deformed cells enlarged cell technique Dey-Mittra Zagorodnov
Q18 web: gravity magnetic forward modelling prisms staircase topography sloping top polyhedral terrain
Q19 web: VIE FFT conformal boundary cells non-cubic voxels near-field correction precorrected staircase
Q20 web: Power particles (de Goes et al.); Lévy 2015 semi-discrete OT 3D
Q21 web: micromagnetic FD sphere cubic cells staircase demag correction smooth surface level set (extended)
Q22 web: P3M dipolar (Cerdà et al. 2008); Farjadpour et al. 2006 subpixel smoothing
Q23 OpenAlex: Xiao&Liu 2004 ECT; Zagorodnov 2003 USC; Dey-Mittra; Aurenhammer-Hoffmann-Aronov; Kottke 2008;
    Smith-Robertson-Milbert 2001; Heck condensation; AIM (Bleszynski 1996); Gallouët-Mérigot
Q24 web: Heck 2003 Helmert condensation; Plouff 1976 polygonal prisms
Q25 web: Meng et al. 2025 JMMM (powder, cites D&M) -> not found by web tool; no abstract in OpenAlex
Q26 web (extended): "edge correction" OR "staircase correction" micromagnetic FFT same volume deformed
    hexahedral polyhedra -> nothing new
Q27 NIST pages: Donahue abstracts + talks list -> 3 talk PDFs downloaded and read

## Donahue & McMichael talk slides (more detail than the 3-page paper)

H17 [FT of slides, not peer-reviewed]
M. J. Donahue, R. D. McMichael, "Micromagnetics on curved geometries using rectangular cells: error
correction and analysis", talk, Bonn 2006 (errata 2007-01-16):
https://math.nist.gov/~MDonahue/talks/bonn2006-staircase-errata20070116.pdf ; and Joint MMM/Intermag
2007 talk GE05: https://math.nist.gov/~MDonahue/talks/joint2007-GE05-staircase.pdf
- Rotated strip: staircase charges "Total charge preserved, but strength reduced, area increased".
- Uniformly magnetised rotated Py square (350 nm x 350 nm x 5 nm): staircase field error shows
  "1st order convergence" in cell size.
- Dilution ("Use full cell so all cells have same geometry, but reduce Ms so far-field demag is
  correct", from Porter & Donahue 2001): far-field error of a 1x1x1 cube with Ms/2 vs a 1x1x0.5 brick
  with Ms: relative error ~0.625/x^2, absolute ~0.05/x^5 (x = offset/cell size).
- "Infinite plate" slide: for the reduced-|M| model the self field is |H_demag| = (t/Δ) Ms instead of
  Ms -> local deficit; fix = local anisotropy field -(1 - t/Δ) Ms (m·e_z) e_z.
  10x10x1 oblate spheroid: RMS demag error 118 % (staircase disk) -> 29 % (thickness corrections).
- "Edge cell adjustment" (analytic model, straight edge): the energy difference between partially
  filled cells and diluted cells is proportional to v(1-v) ("No correction for empty cells (v = 0)",
  "No correction for filled cells (v = 1)"). The prefactor came out garbled in pdftotext; only the
  v(1-v) proportionality is quoted here.
- Edge-mode plot includes "Far field correction" (dilution only) = no improvement, "Far field +
  empirical" and "Far field + local subgrid" = strong improvement.
=> Literature argument AGAINST pure volume-fraction weighting: right far field (moment), wrong local
   self/near field (∝ v(1-v)); the near-field correction does the work.

## Electromagnetics analogue: subpixel smoothing (FOR/AGAINST evidence)

H18 [FT, author PDF https://math.mit.edu/~stevenj/papers/FarjadpourRo06.pdf, 3 pages]
A. Farjadpour, D. Roundy, A. Rodriguez, M. Ibanescu, P. Bermel, J. D. Joannopoulos, S. G. Johnson,
G. W. Burr, "Improving accuracy by subpixel smoothing in the finite-difference time domain",
Opt. Lett. 31(20), 2972-2974 (2006).
- Staircased (pixelised) interfaces reduce FDTD convergence from quadratic to linear in Δx.
- Smoothing changes the problem; to restore quadratic convergence the smoothing must have zero
  first-order effect (perturbation criterion): anisotropic tensor eps^-1 = P<eps^-1> + (1-P)<eps>^-1,
  P = n n^T (surface normal). "The simplest is to use the scalar mean <eps> for all components, which is
  incorrect for the surface-normal fields." All other schemes converge linearly; two published
  anisotropic schemes were "actually worse than no smoothing" for the TE case.
- [own inference] Analogy: volume-fraction weighting of Ms is a scalar-mean smoothing; it fixes total
  moment (zeroth order) but not the normal-field / surface-charge physics at first order. V2 instead
  puts the charge on the true surface with the true normal. No source states this for magnetostatics.
H19 [ABS] C. Kottke, A. Farjadpour, S. G. Johnson, "Perturbation theory for anisotropic dielectric
interfaces, and application to subpixel smoothing of discretized numerical methods", Phys. Rev. E 77,
036611 (2008), doi:10.1103/PhysRevE.77.036611.

## Conformal FDTD / cut cells

H20 [ABS] T. Xiao, Q. H. Liu, "Enlarged cells for the conformal FDTD method to avoid the time step
reduction", IEEE Microw. Wireless Compon. Lett. 14, 551 (2004), doi:10.1109/LMWC.2004.837384.
Small irregular cut cells are enlarged (borrow area from neighbours) for stability. Comment by
Schuhmann, Zagorodnov, Weiland, doi:10.1109/LMWC.2005.861553 [title only].
H21 [ABS] I. A. Zagorodnov, R. Schuhmann, T. Weiland, "A uniformly stable conformal FDTD-method in
Cartesian grids", Int. J. Numer. Model. 16, 127 (2003), doi:10.1002/jnm.488. PEC, 2nd order, no time
step reduction.
- Both: cut-cell (area/volume-fraction) geometry; cell enlargement redistributes area for STABILITY,
  not to keep a prescribed volume per cell. No volume-preserving deformed-cell variant found.

## Fast integral-equation analogues

H22 [ABS] E. Bleszynski, M. Bleszynski, T. Jaroszewicz, "AIM: Adaptive integral method for solving
large-scale electromagnetic scattering and radiation problems", Radio Sci. 31, 1225 (1996),
doi:10.1029/96RS02504. Arbitrary basis functions replaced by grid point sources that reproduce the far
field (moment matching); far field via Toeplitz/FFT, near field exact (Galerkin). Same structure as
pFFT [H4]. [own inference] V2 = AIM/pFFT with the trivial projection "deformed cell -> cube with the
same moment" (dipole-order matching only, monopole = 0).
H23 [ABS] J. J. Cerdà, V. Ballenegger, O. Lenz, C. Holm, "P3M algorithm for dipolar interactions",
J. Chem. Phys. 129(23), 234104 (2008), doi:10.1063/1.3000389 [verified via OpenAlex]; arXiv:0805.4783. Mesh FFT + real-space short-range part for point dipoles under PBC.
Structure analogue only; no boundary-shape treatment.

## Geophysics (prism = Newell-equivalent)

H24 [ABS] D. Plouff, "Gravity and magnetic fields of polygonal prisms and application to magnetic
terrain corrections", Geophysics 41(4), 727-741 (1976), doi:10.1190/1.1440645. Exact fields of
polygonal prisms; faster/more accurate than lamina integration; more general than rectangular prisms.
Magnetic formulas neglect demagnetisation. Geometry-conforming source bodies (not volume-preserving
deformation of a grid).
H25 [SNIP only, abstract not accessible] D. A. Smith, D. S. Robertson, D. G. Milbert, "Gravitational
attraction of local crustal masses in spherical coordinates", J. Geod. 74, 783-795 (2001),
doi:10.1007/s001900000142. Search snippet: flat-topped prisms inadequate near the station (errors up
to 20 mGal with 30" data); bilinear tops much closer. NOT used as evidence (unread).

## Partition into cells of prescribed volume (the V2 partition step)

H26 [ABS] F. Aurenhammer, F. Hoffmann, B. Aronov, "Minkowski-type theorems and least-squares
clustering", Algorithmica 20, 61-76 (1998), doi:10.1007/PL00009187; conference version SoCG 1992,
doi:10.1145/142675.142747 (abstract read): "there always exists a power diagram whose regions
partition a given ... set into clusters of prescribed sizes, no matter where the sites are taken";
power-diagram assignments = least-squares (L2 transport) assignments; extension to continuous sets:
"power diagrams with prescribed region volumes".
H27 [FT, arXiv 1409.1279, Sect. 3 read] B. Lévy, "A numerical algorithm for L2 semi-discrete optimal
transport in 3D", ESAIM: M2AN 49(6), 1693-1715 (2015), doi:10.1051/m2an/2015055. Theorem 3: for a
density mu and points p_i with prescribed masses summing to mu(Omega), weights W exist with
mu(Pow_W(p_i)) = nu_i; the power-diagram assignment is the optimal transport map; weights from a
concave (dual) maximisation with Newton-type solver; tests up to 1e6 Dirac masses.
H28 [ABS] B. Lévy, "Partial optimal transport for a constant-volume Lagrangian mesh with free
boundaries", J. Comput. Phys. 451, 110838 (2022; OpenAlex year 2021), doi:10.1016/j.jcp.2021.110838; arXiv:2106.03936. Laguerre cells
with prescribed volume, free boundaries; applications fluids, astrophysics.
H29 [ABS] F. de Goes, C. Wallez, J. Huang, D. Pavlov, M. Desbrun, "Power particles: an incompressible
fluid solver based on power diagrams", ACM Trans. Graph. 34(4) (2015), doi:10.1145/2766901. Particles = volumetric power
cells with volume control. Fluid, not field computation.
- No source found that uses a prescribed-volume power diagram to build cells for a magnetostatic or
  electrostatic interaction kernel.
- [own inference, design note] With sites on the cubic lattice and the smooth body as domain, the
  capacity-constrained weights are not constant: interior cells are also slightly deformed (perturbation
  decays inward but is not exactly zero). Any near-field window must cover all cells whose shape
  changed, or the partition must be restricted to a boundary shell (interior weights fixed).

### Query log (part 4)
Q28 web: ADDA "-no_vol_cor" volume correction -> led to Yurkin et al. 2006 (read)
Q29 web: micromagnetic Voronoi tessellation finite volume polyhedral cells FFT -> nothing relevant
Q30 web: DDA/coupled dipole boundary cells deformed conformal volume-preserving/equal-volume -> nothing
Q31 web: demag tensor arbitrary polyhedra / interaction of two polyhedra -> Nielsen et al. 2019 (tool)
Q32 Semantic Scholar: citers of D&M 2007 (40) and Lepadatu 2015 (21); OpenAlex citers of Parker 2000 (10),
    Riahi 2013 (4), Wei 2009 (11)
Q33 web (extended): 3D FD micromagnetics sphere staircase artifact correction -> nothing new
Q34 NIST: Miltat & Donahue chapter "Numerical Micromagnetics: Finite Difference Methods" (2007) refs
    [16-18] for staircase corrections = Parker 2000, Porter & Donahue 2001, GC-G-E 2003 (same set)
Q35 arXiv abs 1303.1403 (Hamm, Renn, Hess)

## Further hits

H30 [FT, arXiv 0704.0035] M. A. Yurkin, V. P. Maltsev, A. G. Hoekstra, "Convergence of the discrete
dipole approximation. II. An extrapolation technique to increase the accuracy", J. Opt. Soc. Am. A
23(10), 2592 (2006), doi:10.1364/JOSAA.23.002592.
- "dipole size correction ... to ensure that the cubical approximation of the scatterer has the correct
  volume; this is believed to diminish shape errors" (global volume match, same convention as V2's
  "sphere with exactly the voxel volume").
- Shape errors "depend upon the position of the particle surface inside the boundary dipole", behave
  like "random noise superimposed upon a smooth variation", asymptotically linear in the dipole size.
  -> relevant to the observed strong dependence of results on the staircase realisation.

H31 [ABS] S. Wang, D. Wei, K.-Z. Gao, "Limits of discretization in computational micromagnetics",
IEEE Trans. Magn. 47(10), 3813-3816 (2011), doi:10.1109/TMAG.2011.2157474 (= "Wang 2011").
Improved FDM-FFT: cubic cells in the body plus prism cells at edges (thin film); results similar to
cubic-only model. Edge prisms = cut cells; not V2 [inference from abstract].

H32 [ABS] K. K. Nielsen, A. R. Insinga, R. Bjørk, "The stray- and demagnetizing field from a
homogeneously magnetized tetrahedron", IEEE Magn. Lett. 10, 2110205 (2019),
doi:10.1109/LMAG.2019.2956895; arXiv:2004.11700. Closed-form field tensor of a uniformly magnetised
tetrahedron (MagTense). Tool for the exact field of deformed (polyhedral) cells; cell-cell averaged
coupling would still need an outer integration.

H33 [ABS] J. Hamm, F. Renn, O. Hess, "Dispersive media subcell averaging in the FDTD method using
corrective surface currents", arXiv:1303.1403, DOI 10.1109/TAP.2013.2292073 (IEEE TAP, per DOI
prefix). Interface discontinuities represented by corrective surface currents; ~2nd order restored.

H34 [title/bibliography only] G. J. Parker, C. Cerjan, D. W. Hewett, "Embedded curve boundary method for
micromagnetic simulations", J. Magn. Magn. Mater. 214, 130-138 (2000),
doi:10.1016/S0304-8853(00)00043-3. No abstract available; 10 citers screened (none V2-like).

## Citation-list screening summary
- D&M 2007: OpenAlex 44 / Semantic Scholar 40 citers. Titles screened. Union adds (S2 only, no DOI,
  NOT read): "Simulations of ferromagnetic nano structures" (2008), "Numerical study of magnetic
  processes: extending the LLG approach from nanoscale to microscale" (2010), "Domain wall behaviour in
  ferromagnetic nanowires with interfacial and geometrical structuring" (2013), "Simulação
  micromagnética ... rugosidade em nanofios de níquel" (2015), "Problem descriptions: Non-uniform fast
  Fourier transform" (2015), NIST MCSD FY2007 summary (2008). Likely theses/reports.
- García-Cervera et al. 2003: OpenAlex 64 citers screened by title; none suggests V2.
- Lepadatu 2015: 21 citers (identical in OpenAlex and S2); all applications/Boris papers.
- Parker 2000: 10; Riahi 2013: 4; Wei 2009: 11 citers screened.
- Meng et al. 2025, JMMM 173350 (powder core loss 1-20 kHz, cites D&M): no abstract obtainable; NOT
  read. Same application field as the magnum.np powder model -> worth a manual check.

H35 [ABS, via NIST abstracts page https://math.nist.gov/~MDonahue/pubs/abstracts.html]
D. G. Porter, M. J. Donahue, "Generalization of a two-dimensional micromagnetic model to non-uniform
thickness", J. Appl. Phys. 89, 7257-7259 (2001). Thinner cells: Zeeman/anisotropy scale with thickness,
pair exchange "scales by a weight dependent on the thicknesses of both cells", demag by scaling the
cell moment plus a local out-of-plane correction; 10x10x1 oblate spheroid more accurate than uniform
thickness. -> origin of the "dilution" model. [own inference] Any cut/diluted-cell model must also
re-weight exchange (GC-G-E Sect. 3.1 does this with ghost cells); V2 keeps all cells identical in volume,
so the exchange stencil needs no volume weights (but the cell-centre spacing near the surface is still
the cubic one).

## Evidence quality note
All quoted effect sizes (D&M: ~50 % vs < 5 % edge-mode shift, 26 % vs 3.2 % vortex-exit change; Riahi:
5 nm grey ~ 2 nm black/white; GC-G-E: corrected 100^2 = uncorrected 200^2) are from single
deterministic simulations per configuration, without statistical error bars (D&M: error bars = field
step only). Zhu et al. 60x is from the abstract only.

## Gaps of this search
- Not read in full: Zhu, Liu, Yurkin 2023 (OA, but the Optica page needs a browser); Wei et al. 2009,
  Wang et al. 2011, Parker et al. 2000, Kritsikis 2008, Exl & Schrefl 2014 (abstract/title only);
  Phillips & White 1997 journal version (4-page 1994 conference version read); Plouff 1976 (abstract).
- Smith, Robertson, Milbert 2001: search snippet only -> not used as evidence.
- Meng et al. 2025, JMMM 173350 (powder, 1-20 kHz, cites D&M): no abstract obtainable -> manual check
  recommended (same application as the magnum.np powder model).
- Six Semantic-Scholar-only citers of D&M (theses/reports, no DOI) not read.
- Semantic Scholar rate limit: GC-G-E citers only from OpenAlex (64), not cross-checked.
- No Google Scholar access (only a generic web search tool, OpenAlex, Semantic Scholar, arXiv API);
  no Chinese-language, patent or thesis-database search. FDTD and geophysics coverage is by abstract.
