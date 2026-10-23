# Literature notes: staircase correction compatible with Newell FFT + full PBC

Date: 2026-10-19. Status: search complete; final report given in the chat of 2026-10-19.

Constraint: keep the FFT convolution of M with the cell-centred Newell tensor
of cubic cells and the full periodic boundary conditions (image shells plus
compensation of the outer surface charges). Allowed: local terms, cell weights,
translation-invariant kernel change, sparse local correction.

Legend for "read": FULL = full text read; ABS = abstract / summary page only;
SNIP = only a search snippet (not cited as evidence).

## Queries

- OpenAlex API: citers of Donahue & McMichael 2007 (W2130474120, 44 citers), Lepadatu 2015
  roughness (W2211305628, 21), Garcia-Cervera/Gimbutas/E 2003 (W2128873244, 64),
  Chernyshenko & Fangohr 2015 JMMM (W2039452433, 11), Holt TMAG 2025 (W4414221784, 2),
  Holt arXiv 2412.10466 (0), Parker/Cerjan/Hewett 2000 (W2053776092, 10),
  Lebecki/Donahue/Gutowski 2008 (W2060046398, 95), Yurkin et al. 2006 DDA II (W2117289988, 53).
  Raw lists in scratchpad (not kept).
- Web: "composite voxels FFT-based homogenization laminate staircase"; "precorrected FFT periodic
  near-field correction"; "micromagnetic finite difference curved boundary staircase correction
  sphere FFT"; "gravity magnetic forward modeling prism FFT staircase topography Gauss-FFT";
  "discrete dipole approximation weighted discretization boundary dipoles shape error";
  "Mayo method fast Poisson solver irregular region".

## Hits by family

### F1 Local near-field corrections (exact tensors of cut / polyhedral cells)

Known (no new content found): Donahue & McMichael 2007 (IEEE TMAG 43, 2878). NEW material checked:
- Donahue M.J., McMichael R.D. Slides "Micromagnetics on curved geometries using rectangular
  cells: error correction and analysis", Micromagnetics Workshop Bonn 2006 (errata 2007-01-16),
  https://math.nist.gov/~MDonahue/talks/bonn2006-staircase-errata20070116.pdf. READ: FULL (text
  layer; slides 16-32 also as images). Content beyond what we know from the paper:
  (i) p.16-17 infinite plate: reduced |M| = (t/dx) Ms gives the right far field but the wrong
  near field, |H_demag| = (t/dx) Ms; fix = local anisotropy field
  H_anis = -(1 - t/dx) Ms (m.e_z) e_z ("thickness correction").
  (ii) p.19-20: 10 x 10 x 1 oblate spheroid (thin), uniform in-plane M, cells cut through the
  thickness: RMS error of the demag field 118 % with plain staircase, 29 % with the cell-local
  thickness corrections only.
  (iii) p.21-28 general case: split into Omega_local (near, exact) and Omega_far (FFT with
  volume-modified |M|); near part from the EXISTING demag code on a local refined grid around
  each boundary cell; correction factors K_ij = fine minus coarse field, computed once at
  initialisation; run time H_corr,i = -sum_j K_ij M_j over Omega_local; cost O(N_boundary).
  Cites Garcia-Cervera, Gimbutas & E 2003 for the near/far concept.
  (iv) p.29-30 "pushed" form: if |m_j - m_i| is small (exchange), the correction reduces to a
  local anisotropy K_i = sum_j (|M_j|/|M_i|) K_ij. -> This shortcut ASSUMES m nearly uniform over
  Omega_local; for our wall widths of the order of the cell size the full K_ij form is needed
  (inference from the slide's own condition).
  (v) p.31-32 FMR "ping" test of a thin square, edge and bulk modes. No PBC discussion.
- Donahue M.J., Porter D.G. "High order methods for computing the demagnetization tensor for
  periodic boundaries", MMM 2020 talk P3-02, https://math.nist.gov/~MDonahue/talks/mmm2020-talk.pdf.
  READ: FIRST PAGES only. About accurate PBC Newell tensors (1D/2D periodic sums), not staircase.
- Nielsen K.K., Insinga A.R., Bjork R. (2019). The stray- and demagnetizing field from a
  homogeneously magnetized tetrahedron. IEEE Magn. Lett. 10, 2110205. doi:10.1109/LMAG.2019.2956895,
  arXiv:2004.11700. READ: FULL (text). Analytic POINT-evaluated tensor field of a uniformly
  magnetised triangle face and tetrahedron; validated against COMSOL FEM along three lines
  through one test geometry ("excellent agreement", no error number). Implemented in MagTense.
  Use for us (own inference): a cut cube = polyhedron = union of tetrahedra, so the exact field
  of the cut cell is available in closed form; but Newell is CELL-AVERAGED, so a point-evaluated
  tensor must be averaged over the target cell (numerical quadrature) to be consistent.
- Bjork R., d'Aquino M. (2023). Accuracy of the analytical demagnetization tensor for various
  geometries. JMMM 587, 171245. doi:10.1016/j.jmmm.2023.171245. READ: ABS. Prism and tetrahedron
  tensors: median relative error vs dipole field below 1e-7 at 200 tile radii; below 1e-2 at
  1e4 radii also in single precision. -> analytic tetrahedron tensor is numerically safe.
- Guptasarma D., Singh B. (1999). New scheme for computing the magnetic field resulting from a
  uniformly magnetized arbitrary polyhedron. Geophysics 64(1), 70-74. doi:10.1190/1.1444531.
  READ: ABS. Field outside a polyhedron from surface poles per polygon face: line integrals over
  edges + solid angle.
- Fabbri M. (2008). Magnetic flux density and vector potential of uniform polyhedral sources.
  IEEE Trans. Magn. 44(1), 32-36. doi:10.1109/TMAG.2007.908698. READ: ABS. Closed form, frame
  independent, elementary functions, uses face-edge incidence data.
- Tandon S., Beleggia M., Zhu Y., De Graef M. (2004). On the computation of the demagnetization
  tensor for uniformly magnetized particles of arbitrary shape. Part I: Analytical approach.
  JMMM 271, 9-26. doi:10.1016/j.jmmm.2003.09.011. READ: ABS (unifind page). Shape amplitudes for
  sphere, tube, arbitrary polyhedron and others. Part II: numerical approach, JMMM 271, 27-38,
  doi:10.1016/j.jmmm.2003.09.010. READ: ABS. Fourier-space and real-space numerical methods,
  checked on a cylinder. (Uniform magnetisation only.)
- Liang H., Yan X. (2024). On the equivalence of demagnetization tensors as discrete cell size
  approaches zero in three-dimensional space. J. Appl. Phys. 136, 083903. doi:10.1063/5.0226603.
  READ: ABS. Not staircase; about near-cell tensor contributions and principal value limit.

### F2 Locally refined / non-equidistant grids that keep an FFT

- Ai F., Duan J., Lomakin V. (2024). Periodic micromagnetic finite element method. arXiv:2409.13958.
  READ: FULL (Sect. 3.2). Magnetostatics of FEM nodes via box adaptive integral method (BAIM, a
  variant of precorrected FFT): (1) project non-uniform sources to a uniform Cartesian grid,
  (2) FFT convolution with the PERIODIC Green function (no zero padding), (3) interpolate back,
  (4) sparse near-field correction: subtract the grid result of near sources and add the exact
  contribution. Key point for us: the near correction uses only the FREE-SPACE kernel G0, also
  under PBC, because "only G0 varies drastically within nearby region"; the correction region is
  only extended across the box faces (periodic wrap of neighbour search). Relative error level
  1e-3 used in timing tests; periodic case as fast as non-periodic. Tested on FEM meshes (cube
  unit cell), not on a voxel staircase.
  Transfer (own inference): a sparse near-field correction (Donahue-McMichael type) is
  compatible with our periodic Newell FFT; compute it with the free-space exact kernel and apply
  the minimum-image wrap.
- Phillips J.R., White J.K. (1997). A precorrected-FFT method for electrostatic analysis of
  complicated 3-D structures. IEEE Trans. CAD 16(10), 1059-1072. doi:10.1109/43.662670.
  READ: ABS. Source of pFFT (project to grid, FFT, interpolate, precorrect near interactions).
- Exl L., Schrefl T. (2014). Non-uniform FFT for the finite element computation of the
  micromagnetic scalar potential. J. Comput. Phys. 270, 490-505. doi:10.1016/j.jcp.2014.04.013,
  arXiv:1305.3162. READ: FULL (abstract, intro, Sect. 3 headings and near-field part).
  Single-layer potential of the surface charge m.n by kernel split: smooth part via Fourier
  approximation + NUFFT, singular part as a sparse near-field correction (Gaussian quadrature on
  surface triangles). States that pFFT had not been adapted to micromagnetics "to the authors'
  best knowledge". Free space, FEM. Relevant for F3 (field of a smooth surface charge on a grid).
- Garcia-Cervera C.J., Roma A.M. (2006). Adaptive mesh refinement for micromagnetics
  simulations. IEEE Trans. Magn. 42(6), 1648-1654. doi:10.1109/TMAG.2006.872199. READ: ABS.
  FD + nested, progressively finer rectangular patches that follow walls and vortices; test:
  rectangular sample of infinite thickness. Not aimed at curved surfaces.
- Livshitz B., Boag A., Bertram H.N., Lomakin V. (2009). Nonuniform grid algorithm for fast
  calculation of magnetostatic interactions in micromagnetics. J. Appl. Phys. 105(7).
  doi:10.1063/1.3076048. READ: ABS. Multilevel non-uniform-grid potential representation, linear
  complexity, free space. Different solver.
- Lepadatu S. (2019). Efficient computation of demagnetizing fields for magnetic multilayers
  using multilayered convolution. J. Appl. Phys. 126(10). doi:10.1063/1.5116754. READ: ABS.
  Exact FFT convolution between layers of different thickness / spacing. Refinement only per
  layer (one axis), like magnum.np demag_nonequidistant; no help for a sphere.
- Ai F., Lomakin V. (2024). Fast Fourier transform periodic interpolation method for
  superposition sums in a periodic unit cell. Comput. Phys. Commun. 304, 109291.
  doi:10.1016/j.cpc.2024.109291. READ: ABS. Periodic potential split into near-zone (few images,
  FFT on a fine grid) and far-zone (sparse grid + interpolation); O(N log N); non-uniform sources.
- Kritsikis E., Toussaint J.-C., Fruchart O., Szambolics H., Buda-Prejbeanu L.D. (2008). Fast
  computation of magnetostatic fields by nonuniform fast Fourier transforms. Appl. Phys. Lett. 93,
  132508. doi:10.1063/1.2995850, arXiv:0807.2387. READ: ABS. NUFFT brings N log N to FEM.

### F3 Charge-based corrections (surface charge, face charge)

No paper found that adds the field of the smooth-surface charge minus the staircase face charge
as a sparse term to a Newell FFT (searched: see Queries). Closest building blocks:
- Exl & Schrefl 2014 (see F2): field of the surface charge m.n on a triangulated surface by
  NUFFT + sparse singular near-field correction (free space).
- Ai, Duan, Lomakin 2024 (see F2): same structure with a periodic Green function; near
  correction uses the free-space kernel only.
- Mayo-type embedded boundary integral methods (Mayo 1984, not opened, cited only via search
  summaries; Biros, Ying,
  Zorin "embedded boundary integral solver for Stokes"; Marques, Nave, Rosales 2014 arXiv
  1401.8084 "correction function method" + Mayo). READ: SNIP only (search summaries). They correct
  the right-hand side of a regular-grid FFT Poisson solver locally at the interface from boundary
  integral data; cost about twice an FFT Poisson solve (Marques et al., APS DFD 2011 abstract,
  SNIP). Different solver class (Poisson FD) -> short note only.
- Bruckner F., Ducevic A., Heistracher P., Abert C., Suess D. (2021). Strayfield calculation for
  micromagnetic simulations using true periodic boundary conditions. Sci. Rep. 11, 9202.
  doi:10.1038/s41598-021-88541-9. READ: FULL (text). FD part = differential formulation (potential
  on vertices, central differences) solved in Fourier space; NOT a Newell form-factor kernel; no
  staircase treatment. (magnum.np authors; for reference only.)
- Abert C., Exl L., Selke G., Drews A., Schrefl T. (2013). Numerical methods for the stray-field
  calculation: A comparison of recently developed algorithms. JMMM 326, 176-185.
  doi:10.1016/j.jmmm.2012.08.041, arXiv:1204.4302. READ: FULL (sphere section). Uniformly
  magnetised sphere, cell-centre voxelisation: FFT/tensor-grid energies converge slower than FEM
  "obviously a consequence of the better approximation of the curved surface"; refers to D&M 2007.

### F4 Discrete dipole approximation (DDA)

Transfer note (own inference): DDA on a cubic lattice with FFT is, in the static limit k -> 0,
a cell-based volume-integral method of the same structure as the Newell FFT. It is close to,
not identical with, Newell: "integration of Green's tensor" (IGT) integrates over the SOURCE
cell and evaluates at the target centre, Newell averages over source AND target cell; DDA
also solves a linear polarisability problem, while in micromagnetics M is the state variable.
The DDA "shape error" of a voxelised sphere is the same geometric problem as our staircase.

- Yurkin M.A., Maltsev V.P., Hoekstra A.G. (2006). Convergence of the discrete dipole
  approximation. I. Theoretical analysis. J. Opt. Soc. Am. A 23(10), 2578.
  doi:10.1364/JOSAA.23.002578, arXiv:0704.0033. READ: FULL (arXiv text, Sect. 2.E and WD part).
  Findings: shape error of boundary dipoles adds a term linear in y = |m| k d (d = dipole size)
  to the error; it "reinforces" the surface (linear) error. Analysis of Weighted Discretization
  (WD, Piller 1997): effective susceptibility and effective self-term of each boundary cell
  from the two parts of the cut cube (planar interface). Result: WD keeps a linear term but with a
  smaller constant (no ln y factor). Key statement for us: integrating the Green tensor over the
  two parts of the cut cell separately "would ruin the block-Toeplitz structure of the
  interaction matrix and hinder the FFT-based algorithm" -> they reject it. (In our constraint a
  sparse near-field correction in addition to the FFT is allowed, so this objection does not
  apply to us; own inference, cf. Donahue & McMichael 2007.)
- Yurkin M.A., Maltsev V.P., Hoekstra A.G. (2006). Convergence of the discrete dipole
  approximation. II. An extrapolation technique to increase the accuracy. J. Opt. Soc. Am. A 23(10),
  2592-2601. doi:10.1364/JOSAA.23.002592, arXiv:0704.0035. READ: FULL.
  Findings: Richardson-type extrapolation in y (quadratic fit, data weights y^3). Shape errors
  "can be viewed as random noise superimposed upon a smooth variation" because they depend on
  where the surface cuts the boundary cells; asymptotically linear in y. Recommendation for
  non-cubic shapes: 9 discretisations over [ymin, 4 ymin]; error estimate 2 x SE. Extrapolation
  satisfactory for spheres only for ymin < 0.15 (cubes: < 0.4), at m = 1.5. Volume correction
  (dipole size so that the voxel body has the true volume) used. kD = 10 sphere, 16 dipoles per
  diameter: maximum S11 errors mostly shape errors (ratio max shape/max discretisation error
  about 4), Qext errors mostly discretisation.
  Transfer (own inference): Richardson in dx at a non-locked staircase geometry is unreliable;
  quasi-random shape error must be averaged (many dx values) or removed first.
- Zhu Y., Liu C., Yurkin M.A. (2023). Reducing shape errors in the discrete dipole approximation
  using effective media. Optics Express 31(26), 43401. doi:10.1364/OE.509479.
  READ: ABS (OpenAlex abstract; full text behind JS page, not reached).
  Findings (abstract): effective-medium approximation (EMA) for boundary dipoles; spheres vs
  Lorenz-Mie: integral quantities up to 60 times more accurate when dipole size is only several
  times smaller than the diameter; often comparable to halving the dipole size; for coated
  black-carbon: Qext error 4.7 % -> 0.3 %; 1 % accuracy with 3 times larger dipoles (about 30
  times faster). Regular lattice and FFT kept (cell-local change of polarisability only;
  inference from the abstract: "boundary dipoles").
- Piller N.B. (1997). Influence of the edge meshes on the accuracy of the coupled-dipole
  approximation. Opt. Lett. 22(22), 1674. doi:10.1364/OL.22.001674. READ: ABS. Weighted
  discretisation: convergence "much more regular", accuracy better "by as much as a factor of
  24", no extra computation time. (Method details via Yurkin I: effective susceptibility and
  effective self-term of each boundary cell from a planar cut; lattice and FFT kept.)
- Piller N.B., Martin O.J.F. (1998). Increasing the performance of the coupled-dipole
  approximation: a spectral approach. IEEE Trans. Antennas Propag. 46(8), 1126-1137.
  doi:10.1109/8.718567. READ: ABS. Filtered coupled dipoles (FCD): source with spatial extension,
  "filtered Green's tensors" (1D/2D/3D) reduce aliasing of the discretised scatterer; cross
  section up to 100 times more accurate at the same grid, similar cost. Translation-invariant
  kernel change -> FFT kept. Transfer (own inference): Newell is already a cell-averaged
  (box-filtered) kernel; FCD addresses the kernel/aliasing part, the voxel indicator stays a
  staircase.
- Evans K.F., Stephens G.L. (1991), Lorentz-Lorenz averaging of boundary susceptibility: known
  only via Yurkin & Hoekstra 2007 review (not read).
- Yurkin M.A., Min M., Hoekstra A.G. (2010). Application of the DDA to very large refractive
  indices: Filtered coupled dipoles revived. Phys. Rev. E 82, 036703.
  doi:10.1103/PhysRevE.82.036703. READ: ABS. Compares FCD, lattice dispersion relation (LDR) and
  radiative reaction correction; FCD good for all m tested, LDR bad for m = 10+10i.
  (Description of FCD as a band-limited, translation-invariant lattice Green tensor that keeps the
  FFT is OWN BACKGROUND KNOWLEDGE, not from the abstract; to verify before use.)

### F5 Geophysics (prisms, staircase topography)

Finding: geophysical FFT/BTTB prism modelling uses the same structure as the Newell FFT
(analytic prism kernel, Toeplitz, FFT), but I found NO staircase-correction technique for the
SOURCE body in that literature; "undulating surface" there means observation points.
- Hogue J.D., Renaut R.A., Vatankhah S. (2020). A tutorial and open source software for the
  efficient evaluation of gravity and magnetic kernels. Computers & Geosciences 144, 104575.
  doi:10.1016/j.cageo.2020.104575, arXiv:1912.06976.
  READ: abstract + section scan. BTTB kernels of uniform prisms, FFT forward; padding; no
  partial-cell treatment found.
- Chen L., Liu L. (2019). Fast and accurate forward modelling of gravity field using prismatic
  grids. Geophys. J. Int. 216(2), 1062-1071. doi:10.1093/gji/ggy480. READ: ABS. Analytic prism
  kernel + 2D discrete convolution.
- Wu L., Tian G. (2014). High-precision Fourier forward modeling of potential fields.
  Geophysics 79(5), G59-G68. doi:10.1190/geo2014-0039.1. READ: ABS. Gauss-FFT: shift-sampling with
  Gaussian-quadrature shifts removes aliasing / imposed periodicity / edge effect of
  Fourier-domain forward modelling. Not a staircase method (concerns the spectral kernel).
- Parker R.L. (1973). The rapid calculation of potential anomalies. Geophys. J. R. Astron. Soc.
  (Geophys. J. Int.) 31(4), 447-455. doi:10.1111/j.1365-246X.1973.tb06513.x. READ: ABS. Series of
  FFTs of powers of the interface height for an uneven layer. Transfer (own inference): only for
  a layer z = h(x, y), not for a closed sphere; no direct use.
- Polyhedron formulas (Guptasarma & Singh 1999; Fabbri 2008) listed under F1.

### F6 Exchange-side corrections at curved boundaries

Known: Garcia-Cervera, Gimbutas & E 2003 (exchange and stray field in incomplete cells); Parker,
Cerjan & Hewett 2000 (embedded curve boundary). No new micromagnetic source found.
General cut-cell Poisson/Neumann methods (metadata only, NOT read; not cited as evidence):
- Johansen H., Colella P. (1998). A Cartesian grid embedded boundary method for Poisson's equation
  on irregular domains. J. Comput. Phys. 147, 60-85. doi:10.1006/jcph.1998.5965.
- Papac J., Gibou F., Ratsch C. (2010). Efficient symmetric discretization for the Poisson, heat
  and Stefan-type problems with Robin boundary conditions. J. Comput. Phys. 229, 875-889.
  doi:10.1016/j.jcp.2009.10.017.
- Donahue M.J., Porter D.G. (2004). Exchange energy formulations for 3D micromagnetics. Physica B
  343, 177-183. doi:10.1016/j.physb.2003.08.090 (metadata only).
Transfer (own inference): the exchange operator is local, so a cut-cell (volume fraction + face
apertures) Laplacian with Neumann on the true surface does not touch the FFT at all.

### F7 Averaging strategies (rotation, offsets, Richardson)

- Yurkin et al. 2006 II (see F4): main evidence. Shape error is quasi-random in the cell size
  because it depends on where the surface cuts the boundary cells; Richardson extrapolation
  for spheres is reliable only at fine grids (ymin < 0.15 at m = 1.5) and needs 9 grid sizes
  over a factor 4. -> Richardson in dx WITHOUT geometry locking is weak for staircase error.
- Yurkin & Hoekstra (2007). The discrete dipole approximation: an overview and recent
  developments. JQSRT 106, 558-589. doi:10.1016/j.jqsrt.2007.01.034, arXiv:0704.0038.
  READ: FULL (Sect. 3 shape-error passages). Shape errors concentrate in the internal fields near
  the boundary (Druger & Bronk: average internal-field errors > 30 % for a sphere, m = 1.8,
  y = 0.17); boundary-susceptibility averaging (Evans & Stephens 1991, Lorentz-Lorenz mixing with
  the volume fraction f) and WD (Piller) as remedies; adaptive dipole sizes (Bourrely et al.)
  give about 2 times better accuracy but break the regular lattice.
- Adams M.P., Michels A. (2026). Minimal model for vortex nucleation and reversal in spherical
  magnetic nanoparticles. Phys. Rev. B 113 (22). doi:10.1103/8rkf-j2kn, arXiv:2601.17176.
  READ: FULL (text scan). Context only: remaining difference between their model and mumax3
  for a 40 nm sphere (2 nm cells) "is highly sensitive to the discretization cell size",
  attributed to FD discretisation effects (cites Holt et al.).
- No source found for: (a) ensemble averaging over sub-cell offsets or grid orientations of a
  particle as an error-reduction method in micromagnetics or DDA; (b) Richardson extrapolation
  at a locked (fixed) staircase geometry. One web query each; documented as a gap.
  (Garcia-Cervera 2003 [known] reports that the switching mode depends on grid orientation.)

### F8 Other (form factors, filtered kernels)

FFT-based computational homogenisation (Moulinec-Suquet) is the closest structural analogue:
regular voxel grid, FULL periodic BC, FFT Green operator, staircase interfaces.

- Kabel M., Merkert D., Schneider M. (2015). Use of composite voxels in FFT-based homogenization.
  Comput. Methods Appl. Mech. Eng. 294, 168-188. doi:10.1016/j.cma.2015.06.003. READ: NOT
  DIRECTLY (no abstract in OpenAlex; content known only via the two papers below).
- Lendvai J., Schneider M. (2024). Level-set based composite voxels for FFT-based computational
  homogenization of thermal conductivity. ECCOMAS Congress 2024, Lisbon. PDF (scipedia draft):
  https://scipedia.com/wd/images/b/b0/Draft_Sanchez_Pinedo_909610227pap_645.pdf.
  READ: FULL. Laminate composite voxel: cut voxel gets the conductivity of a planar laminate,
  harmonic mean along the interface normal n, arithmetic mean tangential (their Eq. 20).
  Volume fraction and n from a linear level-set fit to the edge intersection points; volume by
  Mirtich formulas. Periodic cell with ONE SPHERE (17 % volume, contrast 10), Willot
  discretisation, CG solver, reference 1024^3: relative error of effective flux 2 % -> 0.1 %
  (16 ... 256 voxels per edge) without, 0.2 % -> below 0.01 % with composite voxels, i.e. about
  one order of magnitude. Fibre composite: 3 % -> 0.5 % without, 1 % -> below 0.02 % with; "about
  one order of magnitude which roughly amounts to using only one fourth of the resolution per
  dimension"; about 20 % more computing effort.
  Local flux: voxelisation artefacts at the interface "removed to a large extent".
  Transfer (own inference): the magnetostatic analogue is a linear-response problem; in
  micromagnetics M is the state variable, so the laminate rule maps to (a) the cell weight f
  (volume fraction) and (b) a corrected self/near term for the inclined interface charge. It is
  not a drop-in rule for LLG.
- Keshav S., Fritzen F., Kabel M. (2022). FFT-based homogenization at finite strains using
  composite boxels (ComBo). Comput. Mech. 71(1), 191-212 (2022 online), doi:10.1007/s00466-022-02232-4,
  arXiv:2204.13624. READ: FULL (abstract, Sect. 1). Image-based normal detection reduces the error
  of phase-averaged stresses by about 30 % vs Kabel 2015; volume fractions usually limited to
  5 ... 95 % for robustness (Voigt average outside).
- Schneider M. (2025). Accurate and consistent composite voxel methods for digital images in
  computational micromechanics. Comput. Mech. doi:10.1007/s00466-025-02716-z. READ: SNIP only.
- Lendvai J., Schneider M. (2024). Assumed strain methods in micromechanics, laminate composite
  voxels and level sets. Int. J. Numer. Methods Eng. doi:10.1002/nme.7459. READ: title only
  (OpenAlex). Cited in the ECCOMAS paper as the justification of laminate voxels by a kinematic
  assumption.
- Farjadpour A., Roundy D., Rodriguez A., Ibanescu M., Bermel P., Joannopoulos J.D., Johnson S.G.,
  Burr G.W. (2006). Improving accuracy by subpixel smoothing in the finite-difference time
  domain. Opt. Lett. 31(20), 2972. doi:10.1364/OL.31.002972. READ: ABS. Smoothing of the
  dielectric function in boundary pixels designed from a perturbation-theory criterion; "the
  only one that attains quadratic convergence with resolution for arbitrarily sloped
  interfaces"; sharp corners remain difficult. (Different solver: FDTD.)
- Kottke C., Farjadpour A., Johnson S.G. (2008). Perturbation theory for anisotropic dielectric
  interfaces, and application to subpixel smoothing of discretized numerical methods. Phys. Rev.
  E 77, 036611. doi:10.1103/PhysRevE.77.036611, arXiv:0708.1031. READ: ABS. First-order shift of an
  interface = surface integral of continuous field components; gives the smoothing rule.
  Transfer (own inference): same principle as the laminate composite voxel: choose the
  boundary-cell effective parameters so that the FIRST-ORDER error of the interface position
  vanishes. In micromagnetics the analogue for the stray field is a first-order (local
  anisotropy-like) cell term, i.e. the D&M local correction.
- Hamm J., Renn F., Hess O. (2014). Dispersive media subcell averaging in the FDTD method using
  corrective surface currents. IEEE Trans. Antennas Propag. 62(2), 832-838.
  doi:10.1109/TAP.2013.2292073, arXiv:1303.1403. READ: ABS. Volume averaging + corrective
  electric/magnetic SURFACE CURRENTS at interfaces; Mie-scatterer test; "approximately restore
  second order accuracy". Analogue of a surface-charge correction (F3), different solver.
- Hamada S., Kobayashi T. (2007). Analysis of electric field distribution induced by 50 Hz
  magnetic fields utilizing fast-multipole surface-charge-simulation method for voxel models.
  ICCM 2007, Hiroshima. PDF: Kyoto University repository. READ: FULL (text). Voxel-face surface
  charges; fields near the boundary show "sharp fluctuations" from the staircase; remedy is only
  a post-processing "near boundary smoothing" (weighted average over 3x3x3 neighbours).
  Not a physical correction.
- Dawson T.W. (2001). Evaluation of modeling accuracy of power frequency field interaction with
  the human body. ACES J. 16(2), 162-172. READ: ABS (PDF is scanned images, not read). Staircase
  errors identified with analytic sphere tests; reduction methods (not read). Solver: SPFD.
- Durhuus F.L., Insinga A.R., Bjork R. (2026). Including sample shape in micromagnetics with 3D
  periodic boundary conditions. arXiv:2604.08777. READ: FULL (abstract, Sect. I-II start).
  Not staircase: proves that for large samples only the average magnetisation carries the
  sample-shape (outer surface) effect; adds a shape term to 3D PBC. Relevant only to our outer
  surface compensation, which already uses this idea.

## Citer lists checked

- Donahue & McMichael 2007 (44 citers, OpenAlex): only applications, reviews and solver papers;
  no follow-up method that extends the correction (Holt et al. 2025 TMAG is the known
  discretisation-anisotropy paper; Wei, Wang, Ding, Gao 2009 TMAG 45, 3035 = triangular-prism
  cells, known). Two sphere applications (Djuhana et al. 2019, 2021) with no correction study
  in the abstract.
- Lepadatu 2015 roughness (21): applications, Boris papers, multilayer convolution 2019.
- Garcia-Cervera, Gimbutas, E 2003 (64): AMR 2006, NUFFT 2008, "Moving magnets in a
  micromagnetic FD framework" (PRE 2018, not opened), theory papers; no new staircase method.
- Chernyshenko & Fangohr (11; arXiv and JMMM records merged in OpenAlex): Bjork & d'Aquino 2023,
  Liang & Yan 2024, MagTetris+ 2026; no staircase method.
- Holt et al. (TMAG 2025: 2 citers; arXiv: 0): Adams & Michels 2026 (sphere context, see F7).
- Parker, Cerjan, Hewett 2000 (10): reviews and D&M 2007 only.
- Lebecki, Donahue, Gutowski 2008 (95): PBC applications; Bruckner 2021 true PBC; Ai & Lomakin
  2024 CPC; Durhuus 2026 not listed (too new) but found by search.
- Yurkin et al. 2006 II (53): Zhu, Liu, Yurkin 2023 EMA; Smunev et al. 2015 rectangular dipoles;
  Zhu et al. 2019 MDR of spheres (ABS only); Yurkin et al. 2010 FCD.

## Own estimates and calculations (marked; not from sources)

- Staircase surface area: for a cell-centre voxelised sphere the staircase face area tends to
  (|n_x| + |n_y| + |n_z|) integrated over the surface; the mean of |n_i| over a sphere is 1/2,
  so the staircase area tends to 3/2 of the true area for any dx. (Own calculation.) This bias
  does not average out over grid offsets or rotations.
- Number of cut (boundary) cells per sphere: about 6 pi (R/dx)^2 (same argument). D = 100 l_ex,
  dx = 3 l_ex: about 5e3 per sphere; D = 300 l_ex, dx = 3 l_ex: about 4.7e4; D = 300 l_ex,
  dx = 1.5 l_ex: about 1.9e5. Four spheres: 2e4 ... 7.5e5. (Estimate.)
- Sparse near-field correction storage, radius 1 cell (27 neighbours, 9 entries each): about
  250 floats per boundary cell -> 4e4 boundary cells: 1e7 floats (80 MB in float64); radius 2
  (125 neighbours): about 1100 floats per boundary cell -> 2e8 floats (1.7 GB) for 1.9e5 cells.
  (Estimate.)

## Draft table rows (for the report)

| Family | Best sources | Corrects | Validation | Order | Newell FFT + full PBC | Cost (estimate) | Main risk for us |
|---|---|---|---|---|---|---|---|
| F1 near-field correction (cut-cell exact tensors) | D&M 2007 (known) + Bonn slides; Piller 1997 WD; Yurkin I; Nielsen 2019, Fabbri 2008 (exact polyhedron fields) | stray near field (self + neighbours) | D&M: thin films (oblate spheroid, square; FMR); WD: 3D spheres (optics) | WD: still 1st order, smaller constant (Yurkin I); cell-local thickness correction only: demag RMS error 118 % -> 29 % (thin spheroid, D&M slides) | yes, sparse local term; PBC: free-space kernel + minimum-image wrap (precedent Ai 2024) | days to weeks of code; 0.1-2 GB matrix | D&M 'pushed' local-anisotropy form needs m nearly uniform over the near zone (their own condition) -> full K_ij needed; symmetry of K for an energy (inference) |
| F2 refined grids keeping FFT (pFFT/BAIM, NUFFT, AMR, multilayer) | Ai/Duan/Lomakin 2024; Phillips & White 1997; Exl & Schrefl 2014; Garcia-Cervera & Roma 2006; Lepadatu 2019 | geometry of sources (fine surface or fine patches) | FEM cube unit cell (Ai); free-space FEM (Exl) | controlled error level (1e-3 used) | yes for pFFT (periodic kernel + G0 near correction); AMR only with heavy changes; multilayer: one axis only | weeks; new solver path | effectively a second solver; projection/interpolation error near deep steps |
| F3 charge-based (smooth surface charge minus staircase charge) | none direct; Exl & Schrefl 2014; Hamm 2014 (FDTD surface currents); Mayo-type (Poisson) | stray near field | none for FD micromagnetics | Hamm: about 2nd order (FDTD) | yes as sparse local term (same as F1 in another basis) | similar to F1 | needs M on the true surface: extrapolation over 1 cell when M varies on the cell scale |
| F4 DDA boundary corrections (WD, EMA, FCD) | Piller 1997; Zhu 2023; Piller & Martin 1998; Yurkin I/II | shape error of boundary cells; FCD: kernel aliasing | 3D spheres vs Mie, aggregates | WD 1st order; integral quantities x24 (WD), up to x60 (EMA), x100 (FCD cross-sections) | yes: cell weights / self-terms / filtered translation-invariant kernel | low (cell-local) | linear-response rules; DDA metrics are integral quantities, local surface fields may gain less |
| F5 geophysics prisms | Chen & Liu 2019; Hogue 2020; Wu & Tian 2014; Parker 1973; Guptasarma & Singh 1999 | (no staircase correction found); exact polyhedron formulas | gravity/magnetic bodies | - | structure identical (prism kernel + FFT); polyhedron formulas usable inside F1 | - | none beyond F1 |
| F6 exchange at curved boundary | Garcia-Cervera 2003, Parker 2000 (known); Johansen & Colella 1998, Papac 2010 (not read) | exchange, Neumann on true surface | known (D&M/GC tests); CFD papers not read | not read | yes (local, no FFT) | days | small cut-cell fractions -> stiff time step (inference) |
| F7 averaging / extrapolation | Yurkin II (Richardson with shape error); no source for offset ensembles or locked-geometry Richardson | random part of the staircase error; error bar | DDA spheres | quadratic fit; reliable only at fine grids | yes, no solver change | N runs (N = 4-9) | systematic part (3/2 area, own calculation) stays; locked-geometry Richardson converges to the staircase body |
| F8 composite voxels / subpixel smoothing / shape amplitude | Lendvai & Schneider 2024; Kabel 2015; Farjadpour 2006; Kottke 2008; Tandon/Beleggia 2004 | effective boundary-cell properties; uniform-M body field | periodic FFT cell with one sphere (Lendvai); FDTD | sphere: about 10x smaller error (fibre case: = 4x coarser grid); FDTD 2nd order | yes (cell-local; uniform-mode term per particle) | low-medium | rules are for linear media; for LLG they reduce to f weights + D&M-type local term (inference) |

## Gaps of the search

- No source tests any correction on a hysteresis or loss observable of spheres; validations are
  field errors, energies, FMR frequencies or scattering cross-sections.
- No micromagnetic paper combines a staircase correction with FFT + full 3D PBC; the PBC
  precedent for a free-space near-field correction comes from FEM (Ai, Duan, Lomakin 2024).
- Full text not reached: Zhu et al. 2023 (EMA), Kabel et al. 2015, Piller 1997, Piller & Martin
  1998, Tandon et al. 2004, Johansen & Colella 1998, Papac et al. 2010, Dawson 2001 (scanned).
- No source for offset/rotation ensembles or for Richardson extrapolation at a locked geometry.
- Not opened: "Moving magnets in a micromagnetic finite-difference framework" (PRE 2018),
  Hillan & Foss ASEG 2012 (self-demagnetisation corrections, page timed out).
- Semantic Scholar API was rate-limited; citation lists come from OpenAlex only.
