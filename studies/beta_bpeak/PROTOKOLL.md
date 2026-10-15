# Protokoll: β(B_peak) study with magnum.np (fp32, reduced units)

Owner: Christoph Vogler (Chris). Branch: `study/beta-bpeak-pinning` (based on
`perf/memory-optimizations`, the fp32 fork). Folder: `studies/beta_bpeak/`.

Rules: documents in ASD-STE100 Simplified Technical English, SI units only.
**Chris makes the decisions.** Before you change the model, the parameter
matrix or the budget, give him the options with pros, cons, chances and risks,
and wait. Engineering fixes (bugs, logging, speed) need no decision, but tell him.

## 0. Current state (2026-10-12)

Only the latest valid results. Sources: 7.65 (main result, conditioning
test), 7.66 and audit_20261012/REPORT.md (independent audit), 7.63 / 7.64
(mesh), PLAN 4.6 (running).

**Terms.** l_ex = exchange length (3.34 nm at A = 10 pJ/m); K_d = μ0 Ms² / 2
= 895 kJ/m³; w = loop work in units of K_d per volume of the powder core
(voids included); P = f · K_d · w = loss density in W/cm³ per core volume
at f = 30 MHz. L_eff = edge of the
Herzer cell (random-anisotropy cell), K_eff = K_d (l_ex / L_eff)². K_σ =
r_p · K_eff = stress anisotropy per particle. Simulation model = the micromagnetic computation (LLG with exchange,
demagnetisation, Herzer and stress anisotropy, drive field; T = 0; periodic box of 4 spheres). Main fit = the statistical
fit (Gamma GLM) of P over the design parameters of all α 0.1 design runs (before 2026-10-14 also called "design
model"). Range effect = change of an
output from the low to the high end of one factor, the other factors
averaged over the design points. Lever = a range effect that is significant
after the Holm correction over the 8 main tests (5 %). Event cycle: w_loop ≤
0 or |w_dis / w_loop − 1| > 0.20.

**Basis.** Fictitious nanocrystalline powder, d ≈ 1 µm, α = 0.1, B_peak
9 … 150 mT, mesh dx = 3 l_ex (production). Converged virgin relaxation;
steady rule: a stage keeps the measured cycles after its last event cycle,
at least 2 (7.54); finished reruns replace their original. Data: 116
design runs (96 Sobol, 8 centre, 12 corner), mesh pairs dx 3 / dx 1.5 at
the base point (d/l_ex 300, n 4) and at K3 (d/l_ex 212, n 6), 18
conditioning runs. The audit checked 275 run folders (256 current, 19
replaced by reruns): complete, converged, consistent; ≥ 2 calm cycles in
every stage except the 8 runs at α 0.01. Errors: mean ± SE over independent
seeds; CI = confidence interval.

**1. β of the loss curve (power-law fit per run, dx 3, n 116).**
- 9 … 150 mT: β = 2.12 ± 0.02 (sd 0.25, range 1.59 … 2.60); β > 2 in 77
  of 116 runs (fit over the nominal B, as the design analysis; 78 with the
  realised B).
- Local β (3-amplitude windows): 2.15 at 9 … 70 mT, 2.07 at 50 … 100 mT,
  1.98 at 70 … 150 mT; β(70 … 150) − β(9 … 35) = −0.17 ± 0.04.
- One power law gives the trend, not the detail (rms 0.13 in ln P against
  ≈ 1 % cycle error). There is no common break point.

**2. Mechanism.**
- The response of the powder is linear: μ_r = 6.9 (quartiles 5.9 … 8.2),
  constant in B to 0.1 %, set by φ; harmonics < 10⁻³ of the fundamental.
- Thus β − 2 = d ln(sin δ) / d ln B exactly (δ = loss angle, 0.33 …
  0.48°; confirms 7.22 / 7.24). Every lever on β acts through the growth
  of the loss angle with B.
- Where the loss is (dx 3, calm accumulator windows, 7.69): the 10 % of
  the cells with the largest loss carry 69 … 76 % of P (each run ≥ 50 %).
  The particle centre (r/R < 0.5, 12.5 % of the volume) carries 41 % of P
  at 9 mT and 18 % at 150 mT. The particle rim (r/R ≥ 0.8, 49 % of the
  volume) carries 21 % at 9 mT and 50 % at 150 mT. Thus the loss stays in
  a few zones at every B, and with rising B the weight moves from the
  centre to the rim. The zones lie mostly where 𝐦 is non-uniform (61 … 69 %
  of P up to 70 mT, 55 % at 150 mT). The particles do not touch (gap ≥ 1.7
  cells).
- Loss channels (dx 3, α 0.1, 7.70, 7.71): the share of the cell loss
  from the local motion at f falls from 0.88 at 9 mT to 0.40 at 150 mT.
  Motion at other frequencies takes over (0.08 → 0.52). Jumps carry ≤ 12 %
  on average (fine mesh ≤ 7 %). Exactly β = β_f + shift: the fall of the
  share at f adds 0.25 … 0.33 to β, and the A effect on β (+0.41 ± 0.09)
  sits mostly in β_f, the growth of the loss at f (+0.33 ± 0.07). At
  α 0.01 the share at f is only ≈ 0.1 at all B: the split depends strongly
  on the damping. The nature of the other frequencies (harmonics, local
  precession, slow drift) is open; option O2 (7.71) measures it.
- The switched volume per half cycle grows about ∝ B (0.65 % at 9 mT,
  9.8 % at 150 mT) and follows β (r = 0.91; partly by construction, 7.46).
- Events (energy release in the first measured cycles) are the end of a
  transient, not a slow drift (7.55, 7.65).

**3. Levers on β, 9 … 150 mT (dx 3, main result, range effects with 95 %
CI, results/sobol_analysis_steady2_mc2.md).**

| Factor (low → high) | Δβ | Class |
|---|---|---|
| A 10 → 20 pJ/m at d = 1 µm (= d 1.0 → 0.71 µm at A 10 pJ/m) | +0.57 (+0.48 … +0.67) | lever |
| K_eff 1.0 → 6.2 kJ/m³ (L_eff 30 → 12 l_ex) | +0.15 (+0.07 … +0.22) | lever |
| stress, r_p 0 → 3 | +0.14 (+0.02 … +0.26) | at the decision limit |
| packing φ 0.55 → 0.69 | −0.06 (−0.16 … +0.03) | at the decision limit (not a lever without the reruns, 7.65) |
| damping α 0.1 → 0.01 (centre, n 8) | sign depends on the rule | not determined |

- A and K_eff interact strongly (audit bootstrap, sobol_extra.md). A
  effect: +0.93 (+0.80 … +1.05) at K_eff 1.0 kJ/m³, +0.55 at 2.5, +0.16
  (n.s.) at 6.2 kJ/m³ (deck bootstrap: +0.93 (+0.80 … +1.08), +0.55
  (+0.44 … +0.64), +0.16 (−0.02 … +0.33)). K_eff effect: +0.53 (+0.38 … +0.66) at A 10 pJ/m,
  −0.24 (−0.40 … −0.09) at A 20 pJ/m. Lowest β ≈ 1.6 at low A and low
  K_eff; highest ≈ 2.5 at high A and low K_eff. Reading (hypothesis): weak
  anisotropy lets exchange and stray field set the state, so the size in
  units of l_ex acts strongly.
- Real alloy (Herzer coupling K_eff ∝ A⁻³): A × 2 gives Δβ +0.40
  (+0.28 … +0.53); an extrapolation (the design covers K_eff × 6.25, the
  coupling needs × 8).
- By range (secondary outputs, hints): A acts most at 35 … 70 mT (+1.00),
  K_eff at 70 … 150 mT (+0.31); stress at 35 … 70 mT (+0.29, p 0.035
  without a multiplicity correction).

**4. Effects on P (dx 3).** Lever (primary output): A at 50 mT. Effects:
A 10 → 20 pJ/m × 0.33 at 10 mT, × 0.74 at 50 mT, × 1.21 at 100 mT, × 1.32
at 150 mT (crossing at ≈ 70 mT on dx 3); K_eff 1.0 → 6.2 kJ/m³ × 0.73 at
10 mT, × 1.17 at 150 mT; φ 0.55 → 0.69 × 1.34 at 10 mT (× 1.07 without the
dilution).

**4a. A low β does not mean a low loss (7.68, coarse mesh).** The best
material depends on B_max: at 9 mT small particles with low K_eff (β 2.55,
2.0 W/cm³, × 6 below the worst), at 150 mT large particles with low K_eff
(β 1.61, 1170 W/cm³). The switch point is ≈ 70 … 100 mT on the coarse mesh
and below 50 mT for the only fine-mesh pair. β alone does not rank
materials. The loss densities in W/cm³ of all levers, the loss location,
the turn angles and the real-alloy path are in 7.68.

**5. Mesh: dx 3 against dx 1.5 (mesh error = dx 3 − dx 1.5, same seed).**

| Output | dx 1.5 needed? | Error to assume with dx 3 |
|---|---|---|
| P level | yes | ± 0.5 in ln P (× 0.6 … 1.6) at 9 … 70 mT, ± 0.25 at 100 mT (estimate). The sign differs between the 2 tested points (base × 1.48 / 1.40 / 1.25, K3 × 0.63 / 0.70 / 0.87 at 9 / 50 / 100 mT); cause open (PLAN 4.6) |
| β, 9 / 50 / 100 mT | no, with an error bar | base −0.06 ± 0.05, K3 +0.12 ± 0.03; K3 − base +0.18 ± 0.06 (p 0.039). Assume ± 0.2 (estimate). 100 … 150 mT not checked |
| β, 9 … 35 mT and 9 → 50 mT | no | K3 9 … 35 mT −0.01 ± 0.11 (not resolved); 9 → 50 mT base −0.02 ± 0.03, K3 +0.05 ± 0.05; K3 20 … 50 mT −0.16 ± 0.10 |
| β above 50 mT | yes | K3 +0.35 ± 0.11 (50 … 100 mT); base 50 → 100 mT −0.17 ± 0.15 (other sign, n.s.) |
| A effect on β | partly | fine mesh ≈ +0.38 ± 0.11 instead of +0.57 (estimate from one pair that also changes L_eff, r_p, φ) |
| A effect on P at 50 mT | yes (hint) | × 0.79 (dx 3) against × 1.59 (dx 1.5, p 0.034, one uncorrected test, confounded pair) |
| K_eff and stress effects | not known | no fine-mesh check (K2 not run); they come from > 35 mT, where the β mesh error reaches ± 0.35 |

- At K3 (1 point, n 6) the dx 3 grid cannot hold the states of dx 1.5
  (7.63). dx 1.5 itself is not shown to be converged.
- Pending: the mesh check at d/l_ex 212 for the base microstructure
  (PLAN 4.6, running, end ≈ 2026-10-12 23 UTC, estimate). It shows whether
  the sign of the level error follows the size or the microstructure.

**6. Convergence.** 30-cycle test at the centre point only (n 6): the
design values stand at 9 / 20 / 35 mT (drift of P ≤ 3 %, 90 % CI; Δβ
−0.003 ± 0.018 and −0.004 ± 0.021). Reruns are bit-identical to their
originals up to the first extended stage.

**7. Damping (centre; rule since 2026-10-12: steady loops ≥ 2, decision
Chris, 7.66).** α 0.01 gives × 0.62 … 0.81 of P at 35 … 150 mT (n 4 … 6
seeds with ≥ 2 calm cycles at both α). At 9 mT no seed has 2 calm cycles
at α 0.01; at 20 mT only 3 of 8, with a too wide interval (× 0.80, SE 0.61
in ln P). 30-cycle runs at 9 and 20 mT follow (PLAN 4.7). Most of the loss
does not scale with α (hypothesis: irreversible jumps set it). The effect on β is not determined.

**8. Limits.** One frequency (30 MHz); d/l_ex 212 … 300 only; P levels from
dx 3 carry the error of item 5; the mesh share of the K_eff and stress
effects is unknown. Hint (REPORT D.4.6): for r_p > 0.5, particles with the
stress axis at 37° to the drive have × 1.45 … 1.84 the loss per volume of
those at 79° (n 100 runs); no further growth with r_p above ≈ 0.3; at r_p
0.02 … 0.26 the ratio is already × 1.27 … 1.40 (n 9). No r_p = 0 reference
run exists.

**9. Open (decisions of Chris; cost estimates).** Stress step test r_p 0 /
0.1 / 0.3 / 1 (≈ 3 $); fine-mesh pair at L_eff 30 (≈ 23 $); K2 mesh check
(≈ 25 $); α 0.01 with ≥ 30 cycles at 9 … 20 mT (queued, ≈ 2.5 $); dx 1.5
convergence with one dx 1 seed at d/l_ex 212 (≈ 100 V100-h ≈ 10 $). Deck
rebuild after the mesh check at d/l_ex 212.

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
     K_eff/√N ≈ 69 J/m³ (N ≈ 8200 full cubes; 9400 with the partial surface
     cubes), acts only on a uniform rotation of
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
   K_p = r_p K_eff per particle. The easy axes are deterministic: the 4
   body diagonals of a cube (isotropic set, ⟨u uᵀ⟩ = I/3), turned so that
   along the drive ⟨cos²θ⟩ = 1/3 and ⟨cos⁴θ⟩ = 1/5 (cos θ = 0.795 / 0.188 /
   0.188 / 0.795; since 2026-10-07, 7.37; the earlier set was not isotropic,
   7.36). There is no orientation statistics.

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
  (Updated 2026-10-04: the target is the total loss at 30 MHz, Chris
  2026-10-03; the damping loss is part of it and α is a model parameter
  that changes β, 7.17. It is not treated as an artifact.)
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
| Q_eff window | 1.1·10⁻³ … 6.9·10⁻³ | L_eff ≥ 4 cells (resolved); d/L_eff ≥ 10 (walls in the particle, ≥ 500 cubes per particle, small realisation scatter); superseded by D3 (PLAN 4.0): d/L_eff ≥ 7, residual check vs d/L_eff |
| Mesh-identical cubes | box edge a = multiple of 12 l_ex, L_eff = multiple of 6 l_ex | dx = 3, 2 and 1.5 see exactly the same cubes (tested on a common fine grid, with a negative control); d is adjusted by ≤ 2.5 % to keep φ exact |
| Contacts | at least one void cell between particles (no face, edge or corner contact) | no exchange coupling, no extreme stray fields at contacts; at d/l_ex = 150 and dx = 3, φ = 0.69 and 0.71 are not possible |
| Size | d/l_ex = 300 (1 µm), fixed; superseded by D3 (PLAN 4.0): d/l_ex 212 … 300 as a factor (at d = 1 µm: A 10 … 20 pJ/m) | Chris |
| Amplitudes | B_peak = 150 / 100 / 70 / 50 / 35 / 20 / 9 mT (`--b_list`; 9 mT keeps 10 mT inside the range). The drive of each stage comes from b/h of the stage before and is corrected once after cycle 0 with b/h of cycle 0 (factor limited to 0.67 … 1.5) | study range 10 … 150 mT; no time on saturation. b of cycle 0 is within ±4 % of the later cycles (V100 data), thus each stage hits its target to ≈ 4 % |
| Drive | sinusoidal H (not controlled B) | below 150 mT the core is almost linear: B(t) has 0.2 … 3 % harmonics, b varies by 0.1 … 3 % from cycle to cycle (mesh test data). Thus sinusoidal H ≈ sinusoidal B (Steinmetz condition) |
| Initial state | **decided: A, virgin state relaxed to convergence with a hard gate (7.13, 7.15).** Original text: open: decided by the initproto test (section 7.4). Candidates: A = virgin state (random m, full relaxation at H = 0 and α = 1: the T = 0 analogue of the anneal and the cooling without field; a real powder is never magnetized before use), B = AC demagnetization (decaying saturating cycles, as in IEC 60404-6), C = 1 saturating cycle (present) | C gave minor loops around a remanent state (offset up to 0.55 at 23 mT) |
| Protocol | A and B: ascending amplitudes (the classic Rayleigh procedure: every larger loop erases the smaller ones, the loops stay centred). C: descending. 7 cycles each (fixed), cycles 0 and 1 discarded (drive correction) | AC demagnetization. At T = 0 the soft powder does not lock into a periodic cycle (steady test): the mean over cycles is the measurement, not a single steady cycle |
| Loss per cycle | w = mean of w_loop (∮H dB) over the kept cycles; w_dis (LLG dissipation) as a check (`--w dis`) | w_loop of one cycle contains the change of the stored energy when the cycle is not closed; this part averages out over cycles. w_dis scatters less |
| Damping | α = 0.1 (main design); α 0.01 only at the 8 centre replicates (PLAN 4.0, decision 2026-10-06) | α is a model parameter: α 0.01 lowers the loss at 50 mT by ×0.59 (paired n = 4, p = 0.002, 7.31); an α effect on β is not established (paired −0.09 ± 0.17, 7.31; the earlier −0.27 ± 0.08 of 7.17 rests on spliced runs) |
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

- **Withdrawn (7.11, 7.13): all d/l_ex 300 values in this table come from
  runs with an unconverged virgin state. Mesh decision: see 7.28 (dx 3
  necessary and sufficient for β at 1 µm; dx ≤ 1.5 necessary for the
  absolute loss, sufficiency unknown).** Original text: At d/l_ex = 300 (1 µm) dx = 3 is sufficient for β (−0.20 ± 0.23). The
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
- **Validity (withdrawn, see above: open until mesh300):** absolute losses at 1 µm are known to ≈ 30 % (mesh). β and
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
  r_p = 3: Δβ −0.31 ± 0.22, loss ×0.24 (old, non-isotropic stress axes,
  7.36: r_p part not valid). L_eff = 30 at d/l_ex = 96 has
  d/L_eff = 3.2 (< 5, outside the valid window): measure this effect at 1 µm.
- β ≈ 2 between 9 and 50 mT. The dynamic (damping) share at 30 MHz and
  α = 0.1 is not measured; it also gives β = 2. Check with f vs 2 f first.

### 7.9 Frequency test (PLAN step 1, 2026-10-02/03): NEGATIVE

d/l_ex = 300, dx = 3, α = 0.1, virgin state, 7 amplitudes; f (30 MHz) and
2 f, 4 seeds each; base point (L_eff 12, r_p 0) and corner (L_eff 30, r_p 3;
the corner used the old, non-isotropic stress axes, 7.36).
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
  No effect on β shown here (n = 3, 9 … 20 mT only); **superseded by 7.17:
  with all valid runs α changes β.**
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
| P at 10 / 50 / 100 mT [W/cm³] (ln w interpolated, analyze.features) | 9.86 / 192 / 576 | scatter factor ×/ 1.51 / 1.25 / 1.14 | ×/ 1.16 / 1.08 / 1.05 | |

- Energy balance w_dis/w_loop outside 0.9 … 1.1 in 6 of 56 stages:
  s921 at 9 mT (1.40) and 100 mT (1.83), s922 at 20 mT (2.07), s925 at
  9 mT (1.13) and 35 mT (1.21), s926 at 35 mT (1.27); events with energy
  release also in kept or skipped cycles of s923 (20 mT, cycle 1: 6.6):
  events in at least 5 of 8 runs. They coincide with drops of the number of
  large-angle pairs (e.g. 28 → 23, 33 → 20, 29 → 22): rearrangements of the
  structure. (Corrected 2026-10-06, 7.27: they are equally frequent at
  dx 1.5, thus not a coarse-mesh artifact.)
  s921 at 9 mT: kept cycle 2 has balance 4.0 and w ≈ 3× below the other
  seeds (P at 10 mT 3.8 vs 9.4 … 13.8 W/cm³); this one stage gives its
  β_pow 2.11. s922 at 20 mT is its second stage: an irreversible event
  (balance 2.9 / 5.9 / 2.8 in cycles 1 … 3, closure 0.026 > 0.02), not a
  release of the virgin state.
- By PLAN 4.3a the flagged runs stay in; sensitivity: without s921 the
  scatter of β (power law) falls from 0.15 to 0.04. The realisation scatter
  is dominated by events, not by the cube landscape.
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
  β power law −0.34 ± 0.14 (Welch df ≈ 2.6, p = 0.10: not significant
  alone); pooled with the spliced runs (6 vs 11): −0.27 ± 0.08 (p = 0.012).
  3-point windows on 9/20/50/100/150 mT (results/baseline_eval.txt): clean
  −0.43 ± 0.17 / −0.40 ± 0.13 / −0.15 ± 0.11, pooled −0.32 ± 0.09 /
  −0.32 ± 0.10 / −0.17 ± 0.10 (windows 9-20-50, 20-50-100, 50-100-150).
  (The values −0.58 / −0.28 / −0.15 of the quadratic fit are withdrawn.)
  **α changes β** (smaller damping → smaller β, at small and middle B;
  half of the α 0.01 evidence comes from the spliced runs). The earlier
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

Data note (2026-10-05, Chris's go): the 250 snapshots (.vti) of B300_s921 …
s928 and V_s901/902 were deleted on the instance after the evaluation; their
results are in runs/*/domains.json and dissmap.json (local, md5 manifest).
They were never downloaded (Chris 2026-10-04).

Data note (2026-10-06, Chris's go): the 290 snapshots of PH_s921/922 (after
phasemap, phase_extra and slices were fetched and md5-checked) and 67
checkpoints of finished runs (5.8 GB) were deleted on the instance; the 21
init.pt files stay. Disk 3.4 of 40 GB.

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

### 7.21 Jumps in the dissipation (direct test of the jump part, 2026-10-04)

**Revised 2026-10-04 after the independent review (7.24):** the statement
"90 … 99 % of the loss is smooth, i.e. a reversible response" is withdrawn.
The smooth basis contains the constant term, so a dissipation that does not
depend on the phase counts as "smooth". Fit p = c0 + c1 (d⟨M⟩/dt)² per
cycle: the phase-independent share c0/p̄ is 0.87 at 9 mT and 0.52 at
150 mT (α 0.1, n = 8; checked for s923: 0.95 → 0.57 and s928: 0.88 → 0.54),
0.8 … 1.0 at α 0.01; at the field extremes (d⟨M⟩/dt ≈ 0) p is still
0.6 … 0.9 of its mean. The loss does not follow the drive quasi-statically.
"Jumps ≤ 6 %" is a lower bound: the method counts ring-down tails and small
events as smooth. The text below is the original evaluation.

`eval_jumps.py` → results/jumps.md. Per kept cycle: p_dis(t) (256 samples)
= smooth part (Fourier series of the drive phase up to harmonic 8) +
residual; a jump is a sample with residual > 4 robust σ. Clean runs only
(B300_*, α 0.1, n = 8; A010c_*, α 0.01, n = 3); mean ± SE over the runs.

| B_peak | jumps per cycle α 0.1 / 0.01 | jump energy share α 0.1 / 0.01 | non-smooth share α 0.1 / 0.01 |
|---|---|---|---|
| 9 mT | 5.0 / 7.5 | 1.4 / 2.7 % | 3.3 / 7.2 % |
| 20 mT | 9.4 / 6.3 | 4.4 / 2.3 % | 8.2 / 7.8 % |
| 35 mT | 9.8 / – | 5.8 / – % | 11.9 / – % |
| 50 mT | 9.1 / 2.8 | 4.1 / 0.7 % | 9.7 / 5.4 % |
| 100 mT | 6.5 / 1.5 | 2.1 / 0.3 % | 7.3 / 4.4 % |
| 150 mT | 4.3 / 0.9 | 1.0 / 0.1 % | 5.6 / 3.6 % |

- Jumps carry only 0.1 … 6 % of the dissipation; all non-smooth parts
  together ≤ 12 %. Robust against the method (harmonics 4 … 12, threshold
  3 … 5 σ: the shares change by ≤ 2 points).
- Jumps shorter than the sample interval (130 ps) are seen only partly.
  Upper bound for missed energy: the energy balance w_dis/w_loop is 1.00 ±
  0.01 in the clean stages (7.17), thus missed jumps carry ≲ 1 … 2 %.
- **The jump part cannot explain β < 2** (Chris's hypothesis in the form of
  7.20): about 90 … 99 % of the loss is a smooth function of the drive
  phase, i.e. a reversible response. Its dissipation grows slower than b²:
  the loss per unit of reversible rotation falls with the amplitude. This
  is where β < 2 comes from in this model; α changes it (7.17).
- Next possible check (not done): the harmonic content of the smooth part
  vs B_peak, to see which part of the reversible response becomes less
  lossy at large amplitude.

### 7.22 Harmonics of the reversible response (2026-10-04)

**Revised 2026-10-04 (7.24):** two statements below are identities, not
evidence for a mechanism: β = 2 + d ln sin δ / d ln B (H is sinusoidal,
∮H dB contains only the fundamental of M; w_fund/w_dis = 0.96 … 1.08) and
η = μ0 sin δ/(ωχ) (a falling η is the same as a falling sin δ). The
conclusion "internal motion" does not follow from η. Valid: μ_r = 8.0 and
M3/M1 ≤ 0.3 % (linear macroscopic response).

`eval_harmonics.py` → results/harmonics.md. Clean runs (B300_*, α 0.1,
n = 8; A010c_*, α 0.01, n = 3), cycles from 2 on, mean ± SE over the runs.

| B_peak | sin δ α 0.1 / 0.01 | μ_r α 0.1 / 0.01 | M3/M1 | η/η(9 mT) α 0.1 / 0.01 |
|---|---|---|---|---|
| 9 mT | 0.013 / 0.017 | 8.02 / 7.97 | ≤ 0.003 | 1.00 / 1.00 |
| 20 mT | 0.012 / 0.010 | 8.01 / 7.98 | ≤ 0.002 | 1.05 / 0.60 |
| 50 mT | 0.010 / 0.005 | 8.00 / 8.00 | ≤ 0.001 | 0.73 / 0.30 |
| 100 mT | 0.007 / 0.003 | 8.01 / 8.02 | ≤ 0.001 | 0.58 / 0.18 |
| 150 mT | 0.007 / 0.003 | 8.02 / 8.04 | ≤ 0.001 | 0.49 / 0.18 |

(δ: loss angle of the fundamental of ⟨M⟩ against H; μ_r = 1 + |M1|/|H1|;
η = mean dissipation power / ⟨(d⟨M⟩/dt)²⟩)

- The macroscopic response is linear and constant: μ_r = 8.00 ± 0.04 at all
  amplitudes and both α, third harmonic ≤ 0.3 %.
- The loss angle falls with B_peak. With μ const, β = 2 + d ln sin δ /
  d ln B: fit 1.74 (α 0.1, the 8 baseline runs) and 1.38 (α 0.01, the 3
  A010c runs; their β_pow is 1.42 ± 0.13), in agreement with β_pow
  1.76 ± 0.05 and 1.49 ± 0.07 (7.17).
- The dissipation per unit of the macroscopic magnetization rate (η) falls
  with B_peak to 0.49 (α 0.1) and 0.18 (α 0.01) of its 9 mT value. Thus a
  large part of the dissipation does not come from the change of ⟨M⟩ but
  from internal motion of m (precession, ripple, vortex cores) that does
  not add to ⟨M⟩; this part is excited already at small amplitude and grows
  slower than the change of ⟨M⟩. The dissipation is also not in phase with
  (d⟨M⟩/dt)²: p2/p0 = 0.15 … 0.48 instead of 1.
- **Mechanism of β < 2 in this model:** a linear macroscopic response
  (μ = 8) plus an internal dissipation that grows sub-quadratically with
  B_peak; at smaller α the internal part falls faster (η 0.18) and β is
  smaller. Not shown yet: which internal motion it is (a spatial map of
  the dissipation from the snapshots would show it).

### 7.23 Where is the loss? Dissipation map from the snapshots (2026-10-04)

**Revised 2026-10-04 (7.24):** limits of this map: (i) 4 snapshots of one
cycle: non-periodic content (ringing, drift, jumps) adds positively to
|m1|², thus R is biased high; (ii) p_fund is only ≈ 50 % of p_dis at
≥ 50 mT, where the other half sits is not known; (iii) 56 % of p_fund at
9 mT is in the 10 % cells of largest texture, where dx 3 is coarsest
(mesh dependent). The location of the loss in the internal structure is a
hint, not a result.

`dissmap_batch.py` (on the instance) and `eval_dissmap.py` →
results/dissmap.md. From the 4 snapshots of the last cycle per cell the
fundamental m1 of m at the drive frequency and its Gilbert dissipation
p_fund = α μ0 Ms/γ ω² |m1|²/2. Baseline, n = 8, mean ± SE.

| B_peak | R = p_fund/p_dis | share in outer shell r/R_p > 0.8 (vol. 49 %) | share in centre r/R_p < 0.5 (vol. 13 %) | share in 10 % cells of largest texture |
|---|---|---|---|---|
| 9 mT | 0.85 ± 0.03 | 22 % | 37 % | 56 % |
| 20 mT | 0.67 ± 0.04 | 26 % | 33 % | 54 % |
| 50 mT | 0.50 ± 0.01 | 38 % | 22 % | 49 % |
| 100 mT | 0.49 ± 0.01 | 41 % | 18 % | 42 % |
| 150 mT | 0.48 ± 0.01 | 42 % | 18 % | 36 % |

- At 9 mT the local rotation at the drive frequency explains 85 % of the
  loss. The loss sits in the particle centre (37 % in 13 % of the volume)
  and in the regions of strong texture (56 % in 10 % of the cells): the
  curling/vortex-like structure inside the particles moves and carries
  most of the loss at small amplitude.
- With growing B_peak the loss moves outward (outer shell 22 → 42 %, close
  to its volume share) and out of the textured regions (56 → 36 %): the
  bulk rotation of the particle takes over.
- Above ≈ 50 mT only half of the loss is local rotation at the drive
  frequency; the other half is motion at higher frequencies (higher
  harmonics, precession, ringing) that 4 snapshots per cycle cannot
  resolve (or the 3rd harmonic that aliases).
- Interpretation (hypothesis, consistent with 7.21/7.22): the motion of
  the internal structure is lossy and dominates at small B_peak, but it
  grows slower than the bulk rotation; thus the loss per unit of ⟨M⟩
  change falls with B_peak and β < 2. Limits: 4 samples per cycle; one
  base point; the curling/vortex identification rests on the texture and
  F90 measures (7.18, 7.19), not on an image of the structure.

### 7.24 Independent review of 7.13 … 7.23 and corrections (2026-10-04)

A context-free agent (Opus) checked all results against the raw data and
the scripts. Corrections taken over (sections updated): §2, §3, §4, 7.15,
7.17, 7.21 … 7.23 (see the revision notes there). Further findings:
- Phase-independent dissipation (7.21 revision): β < 2 comes from a
  dissipation that does not follow the drive; it grows as B^1.56 ± 0.07,
  the phase-locked part as B^2.28 ± 0.03 (n = 8, α 0.1). Candidates (not
  tested): local responses with spread phases, ring-down of internal modes,
  many sub-threshold events. A broad distribution of relaxation times can
  also give a loss per cycle that does not depend on f.
- The α-independent parts (P at 10 mT ×1.06 ± 0.19 in ln; c1 of 7.20) hint
  at an irreversible part that the jump count does not see.
- **Anisotropy hint from existing data (not reported before):** F_corner_f
  (L_eff 30 and r_p 3) vs F_base_f, n = 4 + 4, 20 … 150 mT (9 mT not
  valid, unconverged relaxation, two factors changed together):
  Δβ = −0.04 ± 0.06; Δ ln P = +0.09 ± 0.09 / +0.02 ± 0.05 / +0.02 ± 0.04 at
  20 / 50 / 100 mT. Hint: no anisotropy lever above 20 mT. (Old,
  non-isotropic stress axes, 7.36: the r_p part of this hint is not valid.)
- β is higher at small B than at 50 mT: paired β(9-35) − β(35-70) =
  +0.24 ± 0.11 (n = 11, p ≈ 0.05); without s921 +0.15 ± 0.07 (p ≈ 0.05).
  A hint; 7.16 "no curvature" came from the quadratic fit.
- Splice check (7.16): n = 3 per α; the 95 % interval of the 20 mT
  agreement is ± 0.14 in ln w (α 0.1) and ± 0.19 (α 0.01), up to ≈ ± 0.1 in
  β at 9-35 mT.
- 3-point windows have 1 degree of freedom and share points (correlated);
  median χ² 3.0 in the 9-35 mT window.
- domains.py (per-particle smoothing) is correct.
- Missing checks before the Sobol design (decided by Chris 2026-10-04:
  Sobol on hold until they are done): (1) mesh300; (2) frequency test
  f/2, f, 2f at 9 … 20 mT with the converged state; (3) ≥ 16 snapshots per
  cycle over several cycles (local in-phase / quadrature parts) and a
  ring-down test (stop the drive at H max, record the decay of p_dis).
- **Linear cross-check (added 2026-10-04, Chris: a check only, do not over-
  interpret):** for a linear, viscously damped response below resonance
  tan δ ≈ α ω/ω_eff (Smit-Beljers / Kittel type small-signal response),
  thus tan δ should scale with α. Measured ratio tan δ(α 0.01)/tan δ(α 0.1)
  = 1.29 / 0.86 / 0.57 / 0.45 / 0.51 at 9 / 20 / 50 / 100 / 150 mT instead of
  0.1 (results/harmonics.md; clean runs, n = 8 vs 3). The loss is not a
  linear viscous damping, least of all at small B. f_eff = α f / tan δ =
  230 … 460 MHz (α 0.1) is a formal number of the linear picture, not a
  measured mode frequency.
- Reviewer errors (not taken over): "mesh300 did not run" (it runs; only
  the relaxation is done so far); "snapshots deleted locally" (they were
  never local; they stay on the instance by decision).

### 7.25 mesh300: dx 1.5 vs dx 3 at 1 µm with the converged virgin state (2026-10-05)

`eval_mesh300.py` → results/mesh300.md. d/l_ex 300, base point, α 0.1,
seeds 921 … 924, one job per amplitude (9 / 50 / 100 mT; 50 and 100 mT
start from the init.pt of the 9 mT job of the same mesh), 5 cycles, kept
from cycle 2. All 24 jobs done, no failure. Relaxation at dx 1.5: 20 500 …
25 500 iterations (9 … 11 h).

| | 9 mT | 50 mT | 100 mT |
|---|---|---|---|
| Δ ln w = ln w(dx 3) − ln w(dx 1.5), paired n = 4 | +0.36 ± 0.24 | +0.26 ± 0.16 | +0.19 ± 0.10 |
| same without s921 (n = 3) | +0.59 ± 0.06 | +0.41 ± 0.09 | +0.28 ± 0.04 |

| segment | β dx 1.5 | β dx 3 | Δβ paired (n = 4) |
|---|---|---|---|
| 9 → 50 mT | 1.92 ± 0.08 | 1.87 ± 0.13 | −0.05 ± 0.06 |
| 50 → 100 mT | 1.81 ± 0.06 | 1.70 ± 0.06 | −0.11 ± 0.10 |
| 9 → 100 mT | 1.89 ± 0.04 | 1.82 ± 0.10 | −0.06 ± 0.06 |

- **No mesh effect on β detected** (Δβ −0.05 … −0.11, each within 1.1 SE;
  n = 4 paired). The 95 % CI (df 3) is −0.24 … +0.14 (9→50 mT), −0.43 …
  +0.21 (50→100 mT), −0.25 … +0.13 (9→100 mT): effects up to ≈ 0.2 … 0.4 are
  not excluded (corrected 2026-10-06 after the audit, 7.32).
- **The loss level at dx 3 is higher than at dx 1.5 in 3 of 4 seeds, not
  significant** (p ≈ 0.16 … 0.23; 95 % CI at 9 mT −0.39 … +1.10 in ln w):
  +20 … +45 % with all seeds (1.5 … 1.9 SE), +32 … +80 % without
  s921 (s921 at dx 3 has events: balance 1.70 at 9 mT, 2.35 at 100 mT). The
  difference is largest at small B. In the mean over all seeds the level
  at dx 3 is +20 … +45 % above dx 1.5 (not significant; the +80 % is the
  value without s921). **This is relative to dx 1.5,
  not to the converged value:** whether dx 1.5 is converged is not known,
  and the true loss can lie below dx 1.5 or above it (corrected 2026-10-06;
  the earlier wording "too high" was not supported).
- Open question (Chris 2026-10-06): at d/l_ex 96 (old relaxation) dx 3 gave
  LESS loss than dx 1.5 (× 0.26 … 0.37, "loss missing" at the pinned vortex
  core), at 1 µm MORE. Hypothesis (not tested): an unresolved vortex core sits
  in the cell grid; if the drive cannot move it (small particle) its loss is
  missing; if it hops from cell to cell (1 µm, more structure in motion) each
  hop over the grid barrier dissipates extra energy (a Peierls-like mesh
  artifact). This fits the largest difference at small B; the event count of
  7.27 does not resolve such small hops.
- Events occur on both meshes (dx 1.5: s922 at 9 mT 1.56, s924 at 50 mT
  1.42; dx 3: s921).
- Rule of PLAN step 1b at 50 mT (same mesh error at d/l_ex 212 and 300):
  +0.26 ± 0.16 vs −0.06 ± 0.14 (7.10), |difference| + 1.645 SE = 0.68 > 0.10:
  not met. In the Sobol design an effect of d/l_ex on the loss LEVEL can be
  a mesh effect; for β the mesh error is small.
- Protocol check at dx 3: one job per amplitude from the virgin state vs
  amplitudes in sequence (baseline, same seeds): Δ ln w = −0.02 ± 0.02 /
  −0.07 ± 0.05 / −0.02 ± 0.01 at 9 / 50 / 100 mT: no history effect.
- §4 updated: the mesh decision rests now on this section.

### 7.26 Checks before the Sobol design: frequency, phase map, ring-down (2026-10-06)

`eval_checks.py` → results/checks.md; `phasemap_batch.py` on the instance.
d/l_ex 300, dx 3, base point, α 0.1, start from the converged init.pt of
B300_s921 … s923 (s921, s922 for the phase map). All 8 jobs done, no failure.

**A. Frequency test at small B** (n = 3 seeds; f = 30 MHz is B300):

| B | w(f/2) | w(f) | w(2f) | w per cycle ∝ f^n |
|---|---|---|---|---|
| 9 mT | 7.5·10⁻⁸ | 2.8·10⁻⁷ | 3.9·10⁻⁷ | n = 1.26 ± 0.19 |
| 20 mT | 4.3·10⁻⁷ | 1.2·10⁻⁶ | 1.8·10⁻⁶ | n = 1.10 ± 0.15 |

β(9 → 20 mT) = 2.26 ± 0.26 (f/2), 1.94 ± 0.15 (f), 1.98 ± 0.22 (2f).
- The loss per cycle at small B grows with f (n ≈ 1): it is dynamic, not a
  rate-independent hysteresis. The growth is stronger from f/2 to f (factor
  3.8) than from f to 2f (factor 1.4): a relaxation-like response whose loss
  peak lies near 30 … 60 MHz (estimate from 3 frequencies).
- β at 9 → 20 mT does not differ significantly between the frequencies
  (n = 3); f/2 is higher (hint).

**B. Phase-resolved local response** (16 snapshots per cycle, last 3 cycles):

| B | R1 = p_fund/p_dis (s921 / s922) | non-periodic share | local phase lag: mean / spread |
|---|---|---|---|
| 9 mT | 0.85 / 0.95 | 0.06 / 0.00 | ≈ 0° / 34 … 59° |
| 20 mT | 0.81 / 0.51 | 0.01 / 0.15 | ≈ 0° / 32 … 37° |
| 150 mT | 0.39 / 0.44 | 0.19 / 0.05 | ≈ −3° / 20 … 24° |

- At small B the local rotation at the drive frequency carries most of the
  loss (R1 0.85 … 0.95, now without the aliasing of 4 snapshots). The cells
  respond with widely spread phases (spread 34 … 59°) whose mean is ≈ 0°;
  the macroscopic lag δ is only ≈ 0.7°. Thus the local motions largely
  cancel in ⟨M⟩, while each dissipates. (Corrected after the audit, 7.32:
  the phase spread is not a loss mechanism, the LLG dissipation does not
  depend on the phase; it describes the internal redistribution of energy.
  The lag uses only the drive component and ignores transverse core motion.
  n = 2 seeds.)
- Non-periodic motion (events, drift) carries ≤ 0.19 of the loss.
- The harmonics 2 … 7 of PH_s921/922 are **not valid** (corrected after the
  third audit, 7.33): the scripts took the drive phase from the file names,
  which hold whole degrees (22, 68, … for 22.5, 67.5, …). The uneven grid
  leaks the static part of m into k = 4 (p₄/p_dis = 44 / 15 at 9 mT). R1,
  the phase lag and R_np are not affected (the leak into k = 1 is exactly 0).
  The earlier statement "at 150 mT about half of the loss is in motion faster
  than T/16" is withdrawn: this share is not determined (the snapshots are
  deleted). Fixed for the Sobol runs (angle from the snapshot index).

**C. Ring-down** (H held at H max after the last cycle, s921 / s922):
- p_dis at the start of the hold is 0.52 … 0.90 of the cycle mean and
  decays over about 0.3 … 0.5 T (≈ 10 … 15 ns): 0.11 … 0.40 at 0.5 T,
  ≤ 0.03 at 2 T. The energy dissipated after the stop is 0.20 … 0.49 of the
  loss per cycle.
- The structure keeps relaxing for a large part of a period after the field
  stops: a slow relaxation (ns scale, not the GHz precession). Together
  with A: a relaxation process with a time near 1/(2π · 30 … 60 MHz)
  (estimate) carries the loss at small B.

**Consequences (hypotheses, base point only, n = 2 … 3; see 7.32):** the
loss at 30 MHz and small B is a dynamic loss of the internal structure
around the vortex core; three hypotheses fit the data and are not separated:
(a) a near-resonant internal mode, (b) α-independent micro-avalanches whose
number grows with dH/dt, (c) relaxation. The earlier "relaxation time close
to the drive period" is not supported (7.27). β then depends on how this relaxation changes with the amplitude,
and probably on f relative to the relaxation frequency (f/2: β 2.26, hint).
For the Sobol design: β is mesh-insensitive (7.25) and the checks show no
reason against dx 3; the frequency is a hidden lever (fixed at 30 MHz).

### 7.27 Further evaluation of the checks (2026-10-06)

Scripts: `phase_extra.py` (instance) → runs/PH_*/phase_extra.json,
slices.npz; `make_phase_slices.py` → results/phase_slices_PH_s92*.png (not
in git, the .gitignore excludes png); `eval_extra.py` → results/extra.md.

**Cut images (PH_s921, PH_s922, plane through particles 0 and 3):** every
particle holds a vortex-like state: two halves with opposite m_x, a broad
curved transition zone between them and a vortex core near the particle
centre (in-plane m curls around it). At 9 mT the change of m over the cycle
is concentrated at the vortex core (a point-like spot); at 150 mT the whole
particle moves.

**Where the loss of the local fundamental sits** (16 phases, s921 / s922):

| B | share in r/R_p < 0.5 (13 % of the volume) | share in 0.5 … 0.8 (39 %) | share in r/R_p ≥ 0.8 (49 %) | share in the 10 % cells of largest texture | lag spread there / elsewhere |
|---|---|---|---|---|---|
| 9 mT | 0.33 / 0.39 | 0.43 / 0.41 | 0.24 / 0.20 | 0.51 / 0.55 | 40 … 62° / 28 … 54° |
| 20 mT | 0.31 / 0.26 | 0.44 / 0.43 | 0.26 / 0.31 | 0.51 / 0.53 | 37 … 38° / 27 … 36° |
| 150 mT | 0.18 / 0.17 | 0.40 / 0.41 | 0.42 / 0.42 | 0.45 / 0.44 | 23 … 29° / 18 … 19° |

- Confirms 7.23 without the aliasing of 4 snapshots: at small B the loss is
  in the particle centre and in the textured region around the vortex core
  (half of it in 10 % of the cells); with growing B it moves outward. The
  phase spread is largest where the loss sits.
- Texture T(2 L_eff) and F90 do not change within a cycle (≤ 1 % at 9 and
  20 mT; F90 2.1 … 3.1 % at 150 mT).

**Ring-down fit** p(t) = p_inf + A exp(−t/τ): τ = 8.5 … 22.5 ns
(0.25 … 0.68 T), p_inf ≤ 0.04 of the cycle mean; ⟨M⟩ changes during the
hold by ≤ 0.6 % of its amplitude. The dissipation goes on after the stop
without a change of ⟨M⟩: an internal relaxation of the structure (vortex
core and its surroundings) on the 10 ns scale.

**Debye fit of the frequency test (7.26 A): not determinable.** One Debye
relaxation per seed does not fit the 3 frequencies (residual 0.2 … 0.75 in
ln w; τ runs to the bound). The statement "loss peak near 30 … 60 MHz"
(7.26) is withdrawn. Per seed the loss changes by × 0.5 … × 11 between
neighbouring frequencies (7.32); only "2 f > f/2 in 3 of 3 seeds" holds.

**Events per mesh** (mesh300 and baseline): kept cycles with balance
outside 0.9 … 1.1: dx 1.5 8 of 36 (22 %), dx 3 8 of 36 (22 %), baseline
dx 3 30 of 280 (11 %); drops of the large-angle pairs ≥ 20 %: 9 / 3 / 5.
**The events are not an artifact of the coarse mesh** (they occur at least
as often at dx 1.5); the statement "mesh dependent" in 7.17 and 7.24 is
corrected: the events belong to the model (unresolved or resolved), not to
dx 3 only.

### 7.28 Mesh: which dx is necessary and sufficient? (synthesis, 2026-10-06)

All mesh evidence with a valid state (newest first):
1. d/l_ex 300, converged virgin state, dx 1.5 vs 3, n = 4 paired (7.25):
   β: Δβ = −0.05 ± 0.06 (9→50 mT), −0.11 ± 0.10 (50→100), −0.06 ± 0.06
   (9→100). Level: dx 3 higher by +20 … +45 % (all seeds), +32 … +80 %
   (without s921), largest at small B.
2. d/l_ex 212, 50 mT, second stage (valid), dx 1.5 vs 3, n = 4 (7.10):
   level −6 % ± 14 %. β not measured with a valid 9 mT value.
3. d/l_ex 96 (330 nm, old relaxation, balance ≤ 1.5; hints, 4.2 / 7.8):
   dx 3 gives 0.26 … 0.37 of the dx 1.5 loss; dx 1 gives 1.2 … 1.6× the
   dx 1.5 loss (dx 1.5 not converged there).
4. Wall measures and F90 on synthetic states: dx 3 = dx 1.5 within 0.3 %
   (7.14); events equally frequent at dx 1.5 and dx 3 (7.27).
5. Mechanism (7.26, 7.27): at small B the loss sits at the vortex core and
   in the textured region around it; the core (size ≈ l_ex) is not resolved
   at dx 3 (3 l_ex) and only marginally at dx 1.5. This explains why the
   level error is largest at small B and why small particles (vortex
   dominated) fail at dx 3.
6. Model rules: the local wall width l_ex/√(Q_eff(1 + r_p)) is 2 cells at
   dx 3 at the Sobol corner (L_eff 12, r_p 3) and the Herzer cubes
   (L_eff ≥ 12 l_ex) must be multiples of dx: dx 3 is the coarsest mesh the
   model rules allow.

**Answer:**
- **β (the Steinmetz exponent) at d/l_ex 300 (1 µm): no mesh effect
  detected between dx 3 and dx 1.5** (|Δβ| ≤ 0.11 ± 0.10, n = 4; 95 % CI up
  to ≈ ± 0.2 … 0.4) **and dx 3 is the coarsest mesh the model rules allow**
  (wall width, cube grid). "Necessary" means only this rule; "sufficient"
  is not proven as an equivalence (corrected after the audit, 7.32).
- **β at the small end of the Sobol range (d/l_ex 212): sufficiency not
  shown** (no valid β comparison; the level error differs from 1 µm, rule of
  step 1b not met). Expected by the mechanism to hold as long as walls and
  the vortex texture span many cells; a check at d/l_ex 212 (≈ 47 V100-h)
  would close it.
- **Absolute loss density: dx 3 is not sufficient** (+20 … +80 % relative to
  dx 1.5 at 1 µm); **dx ≤ 1.5 l_ex (5 nm) is necessary; whether dx 1.5 is
  sufficient is not known** (no dx 1 run at 1 µm, ≈ 300 V100-h per
  realisation; at 330 nm dx 1.5 was not converged). Absolute values from
  dx 3 carry an uncertainty of at least ≈ 40 % (dx 3 vs 1.5: +20 … +45 %
  in the mean, not significant); their
  direction to the converged value is not known;
  ratios between parameters at the same dx are usable if the mesh error is
  similar (for d/l_ex it is not, see step 1b).
- **Small particles (d/l_ex ≲ 100, ≈ 330 nm): dx 3 is not sufficient even
  for the qualitative level** (vortex core carries the process).

### 7.29 Surface roughness with the converged state, α 0.1 (2026-10-06)

`eval_rough.py` → results/rough.md. d/l_ex 300, base point, seeds 921 … 924,
one job per amplitude, 5 cycles, paired by seed. A = staircase 10 nm
(M300_dx3), B = staircase 5 nm (M300_dx15), C = smoothed surface 10 nm
(R300f, `--surface fraction`: surface cells carry Ms, A, K scaled by their
volume fraction; same volume as A, 7.6 % partial cells; the staircase has
1.50× the sphere area).

| comparison (paired, n = 4) | Δ ln w 9 mT | Δ ln w 50 mT | Δ ln w 100 mT | Δβ 9→50 | Δβ 50→100 | Δβ 9→100 |
|---|---|---|---|---|---|---|
| B − A (finer staircase) | −0.36 ± 0.24 | −0.26 ± 0.16 | −0.19 ± 0.10 | +0.05 ± 0.06 | +0.11 ± 0.10 | +0.06 ± 0.06 |
| C − A (smoothed surface) | −0.15 ± 0.25 | +0.08 ± 0.17 | **+0.59 ± 0.10** | +0.13 ± 0.16 | **+0.74 ± 0.20** | **+0.30 ± 0.07** |

β 9→100 mT: A 1.82 ± 0.10, B 1.89 ± 0.04, C 2.13 ± 0.10.

- The smoothed surface gives much more loss at 100 mT in all 4 seeds
  (+0.38 … +0.87 in ln w, balance ≈ 1.0) and a steeper β at 50 → 100 mT
  (2.44 vs 1.70). At 9 and 50 mT the difference is within the scatter.
- Refining the staircase (A → B) moves the result in the other direction
  (−0.19 at 100 mT) and leaves β almost unchanged. Thus C is not the limit
  of a finer staircase: the two surface models differ at high B, and the
  difference is not a mesh convergence effect of the staircase.
- A first hypothesis (partial cells weakly bound because A is scaled with f)
  is **refuted by the code** (magnumnp/field_terms/exchange.py): with A·f,
  Ms·f and the harmonic mean of A to a full neighbour, the exchange field of
  a partial cell is ≈ 4A/((1 + f) μ0 Ms) Δm/dx², as strong as inside; the
  anisotropy and Zeeman fields do not depend on f either. The model is
  energy-consistent. The cause of the extra loss is not known (possibly the
  smeared distribution of the surface charges).
- **Decision (Chris 2026-10-06):** both surface models converge to the exact
  sphere for dx → 0; the staircase is the controlled approximation (its
  refinement is measured: no change of β detected from 10 to 5 nm, n = 4), while the
  smoothed model at 10 nm lies far off in the opposite direction at 100 mT.
  **The fraction model is not usable at dx 3 for high fields; the staircase
  stays.** No new smoothing variant (it would be another unvalidated model).
  At α 0.01 only the staircase runs are made (R300v01_*); the started
  fraction run R300f01_s921_b9 was stopped. Their purpose now (Chris
  2026-10-06): a **paired α comparison** (same relaxed state, mesh, surface
  and protocol as M300_dx3_*, only α differs) and the event statistics at
  α 0.01; they say nothing about the surface or the mesh at α 0.01.

### 7.30 Sobol design v2 prepared (2026-10-06, not started)

PLAN 4.0. Files: `sobol_design.py` → jobs/sobol_main.txt (96 Sobol points +
8 centre replicates, α 0.1), jobs/sobol_corners.txt (4 corners × 3 seeds),
jobs/sobol_a001.txt (first 32 points + 8 centre replicates at α 0.01,
9 cycles, from the init.pt of the α 0.1 run), results/sobol_design.csv.
- (Updated 7.33: the α 0.01 subset is reduced to the 8 centre replicates,
  124 runs in all.) All runs pass the checks of the study geometry: smallest gap 1.71
  cells (contact rule), smallest local wall width 2.00 cells, smallest
  d/L_eff 7.05; 24 main points per L_eff level; realised d/l_ex 209.7 …
  303.9 (box quantization), r_p 0.02 … 3.00, φ 0.550 … 0.690.
- Template = the baseline job (same protocol, relaxation, amplitudes).
- Snapshots: 16 phases of the last cycle at 9 / 50 / 150 mT only (new option
  `--snap_mT`, tested: snapshots only in the selected stage).
- `instance_postproc.py`: per finished run phasemap + phase_extra, check,
  md5 log, then deletion of the *.vti only (7.33: and of the acc_*.npz; tested locally: checkpoint and
  results stay). `instance/chain_sobol.sh`: main + corners in random order,
  then the α 0.01 subset (7.33: the 8 centre runs); not installed on the instance yet.
- Cost: cost model 232 (main) + 31 (corners) + 283 (α 0.01, the model
  overestimates small α ≈ 2×) V100-h; realistic ≈ 400 V100-h ≈ 40 $,
  ≈ 4 days on 4 V100. (Updated 7.33: α 0.01 only 8 centre runs; ≈ 300
  V100-h, rental ≈ 32 … 35 $, ≈ 3.2 days on 4 V100, PLAN 4.0.)

### 7.31 Paired α comparison (2026-10-06)

`eval_alpha_paired.py` → results/alpha_paired.md. R300v01_* (α 0.01,
9 cycles) vs M300_dx3_* (α 0.1, 5 cycles): the same relaxed state (init.pt),
dx 3, staircase, one job per amplitude, seeds 921 … 924; only α differs.
Paired differences α 0.01 − α 0.1, mean ± SE [95 % CI, df 3], p:

| | 9 mT | 50 mT | 100 mT |
|---|---|---|---|
| Δ ln w | −0.25 ± 0.15 [−0.74, +0.23], p 0.20 | **−0.52 ± 0.05 [−0.69, −0.36], p 0.002** | −0.46 ± 0.26 [−1.30, +0.37], p 0.17 |

| segment | β α 0.1 | β α 0.01 | Δβ (paired) |
|---|---|---|---|
| 9 → 50 mT | 1.87 ± 0.13 | 1.72 ± 0.24 | −0.15 ± 0.11 [−0.51, +0.20], p 0.26 |
| 50 → 100 mT | 1.70 ± 0.06 | 1.78 ± 0.36 | +0.08 ± 0.33 [−0.99, +1.14], p 0.84 |
| 9 → 100 mT | 1.82 ± 0.10 | 1.74 ± 0.27 | −0.09 ± 0.17 [−0.64, +0.46], p 0.65 |

- The loss level at 50 mT falls robustly with α (factor 0.59, all 4 seeds
  −0.44 … −0.67). At 100 mT 3 of 4 seeds give −0.64 … −0.78; s921 gives
  +0.32 (events: balance 2.35 at α 0.1, 1.37 at α 0.01).
- **No α effect on β is detected in the paired comparison**; the CI is
  wide (the α 0.01 runs scatter more, sd of β up to 0.7). The earlier
  claim "α changes β" (7.17: −0.27 ± 0.08) rests on the spliced runs; the
  clean data do not confirm it. Status: open.
- Events at α 0.01: 5 of 12 stages have a balance outside 0.9 … 1.1
  (s922 at 9 mT: 6.1); at α 0.1 the same seeds give 2 of 12.

### 7.32 Second independent audit and corrections (2026-10-06)

A context-free agent (Opus) audited the whole study against the raw data
and the scripts; every number it recomputed reproduced. Corrections taken
over (sections updated: 0, 7.25, 7.26, 7.27, 7.28, 7.29, 10, 12; PLAN 4.0):
- Mesh: "β within ± 0.1" and "necessary and sufficient" are equivalence
  claims that n = 4 cannot carry; now "no effect detected, 95 % CI ≈ ±0.2 …
  0.4". The level difference (dx 3 vs 1.5) is not significant (p ≈ 0.16 …
  0.23). s921 is not an event-only outlier: its relaxed state differs
  (texture +5 %, 7.19) and its low 9 mT loss repeats in B300 and mesh300.
- Mechanism (phase map, ring-down, cut images) rests on n = 2 seeds: marked
  as hypothesis. "Spread local phases" is not a loss mechanism.
- Frequency: per-seed ratios at 9 mT f/2 → f: × 0.78 / 11.4 / 7.6, f → 2f:
  × 4.9 / 0.46 / 1.25 (verified); the exponent n ≈ 1.1 … 1.3 hides this.
- §0 numbers corrected (windows are n = 11; 1.61 is the 50-70-100 window).
- §12: "linear, thus β ≤ 2" was a wrong inference (see §12).
- μ depends on the discretisation (verified, 12 cycles each at 50 mT):
  8.00 (dx 3), 7.80 (dx 1.5), 8.48 (smoothed surface); 9 mT: 8.03 / 7.77.
- Eddy-current estimate added (§10).
- Sobol automation (fixed, b437226): 2 snapshot cycles (periodic vs
  non-periodic part, R1); post-processing robust, limited to the Sobol run
  names, a broken run keeps its snapshots (R2, tested); disk guard in
  run_queue (`--min_free_gb 8`, replaced in 7.33 by the guard with
  reservation, margin 3 GB); α 0.01 jobs without the 48 h wait (R3).
- Open, decision by Chris: (a) mesh check at the corners K2 (2-cell walls)
  and K3 (1.7-cell gap) after the main run (the auditor considers the
  stopping rule in PLAN 4.0 backwards: large factor effects arise where the
  resolution is worst); (b) snapshots of some α 0.01 centre runs (test of
  |m1|² ∝ 1/α); (c) model the residual variance over the factors instead of
  F-tests with n = 3; (d) the α 0.01 subset (its power for β is low after
  7.31); daily fetch of the small result files during the design.
  (Decided in 7.33 items 6, 7, 8, 9 and PLAN 4.0.)
- Referee view: no validation against measurements; systematic
  uncertainties of β (surface ≈ 0.3, mesh CI ≈ ± 0.3, f, α) are larger than
  the statistical precision (≈ 0.05 … 0.1); results are relative factor
  effects at a fixed discretisation (dx 3, staircase, α 0.1, 30 MHz).

### 7.33 Third independent audit: plan v2 (2026-10-06)

A context-free agent (Opus) audited PLAN 4.0, the code and the chain before
the start. The findings are handled one by one with Chris; each decision is
recorded here.
- 1. Phase-angle bug in phasemap_batch.py / phase_extra.py (harmonics k ≥ 2
  wrong, R1 / lag / R_np correct; 7.26 corrected). Fix: angle from the index
  θ_j = 2π(j+1)/N, the names are only checked; test with known harmonics.
  No rerun of PH_s921/922 (not needed for the design, Chris). Real-data
  check in the restart test below (d/l_ex 96, 11 / 150 mT): p_4/p_dis 0.001 /
  0.03, no leak.
- 2. Restart safety of the chain (a6581a1): code pinned to a commit
  (/root/sobol_commit), flock lock, start from /root/onstart.sh after an
  instance restart, logs appended. Tested on the instance: second start
  refused by the lock; chain killed (kill -9 of the process group) during
  stage 1, restarted: "[resume] from stage 1", run DONE, post-processing
  done, logs complete.
- 3. cleanup_instance.py: the checkpoint of a run with snapshots is
  deletable only with POSTPROC_DONE, "final" in domains.json and
  domains.json fetched (local md5); unit test with 5 run states.
- 4. Unattended operation: run_queue retries a failed job once (not exit
  code 3); failures of 2 different jobs on one GPU within 10 min pause it
  30 min, the 3rd burst disables it (ALERT); job timeout 12 h; every worker
  exception caught, reservation always released; status_sobol.json every
  5 min; run_loops stops with exit code 4 if m is not finite; the
  post-processing runs under a restart loop and rejects non-finite R1 / R_np
  / p_dis. Unit tests (retry, gate, timeout, exception, burst). Instance
  chain test with one good and one broken job: retry, ALERT, good run DONE
  and post-processed, status file written. The test showed two faults, both
  fixed: one job failing twice disabled a GPU (now: different jobs only, with
  pause), and the log showed exit code 0 instead of 1 (shell $? reset).
- 5. Gate failure (exit 3): one retry with relax_maxiter 200000 and
  --no_resume, then excluded and reported with its factors (PLAN 4.3a).
  relax_maxiter is now a run-control key: an --init_from state that
  converged is accepted with another limit (instance test: accepted; a
  state of another seed is still rejected). Corners K1_0 … K4_0 first in
  the queue.
- 6. Mesh-check rule fixed before the start (Chris: option A): trigger
  |range effect| ≥ 0.15 on β or ln P with the 95 % CI excluding 0; L_eff /
  r_p → K2, d / φ → K3; 6 paired seeds at dx 1.5 and 3, mesh300 protocol;
  cost K2 ≈ 25 $, K3 ≈ 6 $ (the plan said 6.5 $ per corner: wrong); run only
  with a separate go after the regression (PLAN 4.0).
- 7. Estimand (Chris: option C with the others as cross-checks): main fit
  = Gamma GLM with log link on w per amplitude → ln E[w] and β_core of the
  mean loss curve, SE by bootstrap over the runs; cross-checks OLS on ln w /
  per-run β (HC3) and Huber. Reason: a core averages ≈ 10⁹ particles; Huber
  down-weights real rare states (s921), and E[ln w] shifts β by ≈ −0.04
  (estimate, lognormal approximation).
- 8. Decision rules fixed before the start (Chris: as proposed): 2 primary
  outputs (β_core all amplitudes, ln E[w] 50 mT) × 4 factors, Holm 5 %;
  "not a lever" = 90 % CI inside ±0.15 (TOST); scatter-model slopes Holm
  5 %; α 0.01 predictions R ≈ 1 / 10 / 100 (PLAN 4.0). r_p 3 vs 6 pre-test
  noted as an open option.
- 9a. init.pt and checkpoint.pt of all runs (6.15 GB, computed) are fetched
  during the daily fetch before a checkpoint is deleted (Chris: yes).
- Document corrections (Chris: go): PLAN §1 outputs (β at 10 mT itself is
  not delivered, lowest window centre ≈ 18 mT), §2 mesh row, corner seeds,
  box edge 1.45 µm, 4.2 / 4.3 / 4.4 / §7 D1 / §8 updated to 4.0, φ dilution
  note, POSTPROC_FAILED snapshots manual; PROTOKOLL §3 damping row, §9,
  7.30, 7.32, §12.
- Contradiction check of PLAN / PROTOKOLL by a context-free agent (44
  findings, all numbers reproduced): PLAN 4.0 rewritten as one text with a
  terms block; definitions for the analysis fixed (Chris: go): runs in the
  fits = 104 main + 12 corners at α 0.1; w = w_loop; range effect = main
  fit at factor max minus min, other factors averaged over the design
  points; mesh-check trigger in the main fit with an uncorrected 95 % CI;
  classes only for the 2 primary outputs; scatter model for the 2 primary
  outputs (Holm over 4 slopes); |m₁|² for the α test from the accumulators.
  C1: a checkpoint is deletable only when it and init.pt are fetched. C2:
  run_queue keeps the retry state per job (QUEUE_STATE.json) across
  restarts. C3: mesh check runs with a 30 h job limit. C4: analysis code
  written and tested on synthetic data during the runs, before results are
  looked at. Cost: rental ≈ 32 … 35 $ (idle GPUs at queue ends); K3 mesh
  check ≈ 8 $. The α 0.01 test likely gives "mixed" (R ≈ 6 from 7.31,
  estimate).
- Re-check of the rewritten PLAN 4.0 (context-free agent): numbers fixed
  (cost per run 1.2 … 4.2 V100-h, K1 32 … 53 kJ/m³ at L 30, nominal box
  edge 1.47 µm, K1 corner cube count realised 195), terms added (job, Ms,
  B̂, w_dis, μ, r and R, x, u, marker files), mesh-check cost per seed in
  V100-h, 10 mT interpolation rule. run_queue: a later retry after a gate
  retry now resumes from the checkpoint (as after a restart).
- 9b. Per-cell accumulators (Chris: option B): run_loops --acc sums per
  cell over the cycles ≥ 2 of each measured stage the local dissipation, m
  and the fundamental m₁ from all 256 samples; acc_reduce.py reduces them on
  the instance (sum check, R1, shares by class, lag) and the large files are
  deleted. Instance test (d/l_ex 96, 2 amplitudes, with and without --acc):
  sum check p_box/p_dis = 1.000000; R1 0.960 / 0.206 vs 0.956 / 0.207 from
  the 16-phase snapshots; run time +4.9 / +4.1 %; w_dis identical; the two
  runs gave identical slices (md5), so the runs are deterministic. A syntax
  error (decorator) was found by this test; a test now compiles all study
  scripts.

### 7.34 Sobol design v2 started (2026-10-06 19:11 UTC)

Go by Chris after audit 3 (7.33); N = 96 kept (extensible to 128: points
96 … 127 complete the 128-point net). Chain on instance 53995473 (4×V100),
pinned commit 181b9ff, onstart entry, 124 jobs (104 main + 12 corners, then
8 α 0.01). First jobs K1_0 … K4_0 on GPU 0 … 3; post-processing running.
Credit at start 6.16 $ (≈ 14 h at ≈ 0.43 $/h): top-up by Chris needed
before ≈ 2026-10-07 09:00 UTC.

### 7.35 Analysis code of the Sobol design (2026-10-07, before any result is looked at)

`analyze_design.py` implements PLAN 4.0 (main fit Gamma GLM per nominal
amplitude with bootstrap over the runs; fits A (OLS, HC3) and B (Huber);
fit without flagged runs; scatter model; Holm on the 8 primary tests; TOST
±0.15; mesh-check trigger; D3 check; α 0.01 pair test with
R = 10 (R1'/R1)(w_dis'/w_dis)). ln w of each run is interpolated to the
nominal amplitudes (linear in ln B; b is within ≈ 4 % of the target).
Tests on synthetic data on the real design (known answer):
(1) constant scatter: all 8 true primary effects in the 95 % CI; exactly the
2 true levers found (r_p on β_all +0.30, φ on ln E[w] 50 mT +0.40);
(2) scatter growing with r_p, mean independent of r_p: the main fit finds no
r_p effect, fit A is biased (E[ln w]), the scatter model finds the r_p slope.
Loader check on the 13 fetched runs: all parsed, outputs finite (values not
looked at). Runtime ≈ 0.14 s per bootstrap resample (1000 → ≈ 2.5 min).
Environment pinned in requirements-analysis.txt (statsmodels 0.15.0: 0.14.4
does not import with scipy 1.18).

### 7.36 Stress-axis set not isotropic; Sobol design stopped (2026-10-07)

Finding (after Chris's question how the stress is built in): the 4 fixed
particle axes (`anisotropy_model.particle_axes`: cos θ = 1/8, 3/8, 5/8, 7/8
to the drive, azimuths 0 / 90 / 180 / 270°) are not an isotropic set. The
orientation tensor ⟨u uᵀ⟩ has the eigenvalues 0.253 / 0.285 / 0.461
(isotropic: 1/3 each) and the off-diagonal xy = 0.091; along the drive
⟨cos²θ⟩ = 0.328 (correct within 0.005). Cause: the azimuths. With the same
cos θ and the azimuths 137.5 / 208.9 / 66.8 / 299.0° the tensor is nearly
isotropic (eigenvalues 0.328 / 0.335 / 0.337, off-diagonal ≤ 0.002).
Consequence: the r_p effect of the design contains an unknown part of a
fixed net transverse anisotropy ∝ K_p. Not affected: all earlier studies at
the base point (r_p = 0). None of the three audits computed this invariant
(lesson: audits recompute physics invariants of the model inputs).
Decision Chris: stop the design, evaluate what exists, decide later. Chain
stopped 2026-10-07 (process group, onstart entry removed). State: 20 runs
finished and post-processed (fetched, md5-verified, 400 files in
runs/MANIFEST_sobol.md5); S035, S050, S068, S082 interrupted (checkpoints
kept on the instance); 4 α 0.01 and the rest not started.
Exploratory evaluation (`eval_sobol_partial.py` → results/sobol_partial.md,
results/sobol_partial_table.md; main effects only, n = 20, hints):
- β all 2.16 ± 0.25 (mean ± sd, n = 20), range 1.65 … 2.62; higher than the
  base point (β_pow 1.76 ± 0.05, n = 8). Method checked: β all equals the
  earlier β_pow (S005 1.956 / 1.955; B300 runs 2.11 / 1.74 / 1.74 / 1.64 with
  the same method). The difference comes from the factor values (most
  points have r_p > 0, larger L_eff, smaller d/l_ex); not separable at n = 20.
- No main effect on β all with a 95 % CI excluding 0 (largest: d/l_ex
  −0.34, CI −0.79 … +0.07; r_p +0.18, CI −0.19 … +0.53). Hints in windows:
  r_p on β 70-150 mT +0.57 (CI +0.10 … +1.13; contains the axis-set part),
  L_eff (ln Q) and φ on β 9-35 mT (+0.29 / −0.42, CIs just excluding 0),
  d/l_ex on ln w 50 mT +0.57 (CI +0.00 … +1.20).
- Corner replicates: K2 β all 2.31 / 2.47, K3 2.19 / 2.21.
- Accumulators: R1 0.5 … 0.98 at 9 mT, 0.24 … 0.45 at 150 mT; at 150 mT
  ≈ 0.4 … 0.56 of the loss in the shell r/R ≥ 0.8 (volume ≈ 0.49) and
  ≈ 0.67 … 0.79 in the 10 % cells with the largest p.
- The flag rule of PLAN 4.3a ("any flag at any amplitude") flags all 20 runs
  (mostly `*`, energy balance > 2 % at one amplitude; the base-point runs
  too): the rule is too wide for a sensitivity analysis; to be redefined.

### 7.37 Stress axes corrected: isotropic tetrahedral set (2026-10-07)

Decision Chris (option B of two): the 4 particle axes are the body diagonals
of a cube, turned so that the drive lies along [1, t, 0] of the cube frame,
t = (√5 − 1)/2. Then ⟨u uᵀ⟩ = I/3 exactly, and along the drive ⟨cos²θ⟩ = 1/3
and ⟨cos⁴θ⟩ = 1/5 exactly; cos θ = 0.795 / 0.188 / 0.188 / 0.795 (2 particles
each at θ = 37° and 79°). Option A (cos θ = 1/8 … 7/8 with optimised
azimuths) has no exact solution (3 free azimuths for 4 conditions; best
eigenvalues 0.328 / 0.335 / 0.337, ⟨cos⁴θ⟩ 0.190). The test now checks the
full tensor for 3 drive directions (the old test checked only ⟨sin²θ⟩ along
the drive, which is why the fault passed). Runs before this change keep
their axes in config.json (cos_theta_p). Nothing is started (decision on the
continuation open).

### 7.38 Cleanup for the restart (2026-10-07, Chris: go)

Instance: logs of the stopped design fetched (instance_logs/sobol_v2_oldaxes);
the 20 finished Sobol runs checked against the local manifest (all 387
files on the instance present locally with the same md5); the 24 Sobol run
folders (20 finished, 4 interrupted: S035, S050, S068, S082) deleted, also
the shuffled job file, status, logs, pin and lock files. Kept: the earlier
studies at the base point (2.9 GB, r_p = 0, valid). Disk 12 → 4.5 GB of 40.
Onstart entry removed (7.36). Laptop: the 20 runs and their manifest moved to
runs/_archive_sobol_v2_oldaxes (md5 check of all 400 files OK); the small
files of the 4 interrupted runs are there too (no .pt / .vti / .npz).
eval_sobol_partial.py reads the archive. Restart procedure: PLAN 4.0.

### 7.39 Decisions before the restart: r_p relative, flag rule (2026-10-07)

- r_p stays relative (Chris, option A): K_p = r_p K_eff, r_p 0 … 3. The
  ratio of the coherent particle anisotropy to the random Herzer anisotropy
  is the transferable quantity (our K_eff of 1 … 6 kJ/m³ is fictitious; real
  nanocrystalline K_eff is a few J/m³, Herzer; literature value), so the
  absolute K_p does not transfer anyway. Limit: r_p ≫ 3 (possible in real
  materials, rough estimate) is not covered at dx 3 (wall < 2 cells at
  L_eff 12); pre-test option r_p 3 vs 6 at larger L_eff kept.
- Flag rule (set by Claude from the earlier data, Chris: "aus der
  Vorerfahrung"): a run is flagged if |w_dis/w_loop − 1| > 0.20 at any
  amplitude. Data: the median of |w_dis/w_loop − 1| is 0.5 % (140 Sobol
  stages); > 20 % in 7 of 140 Sobol stages (7 of 20 runs; maximum 9.0 at
  K4_0 35 mT) and in 5 of 56 base-point stages (3 of 8 runs). Such a stage
  has a state change during the measured cycles (released stored energy
  appears in w_dis, not in w_loop); its w is not a stationary loss. The flags
  `o` (11 of 20 runs), `#` (9), `!` (8) describe the state and are reported
  only; `c` did not occur. Second sensitivity fit: the main fit with w_dis
  instead of w_loop. Implemented in analyze_design.py (test on the archive:
  a minority flagged).

### 7.40 Physics invariants of the model inputs checked (2026-10-07)

New rule (lesson of 7.36): every claimed property of a model input is
recomputed and has a test. Status (all tests pass, 46 in total):

| Input | Claimed property | Check | Result |
|---|---|---|---|
| stress axes | isotropic set | ⟨u uᵀ⟩, ⟨cos²⟩, ⟨cos⁴⟩ along the drive, 3 drive directions (test) | I/3, 1/3, 1/5 exactly (7.37) |
| cube axes | random, isotropic, uniform on the sphere | existing tests (orientation tensor, uniformity, seed dependence, mesh independence) | pass |
| random initial m | no preferred direction | block mean and ⟨m mᵀ⟩ against their statistical errors (test, 48³ blocks) | pass |
| relaxed virgin state | demagnetised | \|⟨M⟩\|/(φ Ms) of 31 runs (data) | median 0.002, max 0.006 |
| drive | pure sine from phase 0 | zero mean over 256 samples, H(T/4) = +H_amp (test) | pass |
| K_eff, K_p | K_eff = Q_eff K_d, K_p = r_p K_eff (J/m³) | config of 31 runs (data); units test (K_d, l_ex, f_M) | pass |
| dissipation | LLG Gilbert dissipation | run_loops formula = α μ₀Ms/γ ·\|dm/dt\|² (analysis scripts) = μ₀Ms H·dm/dt (test) | equal to 10⁻⁹; magnum.np γ = 2.21·10⁵ m/(A s) includes μ₀ |
| demag | true periodic, k = 0 removed | periodic slab stack: H = ∓Ms/2 inside / outside; uniform full box: H = 0; box mean of H = 0 for random m (test) | pass |
| geometry | φ exact, contact rule, box quantum | existing tests | pass |

Added after a review: the stress term of run_loops,
UniaxialAnisotropyField(Ku="Kp", Ku_axis="Kp_axis"), uses its own parameters
(test: with Ku = 0 the default term is silent; the Kp term gives
h = 2 K_p (m·u) u/(μ₀Ms); energy lower for m ∥ u: easy axis). The cube-axis
test checks the full tensor (diagonal and off-diagonal), not one component.
Earlier runs with r_p > 0 and the old axes (grep of all config.json, local
and instance): P2_rp3 (planning pilot 7.8, d/l_ex 96), F_corner_f / _2f
(frequency-test corner 7.9, L_eff 30, r_p 3) and through them the
anisotropy hint in 7.24; these passages are marked. All other earlier
studies have r_p = 0 and are not affected. Known limit (not fixed, no
trigger now): INIT_KEYS has no axis-set tag, so --init_from would accept
an old-axes init.pt with r_p > 0 (such files exist only in the local
archive).
The initial-state code was moved into `run_loops.virgin_m0` (same random
stream, no change of the states).
Instance test with commit 153f29d (d/l_ex 96, L_eff 12, r_p 3, 2 amplitudes,
--acc, snapshots): relaxation converged (3000 iterations); config
cos_theta_p = 0.7947 / 0.1876 / 0.1876 / 0.7947 (the new tetrahedral set),
K_eff 6217 J/m³, K_p 18651 J/m³ = 3 K_eff; accumulator sum checks 1.000000;
post-processing done; test folder deleted. No Sobol run folder on the
instance.

### 7.41 Sobol design v2 restarted with the isotropic stress axes (2026-10-07 08:53 UTC)

Go by Chris after the checks of 7.37 … 7.40 and a check that no old data
can interfere (no Sobol run folder, QUEUE_STATE, ALL10_DONE, shuffled job
file, status or pin on the instance; the 8 α 0.01 jobs start only from the
new C#/init.pt; locally no Sobol run in runs/). Pinned commit 2da32a8,
onstart entry, 124 jobs; K1_0 … K4_0 first; the configs show
cos_theta_p = 0.795 / 0.188 / 0.188 / 0.795. Same seeds as the stopped
design: the 20 archived points can be compared pairwise (same cubes and
initial state, only the stress axes differ). Credit 50.37 $.

### 7.42 Old vs new stress axes, paired (2026-10-08)

`eval_axes_paired.py` → results/axes_paired.md. The 20 points of the stopped
design recomputed with the isotropic axes (same seeds: same cubes and
initial m). Difference new − old, mean ± SE, n = 20 pairs: β all
−0.04 ± 0.06 (p 0.50; 95 % CI −0.18 … +0.09); β 9-35 +0.06 ± 0.09; β 35-70
−0.14 ± 0.15; β 70-150 −0.07 ± 0.11; ln w 10 / 50 / 100 mT −0.01 ± 0.13 /
−0.03 ± 0.10 / −0.05 ± 0.04; no slope over r_p (β all −0.04 ± 0.07 per
unit r_p). **No systematic effect of the old axis set detected; effects on
β all above ≈ 0.2 are excluded.** The high β (> 2) is not an artefact of
the old axes (mean β all of the 20 points: 2.16 old, ≈ 2.12 new).
Control: at r_p < 0.1 (S005, S039, S057; axes nearly irrelevant) the pairs
still differ by |Δβ all| 0.07 and |Δβ 9-35| 0.17 on average. A run with the
same seed is reproducible only with exactly the same inputs (7.33 9b: two
identical runs gave identical results); a small change of the anisotropy
moves the run into another trajectory. Thus the pairing does not remove the
realisation scatter: the sd of the pair differences (0.28 for β all) is of
the size of two independent realisations (σ ≈ 0.2 per run).

### 7.43 Interim look at the restarted design (2026-10-08, 32 runs, hints only)

`eval_sobol_partial.py --runs runs` → results/sobol_interim_20261008.md
(main effects only, 5 terms, bootstrap; no decision rules; the final
analysis follows PLAN 4.0 unchanged). n = 32 at α 0.1 (24 Sobol, 7 corner,
1 centre).
- β all 2.12 ± 0.26 (mean ± sd), 1.60 … 2.48; 20 of 32 above 2.
- Main effects on β all (range effect, 95 % CI): d/l_ex −0.43 (−0.73 …
  −0.19); ln Q_eff +0.20 (+0.02 … +0.40); φ +0.22 (−0.10 … +0.44); r_p +0.10
  (−0.16 … +0.33). Strongest in the window 35-70 mT: d/l_ex −0.96 (−1.33 …
  −0.56), ln Q_eff +0.51 (+0.14 … +0.91). Reading at d = 1 µm: higher A
  (d/l_ex 212, A 20 pJ/m) and stronger Herzer anisotropy (small L_eff) give
  a higher β. The d/l_ex effect is a mesh-check candidate (K3, PLAN 4.0).
- Loss level: ln w at 10 mT rises with d/l_ex (+0.68, +0.14 … +1.39); no
  clear effect at 50 / 100 mT. P(50 mT) 102 … 466 W/cm³ (median 184).
- Where the loss sits (accumulators, medians): R1 0.82 / 0.50 / 0.40 at 9 /
  50 / 150 mT; core r/R < 0.5 (12.5 % of the volume) holds 0.37 of the loss
  at 9 mT and 0.18 at 150 mT; shell r/R ≥ 0.8 (49 % of the volume) 0.23 /
  0.48; the 10 % cells with the largest p hold 0.72 at 150 mT.
- Flagged runs (|w_dis/w_loop − 1| > 0.20): 15 of 32 (more than the 1/3 of
  the earlier data; with 116 runs ≈ 55 flagged, ≥ 30 left for the sensitivity
  fit).
(results/sobol_partial.md regenerated with the flag rule of 7.39; the
numbers are unchanged, only the flag column.)

### 7.44 Pre-registered analysis of the α 0.1 runs (2026-10-10)

`analyze_design.py --boot 1000` → results/sobol_analysis.md / .json /
_runs.csv. 116 α 0.1 runs finished (0 failed, 0 retries); 115 usable.
Main fit Gamma GLM (15 terms), bootstrap 1000 (0 failed). The scatter
model is significant for ln Q_eff on ln w 50 mT (p 0.006), so the main fit
uses the weights 1/σ²(x). Flagged runs 50 of 115 (43 %; 65 left).

**Primary tests (range effect, 95 % CI, uncorrected p; Holm 5 % over 8):**

| output | factor | main fit | p | class | robust in fit A / B / no flags / w_dis |
|---|---|---|---|---|---|
| β all | d/l_ex | −0.53 (−0.66 … −0.43) | < 10⁻⁴ | lever | yes (−0.56 / −0.57 / −0.63 / −0.73) |
| β all | ln Q_eff | +0.10 (+0.01 … +0.21) | 0.04 | not determined | no (w_dis −0.17) |
| β all | r_p | +0.13 (+0.01 … +0.25) | 0.03 | not determined | no (no flags +0.06, w_dis −0.07) |
| β all | φ | −0.03 (−0.17 … +0.06) | 0.55 | not determined (90 % CI −0.153 … +0.043, misses "not a lever" by 0.003) | |
| ln w 50 mT | d/l_ex | +0.31 (+0.11 … +0.52) | 0.002 | lever | yes (+0.32 / +0.30 / +0.48 / +0.35) |
| ln w 50 mT | ln Q_eff | −0.04 (−0.24 … +0.13) | 0.65 | not determined | |
| ln w 50 mT | r_p | +0.00 (−0.21 … +0.24) | 0.98 | not determined | w_dis +0.36 |
| ln w 50 mT | φ | +0.09 [−0.14 without dilution] (−0.12 … +0.28) | 0.42 | not determined | |

Two of 8 are "lever"; none reaches "not a lever": the design excludes no
factor yet. ln Q_eff and r_p on β all have p 0.03 … 0.04 (do not pass Holm)
and change sign in the w_dis fit: hints that are not robust.

**Mesh caveat (pre-registered, decisive here).** Both levers are d/l_ex
effects. The earlier mesh test 7.10 (d/l_ex 212 vs 300, dx 3 − dx 1.5, 4
seeds, before the relaxation fix of 7.13) gave a mesh-induced difference
212 − 300 of −0.29 ± 0.15 in ln w 50 mT and +0.34 ± 0.30 in β 9 → 50 mT,
i.e. at dx 3 the mesh alone makes ln w rise and β fall with d/l_ex, with
the same sign and size as the measured effects (+0.31, −0.53). With 7.28
(+0.26 ± 0.16 at 300, −0.06 ± 0.14 at 212): +0.32 ± 0.21 in ln w 50 mT.
Thus the d/l_ex effects cannot be separated from a discretisation
artefact without the mesh check at K3.

**Mesh-check trigger (rule PLAN 4.0):** K3 (d/l_ex, 8 outputs; φ on ln w
100 mT) and K2 (ln Q_eff on β 70-150 mT +0.29, one secondary window; −0.03
without the flagged runs). K3 robust, K2 fragile. Decision by Chris.

**Other results.** Secondary (hints): d/l_ex on ln w 10 mT +0.97 (+0.76 …
+1.29), on β 35-70 mT −0.98; φ on ln w 100 mT +0.17 (contains the dilution
0.23). D3 check: no residual difference for d/L_eff < 10 (n = 35; Δ 0.003,
p 0.90). The w_dis fit differs strongly at 9 … 35 mT because w_dis contains
the energy released by state changes (S029 at 9 mT: w_dis/w_loop ≈ −490).

**Rule gap (to be defined, analysis not changed):** S029 was dropped as a
whole because its mean w_loop at 9 mT is negative (−3·10⁻⁸; cycles 2 … 6:
−4.2·10⁻⁷, then ≈ +7·10⁻⁸): a state change in the kept cycles released
stored energy. PLAN 4.3a excludes runs that crash, give NaN or miss an
amplitude; w ≤ 0 is not covered. 1 of 116 runs; it removes an event-type
run that the mean estimand is meant to include.

### 7.45 S029 sensitivity; mesh check K3 started (2026-10-10)

`eval_s029.py` → results/s029_sensitivity.md: S029 enters with the
stationary cycles after the event at 9 mT (cycles 3 … 6, w_loop 6.8·10⁻⁸
instead of a negative mean), n 115 → 116, bootstrap 1000. The two levers
are unchanged (d/l_ex on β all −0.53 → −0.56; on ln w 50 mT +0.31 →
+0.32); the hints on ln Q_eff and r_p change by ≤ 0.03; φ on β all moves
from "not determined" to "not a lever" (90 % CI −0.153 … +0.043 →
inside ±0.15): this class is borderline and is reported as such. The
mesh-check corners are the same. Decision (Chris): S029 is not repeated
with another seed (a replacement would drop event runs from the mean
estimand; PLAN 4.3a: no replacement points). The rule for w ≤ 0 stays
open; the pre-registered result (n = 115) is the main result, the S029
variant a sensitivity check.
Mesh check K3 (go Chris): jobs/meshK3.txt (36 jobs: seeds 2031 … 2036,
dx 1.5 and 3, 9 / 50 / 100 mT, 5 cycles, mesh300 protocol at the design
corner K3), chain instance/chain_k3.sh pinned to be178ae, waits for the end
of the Sobol chain, job limit 30 h; estimate ≈ 80 V100-h ≈ 8 $, ≈ 20 h.

### 7.46 Where the loss and the d/l_ex effect come from (exploratory, 2026-10-10)

`eval_mechanism.py` → results/mechanism.md (115 α 0.1 runs; domains.json,
acc.json; main effects with HC3; hints, 23 measures × 4 factors).
- **Walls exist in most of the design** (unlike the base point, 7.19):
  F90 of the virgin state median 27 % (3.9 … 79 %), above one vortex per
  particle (4.3 %) in 114 of 115 runs. F90 grows with L_eff (ln Q_eff
  −44 %-points over the range; partly by definition, the F90 window scales
  with L_eff), with smaller d/l_ex (−20) and with r_p (+11). The texture
  A_w/V (5.7 … 11.7·10⁻³/l_ex) grows with ln Q_eff and r_p, falls with d/l_ex.
- **The structure does not change with B̂:** A_w/V and F90 at 9 / 50 /
  150 mT equal the virgin state (as 7.19). The switched volume per half
  cycle grows ≈ ∝ B̂ (exponent median 0.98, 0.68 … 1.13).
- **The loss sits in the walls and cores:** the 10 % cells with the largest
  texture hold 0.61 of the local loss at 9 mT and 0.55 at 150 mT (median);
  the 10 % cells with the largest p hold 0.69 / 0.73.
- **The amount of wall does not set β:** adding A_w/V or F90 to the factor
  model leaves the d/l_ex effect on β all at −0.42 … −0.51 (from −0.57).
- **β follows how the switched volume grows with B̂:** adding the V_sw
  exponent shrinks the d/l_ex effect on β all from −0.57 ± 0.07 to
  −0.08 ± 0.04 (measure p ≈ 10⁻⁸⁹; partly by construction, the loss grows
  with the rotated volume). The exponent falls with d/l_ex (−0.20 over the
  range) and rises with ln Q_eff and r_p (+0.11 each). Mechanism of the
  d/l_ex lever: at large d/l_ex (low A at 1 µm) more volume moves already at
  9 mT (V_sw 9 mT +0.31 %-points over the range; ln w 10 mT +0.97), so the
  loss grows slower with B̂ and β is lower.
- Caveat: these are d/l_ex effects at dx 3; the mesh check K3 decides
  whether the switched volume at small B̂ in the smaller particles is a
  discretisation effect.

### 7.47 Absolute loss against β: the curves cross at ≈ 50 … 70 mT (exploratory, 2026-10-10)

`eval_crossover.py` → results/crossover.md, results/crossover.png (115 α 0.1
runs; bootstrap 1000; hints).
- Runs in thirds by β all (means 1.82 / 2.12 / 2.38, n 39 / 38 / 38). Mean
  loss of the high-β third against the low-β third: × 0.35 (0.31 … 0.40) at
  9 mT, × 0.47 at 20 mT, × 0.65 at 35 mT, × 0.90 (0.78 … 1.02) at 50 mT,
  × 1.17 (1.05 … 1.31) at 70 mT, × 1.34 at 100 mT, × 1.48 (1.39 … 1.58) at
  150 mT. **The mean loss curves cross between 50 and 70 mT**: below, the
  high-β configurations have the lower loss, above, the low-β ones. All
  three thirds have ≈ 180 … 200 W/cm³ at 50 mT: the curves turn around a
  pivot near 50 … 70 mT.
- Main fit, range effect of d/l_ex on ln E[w]: +1.01 (9 mT), +1.07 (20),
  +0.69 (35), +0.31 (50), −0.03 (70), −0.21 (100), −0.31 (150 mT); zero at
  ≈ 68 mT (these curves: GLM without the scatter weights; with the weights of
  the main fit, as in the deck, the zero is at ≈ 69 mT and the β end points
  are 1.84 / 2.37). ln Q_eff and r_p raise the loss only at 150 mT (+0.14 / +0.16, CI
  excludes 0); r_p changes sign at ≈ 39 mT, ln Q_eff at ≈ 103 mT (both
  small). φ (with the dilution 0.23) has no sign change.
- Consequence: a low β (large d/l_ex: large particles or small A) lowers the
  loss only above ≈ 70 mT; below ≈ 50 mT it raises it (up to × 2.7 at
  9 … 20 mT). The best configuration depends on the operating flux density.
  The d/l_ex part is subject to the mesh check K3.

### 7.48 α 0.01 centre test (2026-10-10, 8 pairs)

C0 … C7 at α 0.01 (from the init.pt of the α 0.1 run, 9 cycles) against the
α 0.1 runs of the same seed. Fetched and md5-checked (152 files).
- Loss ratio α 0.01 / α 0.1 (mean of ln, ± SE): × 0.34 (± 0.35 in ln) at
  9 mT, × 0.75 / 0.76 / 0.85 at 20 / 35 / 50 mT, × 0.68 / 0.78 / 0.60 at
  70 / 100 / 150 mT (150 mT: −0.52 ± 0.06). Far from × 0.1 (loss ∝ α).
  (The base point gave × 0.59 at 50 mT, 7.31; here × 0.85 ± 0.11.)
- R = |m₁|² ratio (pre-registered test): 9 mT 1.39 (95 % CI 0.93 … 2.09) →
  "viscous R ≈ 1"; 50 mT 1.67 (1.13 … 2.47) → "mixed" (excludes 1 and 10).
  No resonance (R ≈ 100 excluded at both).
- Fundamental share R1 of the loss: α 0.1 0.90 / 0.53 / 0.42, α 0.01
  0.24 / 0.11 / 0.08 at 9 / 50 / 150 mT.
- Reading: the motion at the drive frequency is almost the same at both
  dampings, and its loss scales with α (viscous). At α 0.01 most of the
  loss comes from motion that is not at the drive frequency (fast
  precession, jumps), which depends little on α. Hence the loss falls only
  to × 0.6 … 0.85 for 10 × less damping.

### 7.49 Management deck (2026-10-10)

slide_decks/beta_powder_sobol (make_deck.py, .pptx, .pdf, 14 slides): all
figures and numbers from results/deck_data.json (deck_data.py: the main fit
with the scatter weights and 1000 bootstrap resamples) and
results/sobol_analysis.json. Reviewed by a context-free agent (47 findings,
all taken over): mesh caveat on every result slide, the weighted main fit
for the A curves, R reported as a squared-amplitude ratio, pooled repeat
scatter 0.13 (20 repeats at 5 points), exploratory parts marked as hints,
cost and dates marked as estimates, one word per meaning.

### 7.50 A regime change near 30 … 50 mT (exploratory, 2026-10-10)

`eval_transition.py` → results/transition.md (115 α 0.1 runs, hints).
- The share R1 of the loss in the local fundamental falls from 0.87 (9 mT)
  over 0.72 / 0.56 to 0.49 at 50 mT and then stays at 0.47 … 0.40. The core
  share falls continuously (0.39 → 0.18); the 10 % cells with the largest p
  hold 0.69 → 0.76 (flat above 50 mT). Energy-balance events (state changes
  in the measured cycles, > 20 %) occur mainly at 9 … 35 mT (21 / 13 / 19
  runs), few above 50 mT. Local β between neighbour amplitudes: 2.13 …
  2.17 up to 50–70 mT, then 2.01 (70–100) and 1.94 (100–150).
- The amplitude where R1 falls below 0.6 does not follow the Herzer
  anisotropy field: 27 … 33 mT for H_K = 1.3 … 8.3 kA/m (L_eff levels), and
  not r_p or φ (28 … 32 mT). It follows A: 26 mT (low third, A 9.7 … 12.5
  pJ/m) against 41 mT (high third, 15.8 … 20 pJ/m). The drive field at
  30 mT is ≈ 3 kA/m (μ_r 8).
- Reading (hypothesis): near 30 … 50 mT the response changes from mostly
  linear motion at the drive frequency to a mix with about half the loss in
  other motion (jumps, faster precession); the onset is set on the
  particle / exchange scale (d/l_ex), not by the grain anisotropy. A higher A
  shifts the onset to a higher flux density: the loss stays small longer and
  then rises faster, so β is higher (consistent with 7.46). The K3 mesh
  check shows the largest mesh error of β above 50 mT (50 → 100 mT): this
  regime needs the finer mesh.

### 7.51 Jumps in the dissipation time series, stage 1 of the test of 7.50 (exploratory, 2026-10-10)

`eval_jumps_time.py` → results/jumps_time.md, jumps_time.png. No GPU. Data:
the p_dis(t) and n60(t) samples (256 per cycle) of the 115 α 0.1 runs at
dx 3, and the K3 / base mesh pairs (dx 3 vs dx 1.5, 9 / 50 / 100 mT).
Only the last 3 cycles of each amplitude are used. Reason: in a first
version with cycles 2 … 6, 23 % of the runs at 9 mT had a field-triggered
state change in cycle 2 (p up to 5 × the later cycles; not an incomplete
relaxation, see 7.53). This made false "jumps" in the steady loop at 9 mT.

Terms. P = p_dis samples, Pm = mean cycle, Ps = smooth part of Pm (Fourier
k ≤ 6). Non-smooth share RF = Σ|P − Ps| / ΣP. Non-repeating share
RN = Σ|P − Pm| / ΣP. Sharp repeating share RH = Σ|Pm − Ps| / ΣPm. Event
coupling ρ = Spearman correlation, at the same phase, between the extra
dissipation in a sample interval and the extra change of n60 in it.

A. dx 3, medians of 115 runs (95 % CI of the median by bootstrap over runs
in results/jumps_time.md):

| B (mT) | 9 | 20 | 35 | 50 | 70 | 100 | 150 |
|---|---|---|---|---|---|---|---|
| RF | 0.041 | 0.139 | 0.230 | 0.231 | 0.217 | 0.185 | 0.153 |
| RN | 0.006 | 0.022 | 0.109 | 0.122 | 0.117 | 0.119 | 0.101 |
| RH | 0.031 | 0.119 | 0.172 | 0.181 | 0.165 | 0.139 | 0.109 |
| ρ (runs with ρ > 0) | 0.01 (69/108) | 0.05 (81/106) | 0.10 (91/114) | 0.13 (100/112) | 0.13 (104/115) | 0.11 (108/115) | 0.07 (109/115) |

- CI half-widths: RF ≤ 0.03, RN ≤ 0.05 (35 mT), ρ ≤ 0.03. The snapshot
  share R_np is 0.000 / 0.029 / 0.025 at 9 / 50 / 150 mT.
- Without the run-amplitude pairs with energy balance off (> 20 %, 0 … 8
  per amplitude) the RF medians change by ≤ 0.010.
- At dx 3 the jumps form a window: they rise from 9 to 35 mT, stay at
  35 … 70 mT and fall above. Most of the jump share repeats in each cycle
  (RH > RN). Above 35 mT an extra change of n60 comes with extra loss
  (ρ > 0 in 80 … 95 % of the runs).

B. Amplitude B_pk of the largest RF per run (dx 3, n 115): median 42 mT
(36 … 51). 8 runs have B_pk at 9 or 150 mT. Linear model of ln B_pk on the 4
coded factors (OLS, HC3, R² 0.21), factor on B_pk from the low to the high
end, 95 % CI:

| factor | factor on B_pk | p |
|---|---|---|
| A 20 → 10 pJ/m (d/l_ex 212 → 300) | 0.46 (0.33 … 0.62) | < 0.001 |
| L_eff 30 → 12 l_ex (K_eff × 6.25) | 1.41 (1.07 … 1.87) | 0.015 |
| r_p 0 → 3 | 1.32 (0.91 … 1.93) | 0.15 |
| φ 0.55 → 0.69 | 0.99 (0.72 … 1.36) | 0.96 |

- A thirds: median B_pk 32 (28 … 39) / 36 (32 … 46) / 56 (51 … 65) mT for
  low / mid / high A. This agrees with the R1 onset of 7.50 (26 vs 41 mT).
  In contrast to 7.50, K_eff has a small effect (× 1.41 for × 6.25 in K_eff).

C. Mesh check (one job per amplitude, 5 cycles, last 3 used; medians over
seeds; sign test over the seeds of both points, two-sided, exact):

| B (mT) | RF dx 3 (K3 / base) | RF dx 1.5 (K3 / base) | seeds with RF(dx 3) > RF(dx 1.5) | n60 max dx 3 / dx 1.5 (medians K3 / base) |
|---|---|---|---|---|
| 9 | 0.037 / 0.071 | 0.110 / 0.114 | 1/10, p = 0.021 | 23 / 28 vs 11 / 6 |
| 50 | 0.246 / 0.230 | 0.119 / 0.114 | 7/9, p = 0.18 | 25 / 22 vs 3 / 2 |
| 100 | 0.214 / 0.166 | 0.068 / 0.082 | 9/9, p = 0.004 | 25 / 24 vs 2 / 3 |

- K3 dx 1.5 s2036 at 50 / 100 mT was still running (5 K3 and 4 base pairs
  at 50 / 100 mT). At dx 1.5, 2 of 9 seeds at 50 mT have the energy balance
  off (RF 0.88 and 0.62). This is a state change in the measured cycles or a
  large jump between the samples. The data do not separate the two.
- The dx 3 mesh jobs agree with the dx 3 Sobol runs (RF 0.04 / 0.23 / 0.19).
  Thus the different history (direct start at the amplitude) does not cause
  the difference. The mesh causes it.
- n60 must fall at the finer mesh even for the same structure (the same
  rotation spreads over twice the number of cells, half the angle per pair).
  Thus the drop of n60 is descriptive, not a separate proof. A Bloch point
  would still give pairs > 60° at dx 1.5.

Reading (three levels):
- Shown at 100 mT: the dx 3 jump share is mostly mesh-made (RF 0.17 … 0.21
  → 0.07 … 0.08, 9 of 9 seeds, p = 0.004).
- Indicated at 50 mT: the medians halve (0.23 … 0.25 → 0.11 … 0.12), but
  7 of 9 seeds, p = 0.18, and 2 seeds go strongly the other way.
- Not tested at 20, 35 and 70 mT (no dx 1.5 data). The time series thus do
  not show whether the window at 35 … 50 mT is physical.
- The A dependence of B_pk and of the R1 onset (7.50) is shown only for dx 3.
  R1 at dx 1.5 is not known (the mesh jobs ran without accumulators).
- β difference K3 − base (preliminary K3 evaluation, 34 of 36 jobs, seed
  means ± SE, n 6 / 4 at dx 3 and 5 / 4 at dx 1.5):

  | segment | dx 3 | dx 1.5 |
  |---|---|---|
  | β 9 → 50 mT | +0.46 ± 0.15 | +0.39 ± 0.13 |
  | β 50 → 100 mT | +0.79 ± 0.16 | +0.28 ± 0.21 |
  | β 9 → 100 mT (fit) | +0.54 ± 0.11 | +0.36 ± 0.07 |

  The A effect on β below 50 mT does not depend on the mesh. Above 50 mT
  it shrinks at dx 1.5 to a value that is not different from 0. This agrees
  with the mesh-made jumps at 100 mT.
- Consequence for stage 2 (frequency test f/2, 2f at dx 3): it would test
  the dx 3 jumps, which are at least partly mesh-made. Not started. Options
  go to Chris.
- Update with all 36 K3 jobs (2026-10-10 22:10 UTC, results/jumps_time.md):
  seeds with RF(dx 3) > RF(dx 1.5) 1/10 (p 0.021) at 9 mT, 8/10 (p 0.11) at
  50 mT, 10/10 (p 0.002) at 100 mT. The three levels above do not change.

### 7.52 Mesh check K3: final result, pre-registered comparison (2026-10-10)

All 36 jobs done (21:57 UTC, 0 failed). The last 2 runs fetched without
state files, md5 verified (runs/MANIFEST_k3.md5, 180 entries, all OK).
`eval_k3.py` → results/k3_mesh.md, k3_mesh.png, k3_beta.png. K3: d/l_ex 212
(A 20 pJ/m), L_eff 18, r_p 1.5, φ 0.69, 6 seeds. Base: d/l_ex 300 (A 10
pJ/m), L_eff 12, r_p 0, φ 0.65, 4 seeds (mesh300). Errors: SE over seeds,
Welch tests.

Rule A (PLAN 4.0): the mesh difference at the corner is compared with the
one at the base point. Equal → the effect is physical. Different → the
difference is the artefact share. Mesh difference = value at dx 3 − value
at dx 1.5, paired per seed.

| quantity | mesh difference K3 | mesh difference base | K3 − base (p) | verdict |
|---|---|---|---|---|
| β 9 → 50 mT | +0.05 ± 0.05 | −0.05 ± 0.06 | +0.09 ± 0.07 (0.25) | equal: physical |
| β 50 → 100 mT | +0.36 ± 0.11 | −0.11 ± 0.10 | +0.47 ± 0.15 (0.013) | different |
| β 9 → 100 mT (fit) | +0.12 ± 0.03 | −0.06 ± 0.06 | +0.18 ± 0.07 (0.046) | different |
| ln w 9 mT | −0.48 ± 0.09 | +0.36 ± 0.24 | −0.83 ± 0.25 (0.031) | different |
| ln w 50 mT | −0.38 ± 0.13 | +0.26 ± 0.16 | −0.64 ± 0.21 (0.021) | different |
| ln w 100 mT | −0.14 ± 0.09 | +0.18 ± 0.10 | −0.32 ± 0.13 (0.046) | different |

The contrast K3 − base itself, on each mesh:

| quantity | dx 3 | dx 1.5 |
|---|---|---|
| β 9 → 50 mT | +0.46 ± 0.15 (p 0.030) | +0.37 ± 0.12 (p 0.014) |
| β 50 → 100 mT | +0.79 ± 0.16 (p 0.002) | +0.32 ± 0.18 (p 0.13) |
| β 9 → 100 mT (fit) | +0.54 ± 0.11 (p 0.011) | +0.36 ± 0.06 (p 0.001) |
| ln w 9 mT | −1.03 ± 0.34 (p 0.046) | −0.20 ± 0.18 (p 0.31) |
| ln w 50 mT | −0.22 ± 0.16 (p 0.21) | +0.42 ± 0.17 (p 0.052) |
| ln w 100 mT | +0.33 ± 0.11 (p 0.027) | +0.65 ± 0.08 (p < 0.001) |

Reading:
- β: the higher β at high A (small d/l_ex) is physical. On the fine mesh
  it is +0.36 ± 0.06 (9 → 100 mT) instead of +0.54 ± 0.11: about one third
  of the dx 3 effect is a mesh artefact. Below 50 mT the effect does not
  depend on the mesh. Above 50 mT most of it is mesh-made (+0.79 → +0.32,
  not different from 0 at dx 1.5). This agrees with the mesh-made jumps at
  100 mT (7.51).
- Loss level: the d/l_ex effect on ln w at dx 3 is mostly a mesh artefact.
  On the fine mesh, the high-A point K3 has about the same loss at 9 mT
  and a higher loss at 50 mT (× 1.5, p 0.052) and 100 mT (× 1.9) than the
  base point. Thus on the fine mesh the loss curves cross below 50 mT,
  not at ≈ 57 … 69 mT as at dx 3 (7.47).
- Limit: K3 and base also differ in L_eff, r_p and φ. The Sobol fit at dx 3
  found no lever in these 3 factors (7.44), but their mesh error is not
  known. The K2 corner (L_eff, r_p) is triggered and not run.
- Consequence for the Sobol results (7.44, 7.47): the sign and the
  significance of the d/l_ex effect on β stay. Its size at dx 3 is too
  large by ≈ 1/3 (estimate from one corner). The d/l_ex effect on
  ln E[w] at 50 mT (+0.31 at dx 3) and the crossing at ≈ 57 … 69 mT are
  not supported on the fine mesh.

Fine-mesh fill-in at K3 (go Chris 2026-10-10, option 2 of 7.51): chain
`instance/chain_k3fill.sh`, commit 816b482, jobs/meshK3fill.txt (35, 20,
70 mT; dx 1.5 first, then dx 3; 6 seeds; start from the K3 9 mT init.pt;
with --acc). Local test before the start: a small job (d/l_ex 48, φ 0.55)
with --init_from and --acc, then instance_postproc → acc.json, sum check
1.000000, acc file deleted, other run not touched. Started 2026-10-10
22:03 UTC on all 4 GPUs. The init.pt of the 9 mT jobs loaded. Estimate:
≈ 13 h, ≈ 5.3 $ at 0.41 $/h. Credit 6.19 $ at 21:44 UTC.

### 7.53 Extra dissipation in the first cycles at 9 mT: a field-triggered state change, not an incomplete relaxation (2026-10-10)

Question of Chris: why is a state not relaxed? Check of the 115 α 0.1
Sobol runs, first amplitude (9 mT), from stdout.log, summary.json and the
samples. No GPU.

- The relaxation converged in all 115 runs: 19 500 iterations (median),
  47 000 (max) of the limit 80 000. Last relative energy change ≤ 9.7·10⁻⁷
  per chunk of 500 iterations, two times in a row. Last torque
  |m × H_eff| rms 0.49 A/m (median), 1.53 A/m (max).
- Upper limit of the energy that this rest torque can release at α 0.1:
  p = α γ μ0 Ms ⟨|m × H|²⟩ / (1 + α²) → 2.9·10⁻⁹ K_d per cycle (max run).
  The steady loop loss at 9 mT is ≈ 2·10⁻⁷ K_d per cycle. Thus the rest of
  the relaxation is ≈ 1 % of the loop loss.
- Events (w_dis / w_loop of one cycle > 1.2; cycle 0 > 2.2, because the
  drive ramps there): 187 events in 90 of 115 runs. Per cycle 0 / 1 / 2 / 3
  / 4 / 5: 61 / 66 / 39 / 14 / 4 / 3. The released energy (w_dis − w_loop)
  is 8.8 / 2.0 / 0.9 / 0.7 × the steady loop loss (medians). This is 10³
  × more than the rest torque can release.
- Phase: in an event cycle the largest p occurs at |H| / H_peak = 0.58
  (median). In a calm cycle (cycle 6, 40 runs) it occurs at 0.12, near the
  zero crossing of H, where dH/dt is largest. Thus the events come at a high
  field.
- Example S050: relaxation converged (rel dE 2.8·10⁻⁷, torque 1.4 A/m).
  At 9 mT, w_dis / w_loop = 3.7 / 15 / 144 in cycles 0 / 1 / 2, then
  1.0 … 1.2. Released in cycle 2: 4·10⁻⁵ K_d, 0.5 % of |E| = 7.5·10⁻³ K_d.

Reading: the relaxation at H = 0 ends in a local energy minimum, one of
many metastable states. The drive field changes the energy landscape. At a
critical field a barrier to a lower minimum vanishes, and the
magnetization jumps (vortex core, wall, Bloch point). The energy
difference is dissipated. After some cycles the state does not change any
more under ±H_peak: the steady loop. At T = 0 there is no thermal
activation, so only the field can cause this. 9 mT is the first field
after the relaxation, and its loop loss is very small: a change of
3·10⁻⁴ of |E| is already 9 × the loop loss. The hysteresis literature
calls this change of the first minor loops "accommodation" (no exact
reference checked).

Effect on the results: the analysis drops cycles 0 and 1. The flag rule
(7.39) marks a run with a balance error > 20 % in the measured cycles (50
of 115 runs flagged). 43 runs have an event in cycle ≥ 2 at 9 mT. The
d/l_ex lever without the 50 flagged runs (65 runs, sobol_analysis.json,
key noflag): β all −0.63 (−0.76 … −0.50) against −0.53 in the main fit,
ln w at 50 mT +0.48 (+0.23 … +0.72) against +0.31. Thus the events do not
make the main result. The wording "still relaxed" in the first version of
7.51 was wrong and is corrected.

### 7.54 Steady-tail rule: declared deviation from the pre-registered analysis (rule fixed 2026-10-10, before the result)

Decision Chris 2026-10-10 (option 1 of 3): no run is dropped. Only the
cycles with a state change are dropped. This rule is written down and
committed before the re-analysis runs. Until then only the numbers of
cycles were counted (7.53), no estimate was looked at.

Rule:
- Measured cycles per amplitude: as before (cycles 2 … 6 after the drive
  correction).
- Event cycle: w_loop ≤ 0 or |w_dis / w_loop − 1| > 0.20 in this cycle
  (the threshold of the flag rule 7.39, applied per cycle).
- Kept cycles ("steady tail"): the measured cycles after the last event
  cycle. Count (7.53 data): 5 / 4 / 3 / 2 / 1 / 0 kept cycles in 690 / 55 /
  31 / 16 / 9 / 4 of the 805 run-amplitudes.
- No cycle kept: the loss at this amplitude of this run is missing. The
  main fit at this amplitude uses the other runs. The per-run outputs (fit
  A, B) use the other amplitudes. The run stays in the design.
- All else as pre-registered (Gamma GLM, scatter-model weights, bootstrap
  seed and count, Holm, classes, margin 0.15).
- Report: next to the pre-registered result, as a declared deviation
  (results/sobol_analysis_steady.md / .json). The pre-registered result
  stays valid and is reported too.
- Time-series measures (7.51): the same steady tail. A run-amplitude
  counts only with ≥ 3 kept cycles (the non-repeating share RN depends on
  the number of cycles), and then with the last 3. The numbers left out
  are reported.
- Results are given as loss density P = f K_d w in W/cm³ at 30 MHz (decision
  Chris 2026-10-10); effects on ln P as a factor on P.

Code: `analyze.py` (STEADY_TAIL, steady_tail), `analyze_design.py
--steady_tail` (to_nominal, missing amplitudes in the GLM and in fit A / B),
test `test_steady_tail_rule`. Check of the code: without the option the
analysis gives the stored result again (largest difference 6·10⁻¹³). The
stored sobol_analysis.md now also holds the α section (8 pairs; it was
written before the α 0.01 runs were fetched).

Result (results/sobol_analysis_steady.md; 223 of 4060 measured cycles
dropped; 4 amplitudes missing in 4 runs). S029 now enters: its 9 mT
cycle with the negative loop area is an event cycle and is dropped, 3
cycles stay. Without S029 (115 runs) the numbers change by ≤ 0.013 and no
class changes. Primary tests, range effect (95 % CI), class:

| output × factor | pre-registered (115 runs) | steady tail (116 runs) |
|---|---|---|
| β all × L_eff 30 → 12 (lnQ) | +0.103 (+0.012 … +0.208) not determined | **+0.151 (+0.068 … +0.224) lever** |
| β all × r_p | +0.135 (+0.007 … +0.247) not determined | +0.145 (+0.016 … +0.256) not determined |
| β all × φ | −0.035 (−0.168 … +0.057) not determined | −0.065 (−0.165 … +0.023) not determined |
| β all × d/l_ex 212 → 300 | −0.531 (−0.655 … −0.434) lever | −0.572 (−0.669 … −0.476) lever |
| P 50 mT × L_eff | × 0.96 not determined | × 0.94 not determined |
| P 50 mT × r_p | × 1.00 not determined | × 1.04 not determined |
| P 50 mT × φ | × 1.09 not determined | × 1.06 not determined |
| P 50 mT × d/l_ex | × 1.36 (1.12 … 1.68) lever | × 1.36 (1.10 … 1.66) lever |

- d/l_ex stays the strongest lever; its size does not change within the
  CI.
- L_eff becomes a lever for β all: a smaller cube (higher K_eff) gives a
  higher β. The value moves by +0.05, inside its CI, but now it passes Holm.
  It is at the margin 0.15. Most of it comes from 70 … 150 mT (β@70-150 ×
  lnQ +0.31). There the dx 3 mesh is not reliable (7.51, 7.52), and the
  corner K2 (L_eff, r_p) is triggered but not run. Thus this lever is
  conditional on the mesh.
- Check of the rule: with the steady tail, the fit with w_dis agrees with
  the fit with w_loop (e.g. β@9-35 × lnQ +0.124 vs +0.125; before −0.805 vs
  0.000). In a steady cycle w_dis = w_loop must hold. Thus the rule removes
  what broke the energy balance.
- Cycle-to-cycle scatter of the loss in the steady tail: sd of ln P 2.3 …
  2.7 % (median per amplitude, 75 % quantile ≤ 5.9 %). The scatter between
  realisations is 16 … 35 % (centre, 8 seeds). One steady loop is thus
  enough for the design. Thus a run-amplitude with only 1 … 2 kept cycles is
  usable.

Time-series measures with the same rule (`eval_jumps_time.py`,
results/jumps_time.md; replaces the numbers of 7.51 where they differ):
- A (dx 3): left out 4 / 6 / 11 / 2 / 3 / 3 / 0 run-amplitudes at 9 … 150
  mT. The RF medians change by ≤ 0.012: 0.040 / 0.133 / 0.218 / 0.230 /
  0.216 / 0.184 / 0.153. The window at 35 … 70 mT stays.
- B: changes. Median B_pk 50 mT (42 … 52), before 42 mT. Factor on B_pk
  for A 20 → 10 pJ/m: × 0.67 (0.51 … 0.88), p 0.004 (before × 0.46).
  L_eff: × 1.12 (0.85 … 1.49), p 0.42 (before × 1.41, p 0.015). R² 0.07.
  A thirds 39 / 47 / 54 mT (before 32 / 36 / 56). Thus part of the earlier
  A dependence came from state changes at 20 … 35 mT, which made a false
  peak of RF. The jump peak still follows A, but weaker, and not L_eff.
- C (mesh pairs; the jobs have only 3 measured cycles, so one event cycle
  removes the pair): pairs used 7 / 5 / 9 at 9 / 50 / 100 mT. Seeds with
  RF(dx 3) > RF(dx 1.5): 0/7 (p 0.016), 5/5 (p 0.062), 9/9 (p 0.004). The
  removed cycles are real state changes (per-cycle balance 1.29 … 2.36 with
  a step in w_loop), not sampling noise. The three levels of 7.51 do not
  change: shown at 100 mT, indicated at 50 mT, not tested at 20 … 70 mT.
- Limit: R1 (7.50) comes from the accumulators, which sum all cycles ≥ 2.
  The event cycles cannot be removed there. Thus the R1 onset (26 vs 41 mT)
  can contain state changes. The K3 fill-in jobs (running) have the same
  limit.

Decision Chris 2026-10-11: the steady-tail result is the main result (deck,
reports). The pre-registered result is shown next to it as a check, and
the deviation is marked.

Kept cycles per run-amplitude (116 runs, 812 run-amplitudes): 5 / 4 / 3 / 2
/ 1 / 0 kept cycles in 696 / 55 / 32 / 16 / 9 / 4. Thus 13 run-amplitudes in
13 runs have fewer than 2 steady cycles: 9 mT S033, S035, S055; 20 mT K4_1,
S076, S091; 35 mT K4_0, S020, S026, S047, S067; 70 mT S034; 100 mT K1_0.
Some of them have a late event after clean cycles (e.g. S067 35 mT: only
cycle 5 is an event; S026 35 mT: cycles 2 and 6). Mesh and K3 jobs (3
measured cycles): 4 of 60 jobs have fewer than 2 (M300_dx15_s922_b9,
M300_dx3_s921_b100, M300_dx3_s921_b50, MK3_dx15_s2034_b50).

Check with at least 2 kept cycles (go Chris 2026-10-11, option 2;
`analyze_design.py --steady_tail --steady_min 2` → results/
sobol_analysis_steady2.md; the variant with 1 is reproduced exactly).
232 of 4060 cycles dropped, 13 amplitudes missing in 13 runs. Primary
tests on β all (on P at 50 mT nothing changes: no gap at 50 mT):

| β all × factor | pre-registered | steady, ≥ 1 cycle | steady, ≥ 2 cycles |
|---|---|---|---|
| d/l_ex | −0.531, p < 10⁻⁴, lever | −0.572, p < 10⁻⁴, lever | −0.577, p < 10⁻⁴, lever |
| L_eff (lnQ) | +0.103, p 0.039, not determined | +0.151, p 0.0002, lever | +0.152, p 0.0001, lever |
| r_p | +0.135, p 0.032, not determined | +0.145, p 0.015, not determined | +0.157, p 0.0069, **lever** |
| φ | −0.035, not determined | −0.065 (90 % CI −0.152 … +0.007), not determined | −0.058 (90 % CI −0.144 … +0.019), **not a lever** |

Holm limits at the 4th place: 0.0083 (pre-registered), 0.0100 (steady).

- Robust: d/l_ex (all variants). L_eff on β is a lever in both steady
  variants (p 10⁻⁴); the event cycles hid it in the pre-registered fit
  (p 0.039).
- At the decision limit: r_p (p 0.015 vs 0.0069 against the limit 0.010)
  and φ (90 % CI edge −0.152 vs −0.144 against the margin −0.15). Their
  class depends on the gap rule. Thus they are reported as "at the
  decision limit", not as lever or not a lever.
- The effects of L_eff and r_p on β come mostly from the windows above
  35 mT, where the dx 3 mesh has errors (7.51, 7.52). Their mesh error is
  not known (corner K2 not run).

### 7.55 Events in the measured cycles: the end of a transient, not steady physics (exploratory, 2026-10-11)

Question of Chris: are the events after apparently steady cycles part of
the physics, or states that are not yet settled? `eval_events.py` →
results/events.md. 115 α 0.1 runs × 7 amplitudes, cycles 1 … 6. The stored
energy changes per cycle by ΔE = W_loop − W_dis. In a steady state with
intermittent jumps the energy must go down and up again (ΣΔE ≈ 0, uptake as
often as release). In a transient it only goes down, and the event rate
falls with the cycle number.

| test | result |
|---|---|
| T1 event rate, cycle 1 → 6 | 9 mT: 64 → 36 → 12 → 3.5 → 2.6 → 0 %; 20 … 35 mT: 31 → 13 → 8.7 → 6.1 → 2.6 → 1.7 %; 50 … 150 mT: 7.2 → 3.5 → 2.6 → 1.5 → 0.4 → 0 % |
| T2 energy sign of the events (cycles 2 … 6) | 171 of 173 release energy, 2 take up energy. ΣΔE over cycles 2 … 6: −0.8 / −2.4 / −1.1 × loop loss with an event, −0.06 / −0.007 / −0.001 without |
| T3 after an event (one cycle before / after) | loop loss × 0.82 / 0.90 / 0.76 (median), lower in 63 … 78 %; n60 lower in 63 … 73 %, higher in 15 … 17 % (calm reference: lower 30 %, higher 28 %) |
| T4 closure of the calm cycles | median 7·10⁻⁵ (9 mT), 2·10⁻⁴ (20 … 35 mT), 6·10⁻⁴ (50 … 150 mT); numerical floor 1.5·10⁻⁷ (7.52) |
| T5 late events (cycle ≥ 4 after 2 calm cycles), n 12 | closure of the 2 calm cycles before: median 1.9·10⁻³, against 3.2·10⁻⁴ in calm run-amplitudes (Mann-Whitney p ≈ 10⁻⁵; approximate, the calm cycles are not independent) |
| T6 mesh pairs, event rate cycle 1 → 4 | dx 3: 63 → 13 → 7 → 3 %; dx 1.5: 43 → 7 → 7 → 0 %; all events at dx 1.5 release energy (17 of 17) |

Reading:
- The events are irreversible steps down in energy: the AC field anneals
  the state. Defects (vortex cores, Bloch points, counted by n60)
  annihilate, and the loop loss after the step is lower. The rate falls
  with the cycle number. This is a transient ("accommodation" after the
  start), not a steady intermittent regime (that would need energy uptake).
- A late event after apparently calm cycles has a precursor: the state
  creeps 6 × faster than in calm runs (closure) until a barrier vanishes.
  Thus these states only looked steady.
- The relaxation at H = 0 was converged (7.53). The start state is a
  minimum of the zero-field energy, but not a state that the AC drive keeps.
- The events also occur at dx 1.5, a bit less often. Thus they are not only
  a mesh effect.
- Consequence 1: dropping the event cycles (steady-tail rule, 7.54) is
  correct for the loss of a core in operation (30 MHz: 3·10⁷ cycles per
  second, a real core is conditioned).
- Consequence 2: the calm cycles still creep (closure 10⁻⁴ … 10⁻³, 10³ × the
  floor), and at 20 … 35 mT 1.7 % of the run-amplitudes still have an
  event in cycle 6. The measured state can thus be not fully conditioned.
  After an event the loss is lower. Thus the loss at 9 … 35 mT can be a bit
  too high for a long-conditioned core. The size is not known; a test with
  many more cycles at a few points can show it.

### 7.56 PLAN 4.5 started: completion to 2 calm cycles and conditioning test (go Chris 2026-10-11)

Is 2 calm cycles enough (asked by Chris)? Event in the next cycle after k
calm cycles in a row: 1.6 / 0.8 / 0.4 / 0.1 % for k = 1 … 4 (20 … 35 mT: 2.6
/ 1.7 / 0.7 / 0.3 %). A closure limit (10⁻³) would need 105 of 116 reruns.
Decision: 2 calm cycles, then the conditioning test measures the remaining
drift. The plan was checked by a context-free agent; its findings (56
stages of the 13 runs change, determinism gate, rerun mapping only in the
steady analysis, defined history and decision rule of the conditioning
test, missing init.pt and costs) are worked into PLAN 4.5.

Code (commit b10a7fa; 52 tests pass):
- run_loops.py: `--min_calm N` (stop when the calm tail ≥ N, at most
  max_cycles_per_amp; replaces the closure / dW rule when N > 0),
  `--calm_tol`, `--cycles_last_amp`. An old config without these keys counts
  as the defaults (resume).
- run_queue.py: `--stop_file`: no new job when the file exists.
- check_reruns.py: determinism gate (rerun against reference, cycle by
  cycle, w_loop / w_dis / b_peak to 10⁻⁶, up to the first stage with
  another number of cycles).
- rerun_jobs.py → jobs/rerun_mc2.txt (17 reruns: 13 Sobol, 4 mesh / K3) and
  jobs/cond.txt (18 runs cond_C0 … C5_b9 / b20 / b35). The fill-in reruns
  are made by the chain after the fill-in ends (cap 15 V100-h).
- analyze_design.py `--reruns` (only with `--steady_tail`), eval_k3.py
  `--rule steady|all` and `--no_reruns` (`--rule all --no_reruns` gives the
  tables of 7.52 again, checked), eval_jumps_time.py uses reruns.
- instance/chain_mc2.sh: waits for the fill-in chain, pins b10a7fa, one
  queue (the long mesh job first), the gate runs every 10 min, then the
  post-processing of the accumulators.

Installed on the instance 2026-10-11 ≈ 05:05 UTC: /root/mc2_commit,
onstart line, the chain waits for the lock of the fill-in chain. Fill-in at
04:48 UTC: 8 of 36 done, 0 failed. Credit 28.35 $ (top-up by Chris).

CPU test of the stop rule (small run d/l_ex 48, φ 0.55, 9 mT, 32 samples per
cycle, from the same init.pt, one after another): (a) old rule, 4 cycles;
(b) `--min_calm 2 --max_cycles_per_amp 8`: stops at 4 cycles (calm tail 2),
all numbers equal to (a); (c) as (b) with `--cycles_last_amp 6`: 6 cycles,
the first 4 equal to (a). check_reruns.py on copies (Tmc_mc2 and
cond_Tmc_b9 against Tmc): 4 cycles equal, exit code 0. A first try with the
3 runs at the same time was stopped by my own time limit (one cycle took
44 min with 3 runs on 8 cores); no code fault.

### 7.57 Fine mesh at K3 over 6 amplitudes: no jump window, R1 falls on both meshes (2026-10-11, interim)

Fill-in done 2026-10-11 10:19 UTC (36 of 36, 0 failed, post-processing
exit code 0). Fetched without state files (12 MB), 288 files md5 OK
(runs/MANIFEST_k3.md5, 468 entries). The chain of PLAN 4.5 started at 10:21
UTC (b10a7fa, 37 jobs; 2 fill-in jobs at 35 mT with a calm tail < 2 are in
it: MK3_dx15_s2034_b35, MK3_dx3_s2035_b35). Interim: these 2 jobs are
missing here until their reruns are done.

`eval_k3fill.py` → results/k3_fine.md, k3_fine.png. K3 (A 20 pJ/m), 6 seeds,
one job per amplitude from the relaxed state, steady-tail rule with ≥ 2
calm cycles. Realised B within ±1.5 % of the target in all 72 jobs.
Measurement note: `B_peak_mT` in acc.json is the largest |B| of all
samples of the stage, cycle 0 included (before the drive correction). In
these jobs it is up to 30 % above the amplitude; it is only a label. The
mapping here uses the job name. In the Sobol runs the label is within the
correct nearest amplitude (checked: S005, C0), so 7.50 is not affected.

| B (mT) | 9 | 20 | 35 | 50 | 70 | 100 |
|---|---|---|---|---|---|---|
| P dx 1.5 (W/cm³, geometric mean) | 4.3 | 22.0 | 86.6 | 197 | 503 | 948 |
| P(dx 3) / P(dx 1.5), paired | × 0.63 | × 0.77 | × 0.60 (n 4) | × 0.69 (n 5) | × 0.67 | × 0.87 |
| RF dx 3 | 0.030 | 0.129 | 0.177 | 0.261 | 0.274 | 0.221 |
| RF dx 1.5 | 0.113 | 0.134 | 0.128 | 0.118 | 0.100 | 0.074 |
| seeds RF(dx 3) > RF(dx 1.5) | 0/5 | 2/6 | 2/4 | 4/4 | 6/6 (p 0.031) | 6/6 (p 0.031) |
| R1 dx 3 | – | 0.78 ± 0.05 | 0.62 ± 0.11 | – | 0.49 ± 0.02 | – |
| R1 dx 1.5 | – | 0.71 ± 0.04 | 0.61 ± 0.04 | – | 0.57 ± 0.01 | – |

Local β (per seed, mean ± SE): dx 1.5: 2.05 ± 0.14 / 2.43 ± 0.08 / 2.32 ±
0.25 / 2.54 ± 0.22 / 1.79 ± 0.24 (9 → 20 / 20 → 35 / 35 → 50 / 50 → 70 / 70 →
100 mT); dx 3: 2.31 / 1.99 / 2.87 / 2.43 / 2.54. Paired within the seeds on
dx 1.5: β(20 → 35) − β(9 → 20) = +0.43 ± 0.17 (n 5, p 0.07); β(70 → 100) −
β(50 → 70) = −0.54 ± 0.24 (n 5, p 0.09); β(70 → 100) − β(20 → 35) = −0.43
± 0.15 (n 5, p 0.048). β fitted over 9 … 100 mT: 2.30 ± 0.04 (dx 1.5), 2.37
± 0.02 (dx 3).

Reading (hints, one corner, n 4 … 6, tests not corrected):
- On the fine mesh there is no jump window at 35 … 50 mT. The non-smooth
  share stays at 0.11 … 0.13 from 9 to 50 mT and falls to 0.07 at 100 mT.
  The rise at dx 3 above 35 mT is a mesh effect (6/6 seeds at 70 and 100
  mT). This settles the open level "not tested at 20 … 70 mT" of 7.51.
- R1 falls with B on both meshes (dx 1.5: 0.71 → 0.61 → 0.57 at 20 / 35 / 70
  mT). Thus the loss moves away from the local fundamental also on the fine
  mesh; this part of 7.50 is physical. On the fine mesh the fall is
  gradual, mostly between 20 and 35 mT, not a step at 50 mT. Limit: the
  accumulators sum all cycles from 2 on, events included; the direct start
  at each amplitude gives a lower R1 than the ascending protocol (Sobol
  corner runs K3_0 … 2 at dx 3: 0.89 / 0.83 / 0.49 at 20 / 35 / 70 mT).
- Local β on the fine mesh: about 2.0 below 20 mT, about 2.3 … 2.5 at
  20 … 70 mT, and lower above 70 mT (1.8). The dx 3 values scatter more and
  do not show this shape; above 70 mT dx 3 gives +0.75 ± 0.14 too much
  (7.52).
- At K3, dx 3 gives 13 … 40 % less loss density than dx 1.5 at every
  amplitude (P(dx 3) / P(dx 1.5) = 0.60 … 0.87).

### 7.58 β from power-law fits, both meshes (2026-10-11, interim)

Request Chris: β over the whole range and in parts from real power-law
fits, not from 2 points. Note: the design analysis (7.44, 7.54) already
uses least-squares fits (β all over 7 amplitudes, windows of 3). The
2-point values were only the local β of 7.57 and the segments 9 → 50 /
50 → 100 mT of 7.52. `eval_k3beta.py` → results/k3_beta_fits.md, .png. Fit
P = k · B^β by least squares of ln P on ln B over all amplitudes of a range
(realised B), steady-tail rule ≥ 2 calm cycles. β of the mean = fit to the
arithmetic mean of P over the seeds (the design estimand), 95 % CI by
bootstrap over the seeds (coarse with ≤ 5 seeds). Paired difference on
the seeds complete on both meshes.

| point, range (points) | seeds dx 3 / dx 1.5 / both | β of the mean dx 3 | β of the mean dx 1.5 | β per seed dx 3 | β per seed dx 1.5 | dx 3 − dx 1.5 paired (p) |
|---|---|---|---|---|---|---|
| K3, 9 … 100 mT (6) | 5 / 5 / 4 | 2.38 (2.33 … 2.43) | 2.25 (2.21 … 2.30) | 2.38 ± 0.03 | 2.27 ± 0.03 | +0.11 ± 0.02 (0.005) |
| K3, 9 … 35 mT (3) | 5 / 5 / 4 | 2.18 | 2.18 | 2.19 ± 0.08 | 2.17 ± 0.11 | +0.05 ± 0.11 (0.68) |
| K3, 20 … 50 mT (3) | 5 / 5 / 4 | 2.31 | 2.41 | 2.31 ± 0.07 | 2.39 ± 0.08 | −0.17 ± 0.06 (0.06) |
| K3, 35 … 70 mT (3) | 5 / 5 / 4 | 2.73 | 2.36 | 2.72 ± 0.09 | 2.43 ± 0.20 | +0.21 ± 0.26 (0.48) |
| K3, 50 … 100 mT (3) | 6 / 5 / 5 | 2.46 | 2.21 | 2.49 ± 0.15 | 2.27 ± 0.13 | +0.33 ± 0.13 (0.07) |
| K3, 9 / 50 / 100 mT (3) | 6 / 5 / 5 | 2.37 | 2.21 | 2.37 ± 0.03 | 2.22 ± 0.03 | +0.13 ± 0.03 (0.019) |
| base, 9 / 50 / 100 mT (3) | 3 / 3 / 2 | 1.73 | 1.87 | 1.73 ± 0.03 | 1.88 ± 0.08 | −0.10 ± 0.02 (0.14) |

(The per-mesh columns use all complete seeds of that mesh, the paired
column only the common seeds; thus the paired value is not the
difference of the two columns.)

Reading (interim; 5 jobs wait for their reruns):
- Over 9 … 100 mT the fine mesh gives a smaller β at K3: 2.25 against 2.38
  (paired −0.11 ± 0.02, n 4). The difference comes from above 35 mT (35 …
  70: +0.21, 50 … 100: +0.33 at dx 3); below 35 mT both meshes agree (2.18).
- Shape on the fine mesh: β ≈ 2.2 at 9 … 35 mT, ≈ 2.4 at 20 … 70 mT, ≈ 2.2
  at 50 … 100 mT. The 3-point windows overlap and their CIs are wide
  (± 0.1 … 0.3); the shape is a hint, not a result. dx 3 shows a peak of
  2.7 at 35 … 70 mT that the fine mesh does not show.
- A effect on β with the same amplitudes (9 / 50 / 100 mT), K3 − base: dx 3
  2.37 − 1.73 = +0.64, dx 1.5 2.21 − 1.87 = +0.34 (base n 3; interim). The
  A effect stays on the fine mesh, at about half of its dx 3 size.

### 7.59 Determinism gate stopped the queue: the start path differs (2026-10-11 11:01 UTC)

Finding: check_reruns.py found a mismatch in S026_mc2, stage 0 (9 mT),
cycle 0: w_loop 2.71872144e-07 against 2.71871676e-07 in the original
(relative 1.7·10⁻⁶), w_dis 6·10⁻⁶, b_peak 3·10⁻⁷; the corrected drive h
differs by 3·10⁻⁷. Up to cycle 6 the differences stay ≤ 10⁻⁴ (cycle scatter
2.5 %). Same GPU model and torch as the original. runs/STOP_NEW_JOBS
exists; no new job starts. Running: M300_dx15_s922_b9_mc2,
MK3_dx15_s2034_b50_mc2, S020_mc2, S026_mc2. Done: M300_dx3_s921_b50_mc2
and _b100_mc2, both equal to their originals in all 5 cycles.

Cause: the start path. The originals of S026 (all 13 Sobol runs, C0 … C5
and M300_dx15_s922_b9) relaxed in the same process and then drove the
cycles with the same LLG solver object (run_loops.py: the relaxation steps
with llg._solver at α = 1). The adaptive RKF45 solver keeps its step size
from the relaxation. A rerun with --init_from loads the same m, but its
solver starts with a fresh step size. Thus the step sequence differs and
gives differences at the level of the solver tolerance. Evidence: the
reruns whose originals also started from init.pt (M300_dx3_s921_b50 /
b100) are bit for bit equal.

Consequence: the precondition of PLAN 4.5 B (the rerun is the continuation
of the original) does not hold for the 13 Sobol reruns, for
M300_dx15_s922_b9_mc2 and for the conditioning runs (against C0 … C5). It
holds for the reruns of jobs that started from init.pt (2 M300 dx 3, 1 K3
dx 1.5 at 50 mT, 2 fill-in at 35 mT). Decision by Chris needed.

Decision Chris 2026-10-11 (option 1): every rerun takes the start path of
its original. Done:
- M300_dx15_s922_b9_mc2 (init.pt path, would run 12 … 21 h) stopped: kill of
  its own PID (110567, checked by its stdout file), 0.76 h run.
- S020_mc2 and S026_mc2 (init.pt path) run to their end as a sensitivity
  record: how much does a start difference of ~10⁻⁶ change the transient?
- rerun_jobs.py: `--init_from` only when the original had it; conditioning
  runs without it; `--retire` renames an init.pt variant X_mc2 to
  X_mc2init before a queue (not used by the analysis). Job files new; 9
  tests pass. Commit bb7bb15.
- instance/chain_mc2b.sh installed (/root/mc2b_commit = bb7bb15, onstart
  line): waits for the first mc2 chain, retires the init.pt variants,
  removes the stop file once, runs the queue with the same strict gate.
- Extra cost (estimate): relaxation in the process, ≈ 0.3 h per Sobol and
  conditioning run, ≈ 3 h for M300_dx15_s922_b9: ≈ +12 V100-h.

Sensitivity record S026_mc2 (init.pt variant, done ≈ 12:00 UTC, n 1): the
start difference of ~10⁻⁶ does not change the transient. In the stages 9 /
20 / 35 mT the per-cycle balance is the same to 2 decimals (35 mT: 1.11
1.11 8.31 1.16 1.08 1.07 1.41, the events in the same cycles); the rule
then added 2 cycles at 35 mT (calm tail 2). The later stages differ by
≤ 3 % in w_loop (other start state after the extension). At 35 mT the loss
density of the steady tail is 120 W/cm³, against 159 W/cm³ in the mean of
cycles 2 … 6 of the original (event cycle included). Thus the events are
deterministic at this level, not chaotic; one run only.

### 7.60 Barkhausen analogue from ΔM: jumps of the magnetisation (2026-10-11, exploratory)

Reason (question of Chris): RF uses p_dis as an instantaneous value at 256
samples per cycle; a jump shorter than 130 ps between two samples can be
missed. M is a state: a jump between two samples is fully in their
difference. Measure (eval_jumps_time.py, eval_k3fill.py): V = ΔM_par per
sample interval (∝ dM/dt, the induced voltage of a pick-up coil without
the smooth μ0 dH/dt); the same steps as RF: BF = Σ|V − Vs| / Σ|V|, BN =
non-repeating, BH = sharp repeating (Vs = Fourier part k ≤ 6 of the mean
cycle). A smooth loop gives BF = 0. Same runs and rule (last 3 calm
cycles). This works on all old runs; the integral of p per interval would
need new runs.

Sobol runs, dx 3, medians (95 % CI of the median), n 104 … 115:

| B (mT) | 9 | 20 | 35 | 50 | 70 | 100 | 150 |
|---|---|---|---|---|---|---|---|
| BF | 0.001 | 0.004 | 0.007 | 0.009 | 0.010 | 0.011 | 0.012 |
| BN | 0.000 | 0.001 | 0.003 | 0.004 | 0.005 | 0.007 | 0.007 |

CI half-widths ≤ 0.001. RF against BF over all run-amplitudes: Spearman
0.53 (n 776).

K3, 6 seeds (mean ± SE), BF dx 3 / dx 1.5 and seeds with BF(dx 3) >
BF(dx 1.5): 9 mT 0.001 / 0.007, 0/5; 20 mT 0.005 / 0.010, 0/6 (p 0.031);
35 mT 0.008 / 0.013, 0/4; 50 mT 0.013 / 0.016, 1/4; 70 mT 0.018 / 0.017,
3/6; 100 mT 0.019 / 0.015, 6/6 (p 0.031). Mesh pairs K3 and base: 9 mT
0/7 (p 0.016), 100 mT 9/9 (p 0.004).

Reading:
- In the magnetisation the response is 98 … 99.9 % smooth. The non-smooth
  part grows steadily with B (fine mesh at K3: 0.7 % at 9 mT, 1.6 … 1.7 % at
  50 … 70 mT, 1.5 % at 100 mT). There is no window at 35 … 50 mT, on either
  mesh. The statement of 7.57 (no jump window on the fine mesh) holds with
  this measure too.
- The RF window at dx 3 (35 … 70 mT) is not in BF at dx 3: it is a window
  of dissipation spikes without a matching jump of M (cores that hop from
  cell to cell dissipate but move little magnetisation).
- The mesh error of the jumps goes both ways: below 50 mT dx 3 has fewer
  magnetisation jumps than dx 1.5 (grid pinning: structures stick to the
  cells until a higher field), at 100 mT more. This fits the loss density
  at K3 (dx 3 × 0.63 at 9 mT, 7.57).
- The jumps are a small part of the magnetisation change (≈ 1 … 2 %) but a
  larger part of the dissipation (RF 7 … 13 % on the fine mesh). Both shares
  are L1 shares of different quantities; the ratio is a hint, not a
  measured energy share.
- Limit: a local jump without a net change of M (a turn inside a
  particle) is not seen; the time resolution stays 130 ps.

### 7.61 Loss-density slices overflowed in float16; fixed for future runs (2026-10-11)

Finding (while building the dynamics figure of the deck): in acc_slices.npz
of all runs the slice `<stage>_p` (local loss density in the plane z =
first cell layer) is inf in every particle cell. Cause: acc_reduce.py
stored p in W/m³ as float16; p is 10⁵ … 10⁹ W/m³, above the float16
maximum of 65504. Not affected: acc.json (R1, the "where" shares, the sum
check) is computed from float32 / float64 arrays before the cast; the slice
`<stage>_m1sq` (|m1|², small values) is valid. The accumulator files are
deleted after the reduction, so the old slices cannot be made again.

Fix (go Chris 2026-10-11): `<stage>_p_rel` = p / mean of the particle cells
of the plane (float16) and `<stage>_p_mean` in W/m³ (float64); helper
relative_slice, test test_acc_reduce_relative_slice_is_finite. The running
chain (pinned bb7bb15) still writes the old slice; no evaluation uses it.
The deck's dynamics figure uses the magnetisation cuts (slices.npz) instead.
End-to-end check on the CPU (small run d/l_ex 48, φ 0.55, 9 mT, --acc, then
acc_reduce.py): sum check 1.000000; p_rel finite, mean over the particle
cells 1.0, maximum 5.9; p_mean 5.2·10⁷ W/m³ (the old float16 slice would
be inf).

### 7.62 Deck rebuilt; corrections from its review (2026-10-11, interim)

Deck slide_decks/beta_powder_sobol (18 slides, levers in physical units,
main result = steady loops with ≥ 2 steady loops, dynamics cuts of the
corners K4 / K3). Reviewed by a context-free agent; all "must fix" points
are worked in. Corrections of the study record and of my earlier
statements (steady-loop rule, interim: the base point has 3 complete
repeats on the fine mesh until its rerun is done):
- The A effect on β below 50 mT does depend on the mesh. β(A = 20) −
  β(A = 10) for 9 → 50 mT: +0.59 ± 0.09 coarse, +0.33 ± 0.14 fine; 50 →
  100 mT: +0.82 ± 0.17 coarse, +0.38 ± 0.16 fine (SE, 5 … 6 and 3 … 4
  repeats). The statement "below 50 mT the effect does not depend on the
  mesh" (7.52, 7.58) came from the evaluation with all measured cycles and
  does not hold with the steady loops. On the fine mesh the A effect is
  about half in both segments.
- Loss density, fine mesh, P(A = 20) / P(A = 10): × 0.81 at 9 mT (Welch p
  0.36, 6 and 3 repeats), × 1.48 at 50 mT (p 0.09, 5 and 4), × 1.91 at 100
  mT (p < 0.001, 6 and 4). Thus only 100 mT is resolved, 50 mT is a hint.
  Coarse mesh: the curves cross at ≈ 75 mT for this pair and at ≈ 69 mT
  for all simulations (7.47).
- β of the mean (7.58) against the mean of the per-repeat fits: the deck
  uses the per-repeat means ± SE (A = 20 − A = 10 over 9 / 50 / 100 mT:
  +0.64 ± 0.05 coarse, +0.34 ± 0.08 fine). The fine-mesh A effect on β all
  in the deck, +0.31, is an estimate (coarse effect × the fine / coarse
  ratio of this pair).
- Events with S029 included (116 simulations × 7 amplitudes): 175 events in
  the measured cycles, 173 release energy, 2 take up energy (7.55 had 171
  of 173 for 115 simulations).
- Damping test, all measured cycles: β(α 0.01) − β(α 0.1) = +0.17 ± 0.11
  (SE, n 8), one pair +0.85: not resolved.
- Two different "switched" measures: the switched volume V_sw (7.14, angle
  weighted) shows the A difference at 9 mT (K4 0.67 … 0.82 %, K3 0.46 …
  0.58 % of the particle volume, 3 repeats each), the share of a cut that
  turns by more than 60° shows it at 50 mT (3.0 % against 0.5 %). The deck
  keeps the two apart and calls the mechanism a hypothesis.

### 7.63 The two meshes start from different states; same-state test started (2026-10-11)

Question of Chris: π·l_ex is normally an accepted cell size, so why does dx 3
differ so much from dx 1.5? Check of the measuring method (K3 and base, 9 mT
start jobs):
- Packing on the grid is the same (φ_vox 0.6899 / 0.6902 at K3, 0.6500 /
  0.6499 at the base point).
- The relaxed start states differ: E/K_d −0.00086 … −0.00101 on dx 3 against
  −0.00114 … −0.00118 on dx 1.5 at K3 (≈ 18 % lower on the fine mesh), base
  −0.00078 … −0.00088 against −0.00097 … −0.00098; pairs > 60° 15 … 30
  against 0 … 19.
- The pairing by seed hardly reduces the scatter: sd of the paired ln P
  difference at K3 0.21 … 0.44, against 0.25 … 0.50 expected for unpaired
  runs (base: 0.07 … 0.42 against 0.06 … 0.63).
Reading: from the same start pattern the two meshes relax into different
metastable states. The measured "mesh error" contains this difference of
the states, not only an error of the dynamics. π·l_ex is the wall width
π√(A/K) at K = K_d; the walls here are 100 … 300 nm wide, the vortex cores
≈ l_ex.

Same-state test (go Chris): coarse_init.py makes a dx 3 start state as the
normalised 2 × 2 × 2 mean of the relaxed dx 1.5 state (only fine particle
cells; the dx 3 grid has half the cells per edge, same box; config and init
record of the dx 3 template, so --init_from accepts it). Local check at the
base point: every coarse particle cell has ≥ 4 fine particle cells; mean
alignment in the blocks 0.9998. On the instance for K3, 6 seeds: alignment
0.9997, no coarse cell from the template. Jobs jobs/meshK3cg.txt (cg_jobs.py):
MK3cg_dx3_s2031 … 2036 at 9 / 50 / 100 mT, arguments of MK3_dx3 b9 with
--init_from MK3cg_init/s<seed>.pt, --min_calm 2, --max_cycles_per_amp 12.
Started 2026-10-11 ≈ 15:15 UTC as a second queue next to the PLAN 4.5 chain
(the GPUs are shared), commit d7fd4a4; the initial states loaded. Estimate
≈ 3 … 5 V100-h, ≈ 1.5 h. Decision rule (fixed now, before the result):
P(dx 3 from the fine state) against P(dx 1.5), paired per seed and
amplitude, steady loops: if the mean ln ratio is inside ±0.10 (≈ the loop
and repeat noise of the pairs) at 9, 50 and 100 mT, the dx 3 dynamics is
correct and the mesh error comes from the relaxed state; if it equals the
old dx 3 values, the dynamics itself depends on the mesh.

Result (2026-10-11, all 18 jobs;
eval_cg.py → results/k3_same_state.md). Loss density (geometric mean, steady
loops) and paired ratios per repeat (mean ln ± SE):

| B (mT) | dx 1.5 | dx 3 from the fine state | dx 3 own state | fine state / dx 1.5 | fine state / own state |
|---|---|---|---|---|---|
| 9 | 4.3 W/cm³ | 2.6 | 2.7 | × 0.61 (−0.49 ± 0.16, n 6) | × 0.98 (−0.03 ± 0.09) |
| 50 | 197 | 119 | 148 | × 0.63 (−0.47 ± 0.10, n 5) | × 0.80 (−0.22 ± 0.16) |
| 100 | 948 | 733 | 826 | × 0.77 (−0.26 ± 0.08, n 6) | × 0.89 (−0.12 ± 0.10) |

β over 9 / 50 / 100 mT: 2.22 ± 0.03 (dx 1.5), 2.32 ± 0.03 (dx 3 from the
fine state), 2.37 ± 0.03 (dx 3 own state).
- Decision rule: not met at any amplitude. dx 3 from the fine state gives
  the dx 3 values, not the dx 1.5 values.
- In the first loop the coarse-grained state releases a large energy:
  w_dis / w_loop of loop 0 = 3300 … 6800 at 9 mT, 70 … 120 at 50 mT, 18 …
  29 at 100 mT (own dx 3 state: 0.8 … 19). The fine state is not an
  equilibrium on the dx 3 grid; within one loop it falls into a dx 3 state,
  and from then on the loss is that of dx 3.
Reading: the mesh error at K3 is not a matter of the relaxation path. The
dx 3 grid cannot hold the states of the fine mesh (vortex cores ≈ l_ex
are not representable with 3 l_ex cells), and its own states have about
35 … 40 % less loss at 9 … 50 mT and a β higher by ≈ +0.1 … +0.15. One
point (K3), 5 … 6 repeats.

### 7.64 Old mesh pairs with the steady rule; decision on the size question (2026-10-11)

Question of Chris: "I want to calculate larger particles. Does the mesh
problem become smaller for larger particles?" Answer from the mesh pairs of
the base microstructure (L_eff 12 l_ex, r_p 0, φ 0.65, α 0.1) that exist,
with the steady rule (≥ 2 calm cycles). The pairs at d/l_ex 96 (P3_dx*,
R_dx15, seeds 3 … 5) and 212 (R212_*, seeds 911 … 914) have a relaxation
that did not converge (7.13) and only 9 and 50 mT. `eval_size.py` reads
them in the audit mode (hints only) → results/mesh212.md (table "old
pairs"). Errors: mean ± SE over the seeds; 90 % CI with t(n − 1).

| d/l_ex (d at A 10 pJ/m) | runs | P(dx 3)/P(dx 1.5) 9 mT | 50 mT | 100 mT | Δβ |
|---|---|---|---|---|---|
| 96 (0.32 µm) | old, n 2 … 3 | n 2: no statement | × 0.33 (ln −1.12 ± 0.67, n 3) | – | n 2: no statement |
| 212 (0.71 µm) | old, n 4 | × 0.70 (−0.36 ± 0.34) | × 0.94 (−0.06 ± 0.09) | – | 9 → 50 mT: +0.15 ± 0.23 |
| 300 (1.0 µm) | mesh300 with the reruns of mc2b (7.65), converged, n 4 | × 1.48 (+0.39 ± 0.20) | × 1.40 (+0.34 ± 0.16) | × 1.25 (+0.22 ± 0.06) | 9 / 50 / 100 mT: −0.06 ± 0.05; 9 → 50 mT: −0.02 ± 0.03 |

Also at d/l_ex 96 (old, n 3, 50 mT): dx 1 gives × 1.64 (ln +0.50 ± 0.36) of
dx 1.5, dx 2 gives × 1.30 (+0.26 ± 0.33). There is no dx 3 run at d/l_ex
150.
- Level of P: at 1 µm dx 3 gives more loss (× 1.25 at 100 mT, 90 % CI
  excludes 0; × 1.48 and × 1.40 at 9 and 50 mT, 90 % CI touches 0; values
  updated with the reruns, 7.65). At 212 and 96 the means are below 1
  (dx 3 lower), but all 90 % CIs include 0 (212: −1.15 … +0.43 at 9 mT,
  −0.28 … +0.16 at 50 mT; 96: −3.06 … +0.82). A sign change is thus not
  shown for the base microstructure; dx 3 lower at d/l_ex 212 is shown
  only for K3 (7.63). The level error at 1 µm (+30 … 60 %) does not go to
  zero. No run shows the behaviour above 1 µm.
- β: the mesh error at 1 µm is small: −0.06 ± 0.05 (9 / 50 / 100 mT,
  90 % CI −0.19 … +0.07, n 4, steady rule with the reruns, 7.65). At 212 it
  is +0.15 ± 0.23 (n 4, 9 → 50 mT, old pairs):
  not different from 0 and not different from 0.3. A trend of the β error
  with the size is not shown.
- Corrections of the chat answers of 2026-10-11: "the β error falls with
  the size" used n 2 at 96 and 300; "the level error changes sign" does not
  pass the sign rule of PLAN 4.6 (both 90 % CIs must exclude 0).
- dx 1.5 can have its own error at small particles (hint at 96).
Decision of Chris (2026-10-11): option 1 now (no new size trend; Chris
drops the prepared 150 + 212 design); option 2 into the queue after the
chain mc2b: the mesh check at d/l_ex 212 with a converged relaxation,
9 / 50 / 100 mT, 4 seeds (PLAN 4.6). A test above d/l_ex 300 (larger
than 1 µm at A 10 pJ/m) needs more GPU memory than 16 GB for dx 1.5
(estimate) and stays open.
Installed 2026-10-11 ≈ 17:50 UTC: /root/chain_size.sh with the pinned
commit 96faa32 and an onstart line. It waits for runs/MC2B_DONE. The mc2b
chain (commit bb7bb15, gate watcher, 4 jobs) runs on without change.

### 7.65 PLAN 4.5 done: acceptance, main result with the reruns, conditioning test (2026-10-12)

Chain mc2b: 37 jobs (19 reruns X_mc2, 18 conditioning runs), 2026-10-11
13:20 … 2026-10-12 12:53 UTC, commit bb7bb15. Exit codes: queue 0,
determinism gate 0, post-processing 0. The 3 jobs M300_dx3_s921_b50_mc2,
M300_dx3_s921_b100_mc2 and MK3_dx15_s2034_b50_mc2 were already done in the
first chain mc2. The long rerun M300_dx15_s922_b9_mc2 needed 19.97 h
(20.5 k relaxation iterations, 6 cycles); its cycles 0 … 4 are identical to
the original. Fetched without .pt (75 MB); logs in runs/mc2b_logs/.

**C. Acceptance: met.** Every stage of the steady analysis has a calm tail
≥ 2: 116 design runs (13 reruns used), mesh300 24 (3), K3 36 (1), K3
fill-in 36 (2), same-state test 18 (0). Exception (outside PLAN 4.5): the 8
runs at α 0.01 (C0_a001 … C7_a001) have 23 stages with a calm tail < 2;
the damping test keeps all measured cycles for them (7.62).

**Main result (steady rule, ≥ 2 calm cycles, reruns;
results/sobol_analysis_steady2_mc2.md).** Missing stages 13 → 0.

| output × factor | without reruns (7.54 check) | with reruns (main result) |
|---|---|---|
| β_all × ln Q_eff | +0.152, p 0.0001, lever | +0.150, p 0.0002, lever |
| β_all × r_p | +0.157, p 0.0069, lever | +0.145, p 0.0145, not determined |
| β_all × φ | −0.058, 90 % CI −0.144 … +0.019, not a lever | −0.064, 90 % CI −0.150 … +0.008, not determined |
| β_all × ln d/l_ex | −0.577, lever | −0.572, lever |
| ln w 50 mT × ln d/l_ex | +0.308, p 0.0022, lever | +0.308, p 0.0022, lever |
| ln w 50 mT × ln Q_eff, r_p, φ | not determined | not determined |

- The reruns move each effect by ≤ 0.012. r_p and φ on β lie on the
  decision limits (Holm 5 % over 8 tests: r_p needs p < 0.010 at its rank;
  the margin ±0.15 for φ). Their class changes with this small shift. Thus
  the design cannot decide r_p and φ on β. K_eff and d/l_ex on β and
  d/l_ex on the loss level are robust levers.
- The α 0.01 block of this analysis gives nan (no steady stages); the
  damping test stays with all cycles (7.62).

**K3 mesh check with the reruns (results/k3_mesh.md, k3_beta_fits.md,
k3_fine.md).** n now 6 / 6 / 6 at K3 and 4 / 4 / 4 at the base point.
- β over 9 … 100 mT (6 amplitudes) at K3: dx 3 2.37 ± 0.02, dx 1.5
  2.30 ± 0.04, Δβ +0.07 ± 0.03 (p 0.057; before +0.11 ± 0.02).
- β 9 / 50 / 100 mT at K3: Δβ +0.12 ± 0.03 (p 0.010).
- Base point d/l_ex 300 (n 4): Δβ −0.06 ± 0.05 (p 0.38); ln P ratio
  dx 3 / dx 1.5 +0.39 ± 0.20 / +0.34 ± 0.16 / +0.22 ± 0.06 at 9 / 50 /
  100 mT.
- A effect on β (K3 − base, 9 → 100 mT fit): +0.54 ± 0.11 at dx 3,
  +0.36 ± 0.07 at dx 1.5 (unchanged).
- Jump measures (results/jumps_time.md): changes ≤ 0.01 in the shares.

**D. Conditioning test: the design values stand at all 3 amplitudes**
(eval_cond.py → results/conditioning.md; rule of PLAN 4.5 D; n 6).

| B | mean drift d ± SE | 90 % CI | P(cycles 28 … 29) / P(design) | decision |
|---|---|---|---|---|
| 9 mT | +0.001 ± 0.013 | −0.024 … +0.027 | × 1.001 | design value stands |
| 20 mT | −0.001 ± 0.005 | −0.012 … +0.009 | × 0.999 | design value stands |
| 35 mT | −0.004 ± 0.009 | −0.022 … +0.014 | × 0.996 | design value stands |

- β 9 → 20 mT: design 2.127 ± 0.159, late 2.124 ± 0.156, difference
  −0.003 ± 0.018; β 20 → 35 mT: −0.004 ± 0.021 (n 6).
- Events: 6 in 18 runs, all in the cycles 2 … 4 (C1 9 mT: 2, 3; C5 9 mT:
  2; C5 20 mT: 3, 4; C2 35 mT: 4). No event in the cycles 5 … 29 (450
  cycles). The single-seed drifts are ≤ 0.051 in ln P.
- Reading: after the calm tail, 23 more cycles do not change the loss
  density by more than ≈ 3 % (90 % CI) and β by more than ≈ 0.04. The
  events are the end of a transient (7.55), not a slow drift. 2 calm
  cycles are enough for the design values.

### 7.66 Independent audit of all local data (2026-10-12)

A context-free agent checked all local runs, the analysis code, the result
files and the deck (request of Chris). Report: audit_20261012/REPORT.md,
scripts and tables in the same folder (commit 8e4025c). The main session
re-checked the A × K_eff interaction (own OLS: A lever +0.96 / +0.51 /
+0.15 at L_eff 30 / 18 / 12; K_eff lever +0.54 at A 10, −0.27 at A 20
pJ/m; the same as the OLS table in audit_20261012/beta_fits.md: lnd
−0.555, lnQ +0.134, lnd × lnQ +0.409) and the α 0.01 loops at 9 mT: both
hold.
- Data: 275 run folders (256 current, 19 replaced) complete and consistent (DONE, converged
  relaxation, job line = config, stdout = summary, loop integrals from the
  time curves = summary to 10⁻⁷); every stage has ≥ 2 calm cycles except
  the 8 runs at α 0.01 (exempt); 19 reruns bit-identical up to the first
  extended stage. 112 folders withdrawn (7.13). 63 folders had no md5
  manifest: runs/MANIFEST_mc2b.md5 (506 files, verified).
- Main analysis reproduced to 4 digits with an own Gamma GLM. r_p and φ on
  β lie inside the bootstrap Monte-Carlo noise of their decision limits
  (p(r_p) 0.011 … 0.019 against 0.010; φ 90 % edge −0.150 … −0.161).
- New: A × K_eff interaction on β (§0 item 3); stress-angle hint (REPORT
  D.4.6). Confirmed: the loss-angle identity β − 2 = d ln sin δ / d ln B
  (7.22, 7.24); the loss moves from the particle centre to the rim with B
  (7.23, 7.50).
- Corrections:
  - 7.48: at α 0.01 and 9 mT, 6 of 8 seeds have no calm loop; loop areas
    are negative in 5 of 8 (w_dis / w_loop −180 … +194). Decision Chris
    (2026-10-12): the damping test uses the steady rule (≥ 2 calm cycles)
    like all other results; the ratio × 0.34 at 9 mT and R(9 mT) are
    withdrawn; 30-cycle runs at α 0.01 for 9 and 20 mT (PLAN 4.7, ≈ 23
    V100-h ≈ 2.5 $, estimate; first offered as ≈ 0.5 $, corrected with the
    measured 147 s per cycle). With the new rule: × 0.79 / 0.73 / 0.75 /
    0.81 / 0.62 at 35 / 50 / 70 / 100 / 150 mT (n 5 / 5 / 6 / 6 / 4); with
    all cycles (old rule) × 0.60 … 0.84 at 20 … 150 mT (n 8).
  - 7.63: numbers before the rerun MK3_dx15_s2034_b50_mc2. Current: ln P
    (dx 3 from the fine state / dx 1.5) −0.49 ± 0.16 / −0.58 ± 0.14 /
    −0.26 ± 0.08; (fine state / own dx 3 state) −0.03 ± 0.09 / −0.22 ±
    0.16 / −0.12 ± 0.10 at 9 / 50 / 100 mT; P(dx 1.5, 50 mT) 213 W/cm³
    (n 6); β(dx 1.5, 9 / 50 / 100 mT) 2.25 ± 0.04; the own dx 3 states have
    ≈ 30 … 37 % less loss at 9 … 50 mT (results/k3_same_state.md, updated).
    The conclusion holds.
  - Deck: "higher K_eff raises β" holds only at low A; "61 % of P from 10 %
    of the volume" is the most textured 10 % (the most lossy 10 % carry
    69 %); 36 of 41 items are stale (list in REPORT.md C).
  - Code labels: in the steady mode the column "without flagged runs" of
    analyze_design.py is an unweighted refit (no run is flagged); the α
    test prints "mixed" for nan; make_deck.disp_class cannot show "at the
    decision limit" with the _mc2 analysis. No numbers change.
- No error found in: the steady rule, the rerun mapping, GLM / bootstrap /
  Holm / TOST, eval_k3*, eval_cond, the unit chain, the loop integrals.

### 7.67 Loss decomposition B¹ + B² + B³: no mechanism (2026-10-12, option A of Chris)

Question of Chris: why does β change, what is the physical background? Test
of a three-term model per design run, method fixed before the result
(eval_decomp.py → results/loss_decomposition.md): P = a1 x + a2 x² + a3 x³,
x = B / 50 mT, a_j ≥ 0, fit of the relative residual; hypotheses: x¹
switching events, x² linear viscous loss (∝ α), x³ Rayleigh (progressive
depinning). 116 runs, steady loops, coarse mesh.

- Fit quality: rms of the ln P residual 0.137 (median) with 3 terms against
  0.134 with one power law (2 parameters). The decomposition does not
  describe the curves better. Determined (a_j > 2 SE): x¹ in 17, x² in 91,
  x³ in 38 of 116 runs.
- Structural reason: with a_j ≥ 0 the local β of the model can only rise
  with B. The measured loss curves bend down (mean curvature −0.08 ± 0.02,
  β 2.15 → 1.98 above 70 mT, 7.66). Thus fixed-exponent terms cannot
  describe them. A saturating mechanism is needed.
- Damping check (centre, 4 seeds with ≥ 4 amplitudes at both α): a2 × 0.76
  for α × 0.1, not × 0.1. The x² term is not a viscous Gilbert loss. The
  fits are unstable between seeds (n 4): no further reading.
- The only robust pattern restates β: the share of the steep x³ term grows
  with A (high-A third 0.11 / 0.43 / 0.70 at 9 / 50 / 150 mT, low-A third ≈ 0;
  range effect of A on the x³ share at 150 mT +0.66 ± 0.07, with a strong
  A × K_eff term). This describes the higher β of high A; it does not
  explain it.
- Conclusion: the 30 MHz loss curves alone cannot separate the mechanisms. A
  second axis is needed: the rate dependence (frequency) or the
  rate-independent part (quasi-static loops). Decision of Chris
  (2026-10-12): frequency test C′ on the K1 pair (PLAN 4.8). The old test 7.9
  stays valid only at 20 … 150 mT at the base point (7.13): w(2 f) / w(f) =
  1.02 / 1.44 / 1.83 / 2.06 / 2.09 at 20 / 35 / 50 / 100 / 150 mT; its corner
  used the old stress axes (7.36). Rate-dependent (∝ f) but not ∝ α (7.48,
  7.66): the loss above ≈ 50 mT is dynamic, not a plain Gilbert loss of a
  linear response.
Queue on the instance (2026-10-12 ≈ 19:40 UTC): mesh check d/l_ex 212 →
damping runs at 9 / 20 mT (PLAN 4.7) → pair K1A10 (PLAN 4.8) → frequency
test C′ (PLAN 4.8); every chain waits for the DONE file of the previous one.

### 7.68 Results computed for the deck: loss densities, best material, field cuts (2026-10-12)

Source: deck_data.py → results/deck_data.json (commits c4f2ec8 … d1b7c8f);
main fit (Gamma GLM, 15 terms, steady loops with reruns), coarse mesh, P per
core volume in W/cm³ at 30 MHz. "Small particles" = A 20 pJ/m at d = 1 µm
(= 0.71 µm at A 10 pJ/m), "large" = A 10 pJ/m at 1 µm. Curves of the main
fit: geometric mean over the design points (stress and packing averaged),
95 % intervals by bootstrap over the runs.

**A low β does not mean a low loss density.**
- Four materials (corners of A and K_eff):

| material | β (9 … 150 mT) | P at 9 mT | P at 50 mT | P at 150 mT |
|---|---|---|---|---|
| small particles, K_eff 1.0 kJ/m³ | 2.55 | 2.0 | 137 | 2259 |
| small particles, K_eff 6.2 kJ/m³ | 2.30 | 3.0 | 140 | 1750 |
| large particles, K_eff 1.0 kJ/m³ | 1.61 | 11.9 | 201 | 1170 |
| large particles, K_eff 6.2 kJ/m³ | 2.15 | 4.2 | 174 | 2060 |

- β alone does not rank the materials: small particles with high K_eff
  (β 2.30) lose less than large particles with high K_eff (β 2.15) at all B.
- Simulations in thirds by β (geometric means): low-β third (β 1.82) 7.7 /
  193 / 1363 W/cm³, high-β third (β 2.38) 2.8 / 168 / 2005 W/cm³ at 9 / 50 /
  150 mT; ratio high / low × 0.36 / 0.87 / 1.47; crossing ≈ 59 mT.

**The best material depends on B_max** (lowest P on the A × K_eff grid):

| B_max | best A | best K_eff | β | P | worst / best | fine-mesh pair A 20 / A 10 (K3 / base) |
|---|---|---|---|---|---|---|
| 9 mT | 19.9 pJ/m | 1.0 kJ/m³ | 2.55 | 2.0 W/cm³ | × 6.0 | × 0.82 (n.s.) |
| 50 mT | 19.9 pJ/m | 2.2 kJ/m³ | 2.38 | 134 W/cm³ | × 1.7 | × 1.59 (hint, reversed) |
| 100 mT | 9.9 pJ/m | 1.5 kJ/m³ | 1.69 | 591 W/cm³ | × 1.6 | × 1.91 |
| 150 mT | 9.9 pJ/m | 1.0 kJ/m³ | 1.61 | 1170 W/cm³ | × 1.9 | not checked |

- The optimum lies at the edge of the studied range in every case.
- The mesh error acts also on the ratios between particle sizes: the
  switch point of the best material is ≈ 70 … 100 mT on the coarse mesh and
  below 50 mT for the only fine-mesh pair (which also differs in K_eff,
  stress and packing).

**Loss densities along the levers (W/cm³).**
- A 10 → 20 pJ/m (ends, others averaged): 7.5 → 2.4 at 9 mT, 1423 → 1874
  at 150 mT. Particle size 1 → 0.71 µm at K_eff 2.5 kJ/m³: 9.0 → 3.1 at
  10 mT, 609 → 725 at 100 mT.
- Packing 55 → 69 %: 3.9 → 5.3 at 9 mT (× 1.08 without the dilution), 1630
  → 1831 at 150 mT.
- Damping (centre, steady loops ≥ 2): α 0.1 → 0.01: 173 → 127 at 50 mT,
  1456 → 898 at 150 mT.

**Where the loss is** (median of all runs, per volume of particle
material): particle centre (r < R/2) 19.8 W/cm³ at 9 mT and 3831 W/cm³ at
150 mT; rim (r > 0.8 R) 3.3 and 2654 W/cm³; whole particle 7.0 and 2634
W/cm³. Centre / rim: × 6 at 9 mT, × 1.4 at 150 mT. The loss does not
spread evenly at high B (correction in 7.69).

**Loss angle by A** (thirds, median): high A sin δ 0.0036 → 0.0096 (9 →
150 mT, × 2.7); low A 0.0078 → 0.0068.

**Turn angle of the magnetisation** θ = arccos(𝐦(t₊) · 𝐦(t₋)) between the
positive and the negative field peak (𝐦 = unit vector of the
magnetisation).
- Cut of K1_2 (A 19 pJ/m, K_eff 1.0 kJ/m³): mean 2° / 12° / 32° at 9 / 50 /
  150 mT; 25 % of the cut turns by more than 30° at 150 mT.
- Mean over the particles with 𝐦 averaged over one Herzer cell (= 180° ×
  V_sw / V), median of all runs: 1.2° / 6.0° / 17.7°; low-A third 1.3° /
  6.0° / 16.5°, high-A third 1.0° / 5.9° / 19.6°.

**Per particle, stress axis** (runs with K_σ > 0.5 K_eff, median, per
particle volume): axis at 37° to the drive 7.4 / 3608 W/cm³, at 79° 5.5 /
1882 W/cm³ at 9 / 150 mT (× 1.34 / × 1.92). Hint (no run without stress).

**Real alloy within the design** (K_eff ∝ A⁻³ at a fixed grain size, from
A 10 pJ/m, K_eff 6.2 kJ/m³ to A 18.3 pJ/m, K_eff 1.0 kJ/m³): β 2.15 → 2.42,
Δβ +0.27 (+0.15 … +0.41). The value +0.40 for A × 2 (7.66) needs an
extrapolation beyond the design.

### 7.69 Where the loss is: a few zones at every B (2026-10-12, correction of 7.68 and §0)

Chris saw on slide 8 that the maps do not show "first in the vortex core,
at high B everywhere". Single zones dominate, and at 9 mT a hot zone sits
at the rim of a neighbour particle. Check with the cell loss of all 116
runs (acc.json, calm accumulator windows only, see 7.70; n 70 / 102 / 114
runs at 9 / 50 / 150 mT; coarse mesh; median and quartiles):

| share of the loss | volume | 9 mT | 50 mT | 150 mT |
|---|---|---|---|---|
| 10 % of the cells with the highest loss | 10 % | 69 % (63 … 72) | 76 % (72 … 78) | 73 % (70 … 74) |
| 10 % with the largest texture | 10 % | 61 % (55 … 64) | 69 % (65 … 72) | 55 % (51 … 60) |
| centre, r < R/2 | 12.5 % | 41 % (36 … 46) | 28 % (24 … 32) | 18 % (16 … 20) |
| rim, r ≥ 0.8 R | 48.8 % | 21 % (18 … 24) | 39 % (33 … 42) | 50 % (46 … 52) |

- Texture = |∇⟨𝐦⟩| of the loop-mean 𝐦, smoothed inside each particle
  with σ = L_eff / 2 (acc_reduce.py).
- Each run has ≥ 50 % of its loss in its hottest 10 % (minimum at 9 mT).
- The centre has a higher loss density than the rim in 96 … 100 % of the
  runs at each B (the hottest 10 % have a higher density by definition).
- Loss of the motion at f only (p_fund_share, what the snapshot maps
  show): rim 20 / 30 / 38 % at 9 / 50 / 150 mT against 21 / 39 / 50 % of
  the full loss. Thus the loss at other frequencies (7.70) lies more at the
  rim.
- Map simulation: K1_2 has the least concentrated loss of the runs with
  a calm window at 9 and 20 mT (50 % in its hottest 10 %, rank 1 of 70 and
  of 98). The deck now uses
  S041 (A 12 pJ/m, K_eff 2.8 kJ/m³, K_σ 0.9 K_eff, packing 66 %), close
  to the median curve of the concentration, on slides 7, 8 and 10.
- Cut of S041: 74 / 62 / 48 % of the 10 % cut cells with the highest q lie
  in the 10 % with the largest texture of the cut at 9 / 50 / 150 mT (10 %
  if there were no link). In K1_2 the values were 74 / 71 / 49 %. The hot
  zone of K1_2 at 9 mT, bottom left, lies on the line where the
  out-of-plane component of the loop-mean 𝐦 changes its sign. Its type
  (vortex core or wall) is not determined.
- Method of the maps: q = Σ |Δ𝐦|² over 16 snapshots. q agrees with |𝐦₁|²,
  the loss of the motion at f (correlation 0.93 … 1.00 in the cut of S041,
  0.95 … 1.00 for K1_2). In S041 this motion carries 92 / 46 / 45 % of the
  loss at 9 / 50 / 150 mT. The real cell-loss slices of the old runs
  overflowed in float16 (old format, fixed in 7.61). The runs K1A10
  (queue) give real loss maps.
- Colour scale of the maps: before this entry, 0 … 4 with white at 2 ×
  the mean, so cells with 1 … 2 × the mean showed blue and the text "blue
  = less than the mean" was wrong. Now white = the mean (TwoSlopeNorm).
- Correction: "at 150 mT the loss spreads over the whole particle" (§0,
  7.68) is wrong. 7.23 says that out of the textured regions the bulk
  rotation takes over. This reading is also corrected here: the loss
  stays in a few zones at every B. With rising B the weight moves from
  the centre to the rim. Cause of the error: the conclusion came from rim
  share = volume share, without the concentration number in the same
  paragraph.
- Deck: slide 8 new (maps with the border of the texture zones, shares
  against B with quartiles, terms defined, statements reduced to the
  numbers), slide 25 (conclusion) with the measured result in its own
  bullet. A context-free review (2026-10-12) found 15 points. All points
  on slides 8 and 25 and in this record are done. Its point on slide 10
  (turn angle) is done: one sentence links the two slides.

### 7.70 Loss channels: with rising B the loss moves to other frequencies (2026-10-12)

Question of Chris: does the loss start at the drive frequency f and go
more and more into harmonics? Go for options A + B + C (local, no cost).

**Method correction first (calm windows).** acc_reduce.py sums all
cycles from cycle 2 on (5 … 10 cycles per stage). The steady rule for P
uses calm cycles only. Windows with an event cycle (eval_events, 7.53):
46 / 18 / 25 / 14 / 4 / 7 / 2 of 116 at 9 / 20 / 35 / 50 / 70 / 100 /
150 mT. In these windows R1 is lower (median 0.71 against 0.91 at 9 mT).
- Fix (deck_data.py, acc_stages): all values from acc.json (R1, loss
  shares, region densities) use calm windows only. Individual cycles
  cannot be removed from the sums afterwards. Future runs should store the
  accumulators per cycle.
- R1, median, calm windows: 0.91 / 0.75 / 0.59 / 0.50 / 0.47 / 0.44 / 0.40
  (all windows: 0.86 / 0.72 / 0.56 / 0.50 / 0.47 / 0.44 / 0.40). The fall
  stays.
- jumps_block read acc.json of 13 reruns from the old folder (orig). It
  now reads the folder of the steady analysis.

**Channels.** Per run and amplitude (deck_data.py, channels_block):
- f: R1 = share of the cell loss from the local motion at f (calm window).
- jumps: A = energy of the box dissipation p_dis(t) above its smooth
  periodic curve Ps (mean cycle, Fourier k ≤ 6), last 3 calm cycles,
  256 samples per cycle. A = RF / 2 exactly, because P − Ps sums to zero.
  It counts bursts shorter than ≈ T / 13 ≈ 2.6 ns.
- other: 1 − R1 − A = smooth local motion at other frequencies (harmonics
  3f, 5f …, slow drift). Harmonics and slow drift are not separated.
- C = RN / 2 = energy that differs from cycle to cycle. C is part of A
  (P − Ps contains P − Pm). Also smooth motion at 5f and above gives
  components of p(t) above k = 6 and goes into A. A comes from the box
  sum, where bursts of different cells partly cancel: a lower limit of the
  local jump energy. "other" is the rest, so the 3 shares add to 1 by
  construction (minimum of "other" in a run: 0.014).
- B: energy share Eev and time share ev of the sample intervals in which
  n60 changes (neighbour pairs > 60°: vortex cores, Bloch points, steep
  walls).

Mean channel shares (runs with both values valid; n 70 / 102 / 114 at
9 / 50 / 150 mT):

| group | channel | 9 mT | 50 mT | 150 mT |
|---|---|---|---|---|
| all | f | 88 % | 52 % | 40 % |
| all | other | 8 % | 37 % | 52 % |
| all | jumps | 3 % | 12 % | 8 % |
| low A (9.7 … 12.5 pJ/m) | f / other / jumps | 86 / 10 / 4 % | 51 / 39 / 10 % | 40 / 53 / 7 % |
| high A (15.8 … 20.4 pJ/m) | f / other / jumps | 90 / 7 / 3 % | 52 / 35 / 13 % | 39 / 52 / 9 % |

Growth exponent of the loss density in each channel (fit 9 … 150 mT,
mean ± SE, runs with all 7 amplitudes valid):

| group | n | f | other | jumps | all (β) |
|---|---|---|---|---|---|
| all | 39 | 1.85 ± 0.03 | 2.89 ± 0.08 | 2.76 ± 0.11 | 2.14 ± 0.04 |
| low A | 11 | 1.68 ± 0.06 | 2.52 ± 0.13 | 2.22 ± 0.17 | 1.92 ± 0.08 |
| high A | 17 | 2.00 ± 0.03 | 3.22 ± 0.10 | 3.24 ± 0.13 | 2.33 ± 0.04 |

The β of these subsets agrees with the β of the full thirds (1.94 / 2.34).

- Jumps carry ≤ 12 % on average (≤ 13 % in the high-A third), single
  runs up to 23 %. C is ≤ 7 % on average.
- B: n60 changes in 17 … 80 % of the sample intervals (9 → 150 mT).
  These intervals carry 18 … 82 % of the loss, about the same as their
  time share (ratio 1.08 … 1.11 at 9 mT). B gives no useful limit: on the
  coarse mesh n60 changes also when a core moves smoothly by one cell.
- Spearman correlation of β with ln(R1(150) / R1(9)), calm windows:
  ρ = −0.62 (p 1·10⁻⁸, n 70). Partly by construction, because
  β = β_f − d ln R1 / d ln B (7.71).
- Mesh (K3, n 5 … 6 per mesh and B): R1 at 20 / 70 mT 0.78 / 0.48 (dx 3),
  0.70 / 0.57 (dx 1.5). Jumps (RF / 2) on the fine mesh 4 … 7 % at all B,
  on the coarse mesh up to 15 % at 50 … 70 mT. Thus the jump peak is
  mostly an effect of the mesh. The fall of R1 exists on both meshes and
  is smaller on the fine mesh.
- The source delivers power only at f (H is a sine). Inside the particles
  a growing share of this power dissipates through local motion at other
  frequencies. The statement "its rising share lifts β above 2" was too
  strong: see 7.71 (β = β_f + shift).
- Answer to Chris: yes in the direction, with "harmonics" → "motion at
  other frequencies (not separated)". Jumps are a small channel.
- Deck: new slide 9 (channel shares for the low-A and high-A thirds,
  exponent table). Slide 21 (regime) now shows the jump share as the
  energy share RF / 2 with the same name as slide 9 (before: RF, labelled
  "jump share of the dissipation"). Glossary: "loss channel".
- Open: slides 9 and 21 overlap (R1 and jumps). Decision for Chris: merge,
  or keep slide 21 as the mesh check. The direct per-cell split of the
  other channel into harmonics and slow drift needs per-cell accumulators
  at 3f, 5f … (option D, not started).

### 7.71 Loss channels: local checks O1, second review, pause before C′ (2026-10-12)

Chris asked how secure the channel split of 7.70 is. Go for O1 (local
checks now), and for O2 (extended per-cell accumulators) before C′, with a
pause of the C′ chain. Script: eval_channels.py → results/loss_channels_o1.md.

**Identity β = β_f + shift.** P_f = R1 · P, thus d ln P / d ln B =
d ln P_f / d ln B − d ln R1 / d ln B exactly (check: residual 10⁻¹⁵).
Runs with all 7 amplitudes calm, mean ± SE:

| group | n | β | β_f | shift |
|---|---|---|---|---|
| all | 39 | 2.14 ± 0.04 | 1.85 ± 0.03 | 0.29 ± 0.01 |
| low A | 11 | 1.92 ± 0.08 | 1.68 ± 0.06 | 0.25 ± 0.02 |
| high A | 17 | 2.33 ± 0.04 | 2.00 ± 0.03 | 0.33 ± 0.01 |
| high − low (Welch SE) | | +0.41 ± 0.09 | +0.33 ± 0.07 | +0.08 ± 0.03 |

- The fall of the share at f adds 0.25 … 0.33 to β in all materials.
- The A effect on β sits mostly (≈ 80 %) in β_f, the growth of the loss
  of the motion at f. High A: β_f ≈ 2.0, a constant loss angle of the
  local motion at f. Low A: β_f ≈ 1.7, this loss angle falls with B.
- Factor effects (OLS, HC3, range effect, n 39): A on β_f +0.42 ± 0.08,
  on the shift +0.10 ± 0.03. K_eff on β_f +0.14 ± 0.06, on the shift
  +0.08 ± 0.03. The A × K_eff term acts on β_f (+0.41 ± 0.11).

**Where the rest 1 − R1 sits** (share of the non-f loss per zone, median,
calm windows; in brackets the share of the loss at f):
- 10 % cells with the highest loss: 73 … 88 % (55 … 68 %).
- 10 % with the largest texture: 63 … 84 % (43 … 60 %).
- Rim (49 % of the volume): 28 % at 9 mT → 57 % at 150 mT (20 → 38 %).
- Thus the loss at other frequencies is more concentrated than the loss at
  f, in the zones of non-uniform 𝐦 and with rising B at the rim.

**Damping** (centre C0 … C7, same microstructures, last 3 cycles):

| B (mT) | 10 | 21 | 51 | 101 | 151 |
|---|---|---|---|---|---|
| R1, α 0.1 (median, calm n 6 … 8) | 0.94 | 0.74 | 0.53 | 0.46 | 0.43 |
| R1, α 0.01 (median of 8, all windows) | 0.17 | 0.10 | 0.11 | 0.09 | 0.08 |
| jumps A, α 0.1 / α 0.01 | 0.012 / 0.13 | 0.08 / 0.12 | 0.11 / 0.12 | 0.09 / 0.09 | 0.07 / 0.09 |

- At α 0.01 only 0 … 4 of 8 windows are steady (calm R1 0.09 … 0.12,
  n 1 … 3). Jumps and loop-to-loop changes are ≈ 10 %. Thus ≈ 80 % of the
  loss at α 0.01 is smooth motion at other frequencies, already at 20 mT.
- Hypothesis: at small damping the local motion rings at its own GHz
  frequencies (precession), and this energy is lost to the source whatever
  α is. This agrees with P falling only × 0.6 … 0.7 for α × 0.1 (7.68).
  Not shown: the data cannot separate harmonics of f from precession.

**Second context-free review** (21 points) of slides 7 … 10, 21, 25, 27
and of 7.69 / 7.70. Done:
- Slide 9: the statement "the rising share of other frequencies lifts β
  above 2" replaced by the identity table (low A, high A, difference).
  Jumps defined as "energy of the total loss above its smooth periodic
  curve: fast bursts and loop-to-loop changes, a lower limit". Other
  frequencies "harmonics, local precession, slow drift: not separated".
  Cycle sets, n per B, single-run maximum (23 %), mesh effect on the jump
  peak, α 0.01 limit stated. Table without the red highlight (red = jumps).
  n60 and Spearman bullets removed from the slide.
- Slide 8: subtitle "all samples of loops 2 … 6, windows with an event
  left out". "The pattern of q agrees with the loss at f" (q also contains
  3f and 5f).
- Slide 21: "share of the loss at f" (axis), (a) rest = other frequencies
  and jumps, (b) the same jump definition as slide 9.
- Slide 7: caption for S041 (core in the centre particle, corner particles
  with 𝐦 mostly out of the cut plane).
- Slide 25: measured result with α 0.1 and the identity. Slide 2: new row
  "Loss mechanism" (the cells do not respond linearly; depends on the
  damping). Glossary: "loss channel" with the jump definition.
- Not done: slide 8 title keeps "mostly where 𝐦 is non-uniform" (55 % at
  150 mT is still the larger part; the numbers are on the slide).

**Pause before C′.** /root/chain_freq.sh now waits for runs/PAIR_DONE and
for runs/FREQ_GO (instance/chain_freq.sh). The old waiting process was
stopped by its PID. Its child sleep held the lock and was also stopped
by its PID. The new process (20:57 UTC) waits. FREQ_GO is written only
after O2 is tested and /root/freq_commit points to the O2 commit.

### 7.72 Option O2: per-cell harmonics in the accumulators (2026-10-12, method fixed before the result)

Go of Chris: pause before C′, then O2 first, then C′ with O2. Question of
Chris: "das wäre ja kein Risiko, oder?". Answer: three small risks, each
covered before FREQ_GO: an error in the new accumulators (it would cost
the 12 C′ runs, ≈ 3.5 $) → test run with fixed criteria; disk → only
|m_k|² per cell; GPU memory → all buffers at the start of the job.

**Method** (run_loops.py --acc --acc_harm K, class HarmAcc; acc_reduce.py
harmonic_loss; commit dc95322):
- Per cycle and cell: sums of p, m, m e^{−iωt} (float64, as before) and
  m e^{−ikωt} for k = 2 … K (float32), over the 256 samples of the cycle.
- At the end of each cycle: |m_k|² = |2/n Σ m e^{−ikωt}|², summed over
  the 3 components, k = 1 … K. Box sums over the particle cells go into the
  cycle record (cyc["acc"]: calm, n, p_cells, mk2_cells).
- Only a calm cycle (w_loop > 0, |w_dis / w_loop − 1| ≤ 0.20, the same
  rule as min_calm and the event rule of 7.53) enters the stage sums.
  This is "calm cycles", not the steady-tail rule: a calm cycle before a
  later event is kept. A stage with no calm cycle writes no acc file.
- Saved per stage: p, m0, m1 (coherent over the kept cycles, as before),
  mk2 = mean of |m_k|² over the kept cycles (float32, cells × K),
  cycles_kept, cycles_acc.
- p_k = α μ0 Ms / γ (kω)² |m_k|² / 2 per cell: the same coefficient as R1
  (Gilbert damping: p = α μ0 Ms / γ |dm/dt|²). R_k = Σ p_k / Σ p, rest
  R_rest = 1 − Σ_{k ≤ K} R_k = motion above Kω plus motion that does not
  repeat within a cycle. Zone shares of the harmonic sum (k ≥ 2) and of the
  rest; slices relative to their plane mean (7.61).
- Without --acc_harm the code path is the old one (unchanged lines).
- Unit tests (test_study.py): known harmonics are recovered without leak
  into other k, an event cycle is not added, Parseval: Σ p_k equals the
  mean Gilbert loss of a motion with harmonics ≤ K (to 10⁻⁵). All 60 study
  tests pass (CPU).

**Reading rule (fixed now):**
- Large rest → fast or non-repeating local motion (for example precession
  at GHz frequencies; it lies above the Nyquist frequency 128 f and is
  folded: O2 cannot give its frequency).
- Small rest → the "other frequencies" of 7.70 are low harmonics (k ≤ 7).
- C′ (f/2, 2 f): if R_k depend on B only, the response is a static
  nonlinearity; if they move with f, a rate effect.

**Test on the instance before FREQ_GO** (production job line of K1_0_fhalf,
own worktree /root/o2test at dc95322, outputs /root/o2test_runs, one run
with and one without --acc_harm 7, start from runs/K1_0/init.pt, 9 and
50 mT, 5 cycles each). Criteria:
1. exit code 0 for both runs;
2. the samples files of both runs are equal (the accumulators only read);
3. p_box / p_dis = 1.000 (CHECK_TOL 10⁻³) over the kept cycles;
4. when all accumulated cycles are calm: p, m0, m1 and R1 of both runs
   agree (relative 10⁻⁵);
5. Bessel: Σ_k R_k ≤ 1, and the share of cells with Σ_k p_k > 1.001 p
   below 1 % (a cycle is not exactly periodic, so a small excess in single
   cells is possible);
6. GPU memory and wall time per cycle against the run without harmonics;
7. acc_reduce.py and the post-processing run without error on the test
   folder.
If one criterion fails: no FREQ_GO, report.

**Production change:** pair_jobs.py freq_lines adds --acc_harm 7; the
regenerated jobs/freq_K1pair.txt is equal to the old one apart from this
option (checked).

**Test result (2026-10-12, 22:20 UTC, o2_check.py):**
- Pass: 1 (exit 0), 2 (samples files equal byte for byte), 3 (p_box / p_dis
  = 1.000000), 4 (p, m0 equal, m1 to 4·10⁻¹⁴, R1 0.898664 / 0.588929 in
  both runs), 6 (wall per cycle × 0.94 / 0.97 against the run without
  harmonics; GPU memory + ≈ 0.3 GB), 7 (acc_reduce clean).
- Fail: 5 at 9 mT. Cells with Σ_k p_k > 1.001 p: 8.4 % (limit 1 %), max
  ratio 2.22. At 50 mT 0.01 % (max 1.002): pass.
- Cause analysis: the run without O2 (old R1 method, coherent m1) shows the
  same: 4.4 % of the cells with p_f > 1.001 p at 9 mT, max 2.21, none at
  50 mT. The excess is not in the harmonics (p_2 … p_7 = 0.2 … 0.3 % of p in
  these cells) but in k = 1. The cells lie mostly in the vortex-core region
  (33 % of them in r < R/2, 12 % of the volume). Energetically the excess
  is Σ max(p_f − p, 0) / Σ p = 0.0003: R1 at 9 mT is biased by ≤ 0.0003.
  The energy balance of the calm cycles at 9 mT is w_dis / w_loop = 1.08 …
  1.10 (at 50 mT 1.00): the state still releases stored energy after the
  start. Hypothesis: fast local motion between the samples (256 per cycle)
  that the sampled p misses. Not shown.
- Values at 50 mT (K1_0 at f/2, one run): R1 0.589, R_2 … R_7 0.027 / 0.103
  / 0.028 / 0.026 / 0.012 / 0.015 (sum 0.21), rest 0.20.
- Decision Chris (option A, recommendation): O2 is released although
  criterion 5 failed, because the failure is not in the new code (the same
  in the old method) and is energetically 0.03 %. At 9 mT the per-cell
  values in the core region are read with this caveat. In parallel option
  B: a test at 9 mT with 1024 samples per cycle (same job line otherwise).
  Reading rule: if the share of cells with p_f > 1.001 p falls clearly
  (below 1 %), the cause is fast motion between the samples; if not, it is
  the slow release of energy after the start.

### 7.73 Decisions of Chris: O2 also at f, slide 21 kept, check of all slides (2026-10-12)

1. The pair runs K1A10_0 … 2 (at f) also get --acc_harm 7, and 3 runs
   K1_0 … 2 at f with --acc_harm 7 go into the C′ queue as K1_k_f1 (the
   same job lines as K1_k with the steady rule, no snapshots). Then the
   split into f, harmonics and rest exists at f/2, f and 2 f for both A.
   Cost estimate: ≈ 1 $ more. pair_jobs.py writes both job files; the
   old lines are unchanged apart from --acc_harm 7 (checked by diff). The
   C′ chain post-processes the names K1(A10)_k_f(half|1|2). The pair chain
   reads /root/pair_commit only after ALPHA_DONE: both commit files are set
   after the O2 test (7.72) has passed.
2. Slide 21 (regime) stays as the mesh check of slide 9.
3. All slides are checked against the rule for reading guides (1 … 2
   sentences per plot, no explanation in titles, axis labels and legends).
   A list of changes goes to Chris before any change.

### 7.74 Mesh check at d/l_ex 212, base microstructure: result (2026-10-12, PLAN 4.6)

Chain M212 finished 2026-10-12 22:27 UTC: 24 of 24 jobs, none failed, all
stages steady (≥ 2 calm cycles). Results fetched without *.pt (120 files,
7.2 MB, md5 equal to the instance; runs/MANIFEST_m212.md5, local).
eval_size.py → results/mesh212.md (rules of PLAN 4.6 B, fixed before the
result). Base microstructure: L_eff 12 l_ex (K_eff 6.2 kJ/m³), r_p 0,
φ 0.65, α 0.1, A 10 pJ/m. n 4 seeds per size and mesh; ± SE, 90 % CI (t).

Mesh error P(dx 3) / P(dx 1.5), paired per seed:

| size | 9 mT | 50 mT | 100 mT | Δβ (9 / 50 / 100 mT) |
|---|---|---|---|---|
| d/l_ex 212 (0.71 µm) | × 0.43 (CI 0.28 … 0.67) | × 0.91 (0.68 … 1.22) | × 1.54 (1.06 … 2.23) | +0.51 ± 0.13 (CI +0.20 … +0.81) |
| d/l_ex 300 (1.0 µm) | × 1.48 (0.92 … 2.36) | × 1.40 (0.97 … 2.03) | × 1.25 (1.07 … 1.45) | −0.06 ± 0.05 (−0.19 … +0.07) |
| K3 (d/l_ex 212, other microstructure) | × 0.63 | × 0.70 | × 0.87 | +0.12 ± 0.03 |

Size lever 212 → 300 (value at 1 µm minus value at 0.71 µm, Welch):

| output | fine mesh dx 1.5 | coarse mesh dx 3 | mesh share |
|---|---|---|---|
| β (9 / 50 / 100 mT) | −0.09 ± 0.12 (CI −0.35 … +0.17) | −0.65 ± 0.13 (−0.91 … −0.40) | −0.56 ± 0.14 (−0.86 … −0.26) |
| ln P 9 mT | −0.23 ± 0.22 (−0.69 … +0.23) | +1.00 ± 0.30 (+0.34 … +1.67) | +1.24 ± 0.27 (+0.70 … +1.77) |
| ln P 50 mT | −0.53 ± 0.09 | −0.09 ± 0.17 | +0.44 ± 0.20 |
| ln P 100 mT | −0.38 ± 0.12 | −0.59 ± 0.11 | −0.21 ± 0.17 (n.s.) |

Facts (no decision; PLAN 4.6: Chris decides):
- At this microstructure the coarse mesh makes the size effect on β
  larger: smaller particles raise β by 0.65 on dx 3, by 0.09 ± 0.12 (not
  different from 0) on dx 1.5. The mesh share −0.56 has a 90 % CI that
  excludes 0.
- At 9 mT the coarse mesh gives the small particles × 0.43 of the fine-mesh
  loss and the large ones × 1.48. On dx 3 the large particles lose × 2.7
  more at 9 mT (ln +1.00), on dx 1.5 about the same (ln −0.23 ± 0.22).
- The mesh error at 212 depends on the microstructure: Δβ +0.51 here,
  +0.12 at K3 (L_eff 18, r_p 1.5, φ 0.69).
- Limits (PLAN 4.6 C): one microstructure at high K_eff and without stress;
  the gap (3.1 cells at 212, 4.4 at 300 on dx 3) changes together with the
  size; dx 1.5 is not shown to be converged; 4 seeds.
- Post hoc (request of Chris, not pre-registered; eval_size_no9.py →
  results/mesh_no9.md): β from 50 → 100 mT only. Δβ (coarse − fine, paired):
  base 0.71 µm +0.73 ± 0.20 (CI +0.27 … +1.19, n 4), base 1 µm −0.10 ± 0.11
  (n.s.), K3 0.71 µm +0.36 ± 0.11 (n 6; over 20 … 100 mT +0.05 ± 0.07,
  n.s.). Size lever 50 → 100 mT: fine +0.15 ± 0.13 (n.s.), coarse
  −0.68 ± 0.16, mesh share −0.83 ± 0.23 (CI −1.29 … −0.37). Thus the mesh
  effect at the base point does not come from the 9 mT point alone: at
  0.71 µm the coarse mesh also gives too much loss at 100 mT (× 1.54). My
  first statement to Chris ("almost all at 9 mT") was wrong. Pattern over
  the 3 points: the β error sits at 0.71 µm between 50 and 100 mT (K3 and
  base), not at 1 µm. Hypothesis (not shown): jumps that the coarse grid
  makes at high B (jump peak only on dx 3, slide 21; 21 … 32 neighbour pairs
  > 60° on dx 3 against 2 on dx 1.5 at 100 mT).
- Statements that this touches (coarse mesh): the A / size lever on β (main
  result +0.57; fine-mesh estimate "≈ 67 %" from the confounded pair K3 /
  base, 7.63), "at 9 mT small particles lose least" and the switch of the
  best material with B (deck slides 2, 11 … 13, 16 … 19, 25). The main fit
  gives the smallest A effect at high K_eff (+0.16 at 6.2 kJ/m³, averaged
  over stress and packing); the large effect at low K_eff (+0.93) has no
  fine-mesh check.

### 7.75 C′ with the same sample interval at all frequencies (2026-10-12, decision of Chris)

With 256 samples per cycle the sample interval is 0.26 / 0.13 / 0.065 ns at
f/2 / f / 2 f (Nyquist 1.9 / 3.8 / 7.7 GHz). Fast local precession would be
caught differently at the three frequencies and could fake a frequency
dependence of the rest channel of O2. Decision of Chris: more samples and the
same interval T(f) / 512: 1024 samples at f/2, 512 at f, 256 at 2 f
(pair_jobs.py; the regenerated jobs/freq_K1pair.txt differs from the old one
only in --samples, checked by diff). The pair runs keep 256 (comparison with
the K1 runs of the design).
- eval_jumps_time.series reads the samples per cycle from the file, and
  deck_data.acc_stages from config.json (both assumed 256). New test
  test_jump_measures_any_samples_per_cycle.
- The test test_pair_and_freq_jobs still expected the 12 jobs before 7.73
  (I did not run the full test file at commit 439c878). It now checks 15
  jobs, O2 and the samples. All 61 study tests pass.
- Test B (7.72: 9 mT, f/2, 1024 samples) is also the production test of this
  setting. /root/freq_commit is changed only after it passes.
- Test B result (2026-10-12, 23:20 UTC; K1_0 at f/2, 9 mT, the same start
  state as the 256-sample run): exit 0, p_box / p_dis = 1.000000, cycles
  2 … 4 calm. Cells with p_f > 1.001 p: 3.1 % (256 samples: 4.4 %); with
  Σ_k p_k > 1.001 p: 4.8 % (8.4 %); max ratio 1.18 (2.22). R1 0.910 (0.899),
  R_2 … R_7 unchanged (0.006 / 0.015 / 0.001 / 0.001 / 0.000 / 0.000).
  w_dis / w_loop of the calm cycles 1.07 … 1.08 (1.08 … 1.10). Steps per
  cycle + 7 %, wall per cycle + 27 % (the GPU was shared with another job in
  both runs).
- Reading (rule of 7.72): the excess falls but stays above 1 %. Thus the
  sampling explains a part (the extreme cells: max ratio 2.22 → 1.18; R1 at
  256 samples is ≈ 0.01 too low at 9 mT), and a part remains, which agrees
  with the slow release of energy after the start (w_dis / w_loop > 1).
- Production: the wall time rises by less than 30 % (the limit given to
  Chris), so /root/freq_commit = 67cdc74 (run code equal to the tested
  4bb1900; changed are the samples of the C′ jobs and the analysis scripts).
  /root/pair_commit stays 4bb1900. FREQ_GO is set.

### 7.76 Decisions of Chris (2026-10-12, night) and a finding of the 30-loop test

- Mesh question (7.74): decision tomorrow.
- The instance keeps running after C′.
- The test folders /root/o2test (git worktree at dc95322) and
  /root/o2test_runs were deleted (go of Chris). Before that, the small
  files (JSON, logs, samples files, acc_slices; no *.pt, no raw accumulators)
  were saved to runs/o2test/ (24 files, 2.5 MB, md5 equal to the instance,
  runs/o2test/MANIFEST.md5, local).
- The pair runs K1A10_0 … 2 get 512 samples per cycle (decision of Chris):
  in C′ they are the f point of K1A10, and C′ uses the interval T(f) / 512
  (7.75). The field images keep their 16 phases. The partner runs K1_0 … 2
  of the design have 256 samples; test B showed small changes with 4 × the
  samples (R1 + 0.011 at 9 mT). pair_jobs.py; jobs/pair_K1A10.txt differs
  from the old one only in --samples (diff). 61 study tests pass.
- Finding (record only, decision of Chris): in the 30-loop test of the
  centre point (results/conditioning.md, 6 repeats at 9 / 20 / 35 mT), 3 of
  18 lines have a loss that repeats only every few loops, without an event:
  20 mT repeat 2: every 3rd loop + 5.8 % (exactly the same value); 35 mT
  repeat 4: every 4th loop + 9.0 %; 35 mT repeat 0: about every 5th loop
  + 7.5 … 9.1 %. The other 15 lines stay within ± 5 % of their late median
  from loop 6 on. A value from 2 calm loops can then be off by a few percent
  (within the drift limit ± 5 % of the conditioning test). To look at again
  after C′ with the per-cycle data of O2.

### 7.77 Damping test at small B: result (2026-10-13, PLAN 4.7)

Chain α finished 2026-10-13 04:47 UTC: 16 of 16 jobs, none failed; the
determinism gate (check_reruns.py) passed (exit code 0). Results fetched
without *.pt (136 files, 32 MB, md5 equal to the instance;
runs/MANIFEST_cond_a001.md5, local). eval_alpha_cond.py →
results/alpha_cond.md (rules of PLAN 4.7, fixed before the runs). Note: the
first version of the script did not set the steady rule for the reference
values; corrected before any use.

| B | seeds with a calm tail (≥ 2 of 30 cycles) | P(α 0.01) / P(α 0.1) | 90 % CI | SD of ln |
|---|---|---|---|---|
| 9 mT | 4 of 8 (C2, C3, C5, C7) | × 0.13 | × 0.11 … 0.16 | 0.14 |
| 20 mT | 3 of 8 (C0, C3, C4) | × 1.06 | × 0.38 … 2.91 (not determined) | 0.60 |

- At α 0.01 most runs do not settle at 9 and 20 mT even in 30 cycles:
  event cycles continue up to cycle 29 in several seeds.
- At 9 mT the 4 settled seeds give × 0.11 … 0.15 each: the loss falls
  almost in proportion to α (× 0.10). Above 35 mT (design runs, 7.68):
  × 0.62 … 0.81. Possible selection: the settled seeds can be the calmer
  ones.
- R1 at α 0.01 (acc windows, mostly with events): median 0.28 at 9 mT
  (range 0.003 … 0.85), 0.12 at 20 mT.
- Consequence for β: if the loss at small B scales with α and at high B
  only weakly, β is higher at α 0.01. Not measured: only C2 has steady
  values at both 9 and 150 mT (Δβ +0.65, n 1). No claim.
- Deck: slide 23 with the 30-cycle values at 9 and 20 mT (deck_data
  alpha_block, source per amplitude), text and callout without a β number.

### 7.78 Grid-pinning check with the data at hand (2026-10-13, post hoc)

Hypothesis of Chris: at 9 mT the core moves less than one cell, so the grid
holds it. eval_pinning.py → results/pinning_check.md (mesh pairs base 0.71 /
1.0 µm and K3, last 2 calm cycles, paired per seed). x = core shift estimate
ΔM_p / Ms · R (rigid vortex, order of magnitude only).
- x at 9 mT: 2.8 … 4.1 nm = 0.28 … 0.41 coarse cells (also < 1 fine cell of
  5 nm); 20 mT 0.6 cells; 35 mT 1.1; 50 mT 1.6 … 2.3; 100 mT 3.2 … 4.5.
- ΔM_p (particle magnetisation amplitude) coarse / fine = 1.00 … 1.01 at all
  points and B: the mesh changes the loss, not the reversible response
  (which the demagnetisation of the powder largely sets).
- sin δ coarse / fine: base 0.71 µm × 0.45 / 0.97 / 1.68 (9 / 50 / 100 mT);
  base 1.0 µm × 1.48 / 1.49 / 1.31; K3 × 0.67 / 0.81 / 0.64 / 0.71 / 0.68 /
  0.91 (9 … 100 mT).
- Jump share (RF / 2) fine / coarse at 9 mT: 0.064 / 0.020 (base 0.71),
  0.057 / 0.083 (base 1.0), 0.056 / 0.019 (K3); at 70 … 100 mT the coarse
  mesh has more: 0.09 … 0.14 against 0.04 … 0.06.
- Reading: for the hypothesis: at 9 mT (x < 1 cell) the coarse mesh has
  less loss and fewer bursts at 0.71 µm; at high B more bursts. Against it:
  at 1 µm the coarse mesh has more loss at 9 mT (x also < 1 cell), and the
  mesh error stays at x = 2 … 4 cells (base 1 µm × 1.3 … 1.5, K3 × 0.64 …
  0.71). Thus "shift below one cell" does not explain the pattern in
  general. Sign and size of the mesh error depend on the microstructure at
  every B.

### 7.79 Damping α 0.01 against α 0.1: β fits with and without 9 mT (2026-10-13, request of Chris)

eval_alpha_fits.py → results/alpha_fits.md / .png; deck_data alpha_fits_block; deck slide 23 (means with 90 %
intervals and the number of settled runs, no single seeds). Centre point, 8 runs per B.
- Variant converged: stages with ≥ 2 calm cycles (9 and 20 mT from the 30-cycle runs, 35 … 150 mT from the 9-cycle
  design runs). Variant all: every run, time mean of w_dis of the last cycles (10 at 9 / 20 mT, 3 above).
- Runs that settle at α 0.01, of 8, per B (9 / 20 / 35 / 50 / 70 / 100 / 150 mT): 4 / 3 / 5 / 5 / 6 / 6 / 4. At α 0.1 all 8 settle at
  every B.

Loss ratio P(α 0.01) / P(α 0.1), geometric mean (90 % CI), per B 9 … 150 mT:
- converged: × 0.13 (0.11 … 0.16) | × 1.06 (0.38 … 2.91) | × 0.80 (0.44 … 1.44) | × 0.73 (0.50 … 1.07) | × 0.75 (0.54 … 1.05) | × 0.81 (0.69 … 0.95) | × 0.62 (0.49 … 0.77)
- all: × 0.23 (0.15 … 0.37) | × 0.60 (0.35 … 1.04) | × 0.85 (0.59 … 1.22) | × 0.84 (0.65 … 1.09) | × 0.69 (0.51 … 0.93) | × 0.79 (0.66 … 0.93) | × 0.61 (0.55 … 0.69)

Δβ = β(α 0.01) − β(α 0.1), one offset per run (± SE, 90 % CI):

| fit | variant | Δβ |
|---|---|---|
| with 9 mT: 9 … 150 | converged | +0.47 ± 0.13 (+0.25 … +0.70), 33 values |
| without 9 mT: 20 … 150 | converged | -0.17 ± 0.11 (-0.37 … +0.03), 29 values |
| 9 / 20 / 35 | converged | +1.75 ± 0.52 (+0.53 … +2.98), 12 values |
| 20 / 35 / 50 | converged | -0.10 ± 0.25 (-0.61 … +0.41), 13 values |
| 35 / 50 / 70 | converged | -0.32 ± 0.24 (-0.77 … +0.14), 16 values |
| 50 / 70 / 100 | converged | +0.14 ± 0.32 (-0.46 … +0.73), 17 values |
| 70 / 100 / 150 | converged | -0.17 ± 0.32 (-0.78 … +0.45), 16 values |
| with 9 mT: 9 … 150 | all | +0.32 ± 0.09 (+0.17 … +0.46), 56 values |
| without 9 mT: 20 … 150 | all | -0.01 ± 0.10 (-0.17 … +0.16), 48 values |
| 9 / 20 / 35 | all | +0.97 ± 0.23 (+0.57 … +1.38), 24 values |
| 20 / 35 / 50 | all | +0.40 ± 0.24 (-0.02 … +0.81), 24 values |
| 35 / 50 / 70 | all | -0.30 ± 0.17 (-0.61 … +0.00), 24 values |
| 50 / 70 / 100 | all | -0.09 ± 0.25 (-0.54 … +0.35), 24 values |
| 70 / 100 / 150 | all | -0.16 ± 0.21 (-0.53 … +0.20), 24 values |

- Both variants: without 9 mT no change of β (CI includes 0); the effect sits in the window 9 / 20 / 35 mT.
- Reading: at 9 mT the loss is mostly motion at f (R1 0.91 at α 0.1), which scales with α; above 35 mT motion at
  other frequencies dominates and scales only weakly. Thus α acts on β only through the low-B end. Centre point,
  coarse mesh; whether the coarse mesh gives the same ratio as the fine mesh is open.
- The SE assume independent values within a run; they are rather too small.
- Independent check of the 9 mT value (context-free agent, 2026-10-13): arithmetic and steady tails confirmed
  (w_dis SD 0.8 … 2.2 %, no drift); one settled run is in another state after large events; the runs that do
  not settle give a higher time-mean loss (variant all); one run is a steady period-2 state, which the single-cycle
  rule rejects (open method point, as the 3 … 5 cycle patterns of 7.76).

### 7.80 Memory and time of fine meshes; audit of all influences (2026-10-13)

**Memory and time test** (go of Chris; GPU 3 of the instance, base microstructure, short relaxation of 1000
iterations, 1 cycle at 9 mT, outputs in /root/convtest, outside runs/):

| case | cells | GPU memory | relaxation per iteration | one cycle |
|---|---|---|---|---|
| d/l_ex 212, dx 3 (production) | 108³ = 1.3 M | ≈ 1 GB | 0.04 s | ≈ 85 s (at f) |
| d/l_ex 212, dx 1.5 (production) | 2.1 M particle cells | ≈ 3 GB | 0.40 s | ≈ 0.4 h (estimate from the jobs) |
| d/l_ex 300, dx 1.5 (production) | – | – | 1.59 s | – |
| d/l_ex 212, dx 1 (test) | 312³ = 30 M | 8.4 GB | 4.7 s | not reached (stopped at 07:40 for the next case) |
| d/l_ex 96, dx 0.75 (test) | 192³ = 7.1 M | 2.3 GB | 1.07 s | 1.28 h (59 300 steps, dt 0.56 ps) |

Production relaxations need 11 000 … 25 500 iterations (M212, M300). Thus a full relaxation takes ≈ 3 … 8 h at
d/l_ex 96, dx 0.75, and ≈ 14 … 33 h at d/l_ex 212, dx 1 (estimates). A run with 7 amplitudes has ≈ 50 … 60 cycles.

**Audit of all influences** (context-free agent; numbers of the agent, the frequency number recomputed by me):
- The mesh is not the only open influence of high weight. Frequency: β(2f) − β(f) over 20 … 150 mT from the old
  frequency test 7.9 = +0.35 ± 0.03 (base, n 4) and +0.52 ± 0.06 (corner L_eff 30, r_p 3, old stress axes, n 4)
  (steady rule there: ≥ 1 calm cycle). Damping: ≈ 80 % of the A lever sits in the loss at f (7.71), which carries
  only ≈ 10 % of the loss at α 0.01; A × α was never tested.
- Medium weight: box of 4 particles, fixed stress axes and drive direction; the form of the Herzer cells; T = 0;
  the initial state at 1 µm; solver tolerance and float32 under production conditions.
- Checked and robust: event threshold 0.10 / 0.20 / 0.30 (A lever +0.57 / +0.57 / +0.59); eddy currents; determinism;
  relaxation rest torque; w_loop against w_dis. Without 9 mT the A lever is +0.69 (20 … 150 mT), in the range where
  the mesh error is largest. "A few zones" holds on both meshes (K3); the shift from the centre to the rim is partly a
  mesh effect.

### 7.81 B0 prepared: state transfer between meshes, jobs, analysis rules (2026-10-13, go of Chris for the preparation)

Goal: step 1 of variant B of PLAN 4.10. All runs start from the same physical state, so that the paired mesh
difference measures the discretisation and not the state selection. Details and the pre-registered reading: PLAN 4.10
"B0 as prepared". Not started; the run needs a separate go.

- Code: `run_loops.py --init_transfer <fine run>/init.pt` (function `transfer_m`): block average of m over the fine
  sub-cells of the same particle (or void), normalised; then the production relaxation. Checks: same INIT_KEYS except
  dx_lex and relax_maxiter, converged source relaxation, voxel surface, dx of the source divides dx and a.
- Unit tests (4 new, pass): a constant field stays constant; a block-constant field is recovered exactly also with
  other vectors in the void sub-cells; norm 1; the fallback is used and counted; a non-integer cell ratio is rejected.
  On the B0 grids (d/l_ex 96, 48 / 96 / 192 cells) every coarse particle cell has a fine sub-cell of the same particle
  (0 fallbacks, 0 zero vectors). The cube axes agree on dx 1.5 and dx 0.75 (extended test).
- Run test on the CPU (d/l_ex 18, φ 0.35 for a box that fits on dx 3; code path only, not physics): dx 1.5 → dx 3
  transfer, the relaxation converges, the transfer record is in summary.json and init.pt. A source with another seed
  and a source with another relax_tol_E are rejected (exit code 1).
- Analysis `eval_conv.py` with a synthetic test (known answers: go with transfer, stop, go with small scatter, not
  determined; the paired CI equals the scipy t interval). Jobs `conv_jobs.py` → `jobs/conv_b0.txt` (20 jobs, test: only
  the planned arguments differ from the mesh212 template). Chain `instance/chain_conv.sh` (not installed).
- Cycle noise at 9 mT (earlier mesh runs, pooled SD of ln w over the calm cycles, 4 runs per group): M212 dx 1.5
  0.058, M212 dx 3 0.063, M300 dx 1.5 0.077, M300 dx 3 0.029 (at 50 / 100 mT: 0.087 / 0.028). The calm tails have
  2 … 3 cycles. Thus a paired difference of two tail means has a noise of ≈ 0.024 … 0.077 in ln P.
- The measured relaxation of the dx 0.75 test (7.80) stopped at its cap of 1000 iterations (not converged). Thus the
  iteration count of a full relaxation on dx 0.75 is still an estimate (11 000 … 25 500 from the production meshes).
- Context-free review (one agent, ≈ 210 000 tokens). No problem in: the transfer arithmetic, the cube axes on all three
  meshes (also for the seeds 941 … 944), T against I (only the start state differs), the job lines, the statistics, the
  units. Findings, all fixed after the go of Chris: (1) the queue would start an unconverged dx 0.75 relaxation again
  with 200 000 iterations (≈ 59 h, up to ≈ 39 $ in total, estimate) → `--gate_maxiter 0`; (2) the reading had no rule
  for fewer than 3 independent pairs → new result "go, benefit not determined"; (3) the steady rule was set only in
  main() → set inside run_record; (4) after a reboot the B0 chain could take the lock of the C′ chain before its
  post-processing → the B0 chain waits for the end of the C′ process and takes no foreign lock; (5) wrong noise range
  and missing n 3 error rates in PLAN → corrected; (6) contradictory and too low cost numbers (dx 3 at d/l_ex 96 is
  44 s per cycle, not 6 s) → corrected; (7) the cleanup did not protect the F states → it reads --init_transfer;
  (8) b_peak / b_target per run added to the report; status names, file path, tests with real run folders.
- Decision of Chris 2026-10-13: option A for the go rule (threshold 0.10, 4 seeds, ≥ 2 calm cycles).

### 7.82 Pair K1 / K1A10: the A lever measured directly at one corner (2026-10-13, request of Chris)

Runs: corner K1 (L_eff 30 l_ex, K_eff 1.0 kJ/m³, K_σ 1.5 K_eff, φ 0.62, α 0.1, dx 3, f 30 MHz), seeds 2011 … 2013, at
A 20 pJ/m (d/l_ex 216, design runs K1_0_mc2, K1_1, K1_2) and A 10 pJ/m (d/l_ex 304, K1A10_0 … 2; 3 × 93 MB fetched, md5
identical with the instance, runs/MANIFEST_pair_K1A10.md5). PLAN 4.8 fixed no analysis for the pair; the rules of
`eval_pair.py` were written before the K1A10 values were read. Result: results/pair_K1A10.md. A lever = β(A 20) −
β(A 10), paired by seed, n 3, 90 % CI with t; all 7 amplitudes settled (≥ 2 calm cycles) in all 6 runs.

| quantity | measured (pair) | main fit at the same factors (90 % bootstrap interval) |
|---|---|---|
| A lever on β, 9 … 150 mT | +0.64 (+0.53 … +0.76; per seed +0.56 / +0.70 / +0.67) | +0.92 (+0.82 … +1.04) |
| A lever on β, 20 … 150 mT | +0.65 (+0.42 … +0.87) | +1.00 (+0.86 … +1.14) |
| P(A 10) / P(A 20) at 9 / 20 mT | × 3.11 (2.32 … 4.17) / × 1.61 (1.03 … 2.51) | × 5.61 (4.33 … 7.36) / × 3.27 (2.64 … 4.19) |
| P(A 10) / P(A 20) at 50 / 100 / 150 mT | × 1.26 (0.68 … 2.33) / × 0.58 (0.49 … 0.68) / × 0.51 (0.50 … 0.52) | × 1.34 / × 0.61 / × 0.50 |

- The A lever at K1 is real and does not depend on the 9 mT point (+0.64 with, +0.65 without 9 mT). β is 2.45 ± 0.12
  (SD, n 3) at A 20 and 1.81 ± 0.05 at A 10.
- Correction 2026-10-13 (review): the intervals of the A lever are now per seed (t with 2 degrees of freedom). The
  first version used a model with one offset per seed on all 21 values; it treats the amplitudes as independent and
  gave +0.55 … +0.74. eval_pair.py keeps it as a second column.
- The main fit over-predicts the A lever at this corner by ≈ 0.3 (the intervals do not overlap). The difference
  comes from small B: at 9 and 20 mT the model predicts about twice the measured loss ratio; at 50 … 150 mT the
  model and the pair agree within the intervals. The coded ln d/l_ex of the pair is −0.89 / +1.07 (A 10 is 0.07 coded
  units outside the design range).
- Why × 0.51 at 150 mT: (216 / 304)² = 0.505. At 150 mT the loss density at fixed d scales with A, i.e. with the
  exchange energy density A / d² of the magnetisation textures. At 9 mT the smaller A gives × 3 more loss.
- Loss channels (paired difference A 10 − A 20, n 3, 90 % CI; n 2 at 35 and 100 mT: no statement): at 9 mT the loss
  is more concentrated at A 10 (share of the top 10 % of the cells +0.17, +0.09 … +0.25); at 150 mT more of the loss
  is at f (R1 +0.12, +0.08 … +0.15) and the jump share is smaller (−0.04, −0.05 … −0.03). Other B: no clear
  difference.
- Limits: one design point (low K_eff, high stress), 3 seeds, coarse mesh dx 3 (the mesh question is open, B0);
  the cell-by-cell axes of the pair differ; the protocols differ (K1_1, K1_2: ≤ 7 cycles per stage without
  --min_calm; K1_0_mc2: ≤ 15 cycles, --min_calm 2; K1A10: ≤ 12 cycles, --min_calm 2, 512 samples per cycle).
- Field images (deck slide 10, the purpose of the pair, PLAN 4.8): turn angle θ = arccos(𝐦(t₊) · 𝐦(t₋)) between
  the two field peaks in the mid cut, seed 2013 (the seed with the median A lever: +0.56 / +0.70 / +0.67; rule fixed
  before the images were drawn). Mean over the mesh cells of the cut, mean of 3 seeds, 9 / 50 / 150 mT: 2 / 12 / 32°
  at A 20 and 3 / 12 / 32° at A 10: 𝐦 turns as far locally at both A. The cell mean ⟨𝐦⟩ (Herzer cell, whole volume,
  1.8 × the switched volume per half loop) turns 1 / 6 / 21° at A 20 and 1 / 6 / 16° at A 10: at 150 mT 25 % less at
  A 10 (seed ranges do not overlap).
- Correction 2026-10-13 (context-free review of the deck): the sentence "the same local motion gives × 3 more loss at
  9 mT and half the loss at 150 mT" was wrong. θ compares only the two field peaks, not the path of 𝐦. The motion
  measure q = Σ |Δ𝐦|² over the 16 snapshots (mean of the cut, eval_pair.py section 6), ratio A 10 / A 20: × 3.3 at
  9 mT (per seed 4.5 / 1.5 / 4.1; loss × 3.1), × 1.6 at 50 mT (loss × 1.3), × 0.93 at 150 mT (loss × 0.51). Thus at
  9 mT 𝐦 moves more at A 10, and the loss follows. At 150 mT the motion on the time scale T/16 is the same, and the
  loss halves: the loss at f falls to × 0.67, the loss outside f to × 0.42 (R1 · P, n 3). The difference at 150 mT
  sits in fast motion that the 16 snapshots do not resolve.
- Real cell-loss maps: the snapshot runs K1_0 … 2 (A 20) have no valid cell-loss slices (float16 overflow, old format,
  7.68: p = inf); only |𝐦₁|² (loss at f) is valid. Real loss maps at A 20 with the same protocol as K1A10 come from the
  C′ runs K1_k_f1 (A 20 at f, harmonics, 512 samples; no snapshots). The paired loss maps for slide 8 wait for them.

### 7.83 What the K1A10 runs teach: scaling of the loss with A, and the "other frequencies" (2026-10-13, items 3 and 1 of Chris)

**Item 3: how the loss scales with A** (`eval_scaling.py` → results/scaling_A.md; reading fixed in the script before
the result). Exponent s(B) = ∂ ln P / ∂ ln(d/l_ex) at fixed K_eff, stress and packing; at a fixed particle size
P ∝ A^(−s/2). Reading: s = −2 → P ∝ A (exchange energy A/d² of textures on the particle scale); s = −1 → P ∝ √A
(domain walls); s = 0 → no A dependence; s > 0 → the loss rises when A falls, more than any energy scale gives.
Design model (116 α 0.1 runs, Gamma GLM, 90 % bootstrap interval) and the pair (n 3, 90 % CI):

| B (mT) | L_eff 30 l_ex (low K_eff) | L_eff 19 | L_eff 12 (high K_eff) | pair K1 / K1A10 (L_eff 30) |
|---|---|---|---|---|
| 9 | +5.18 (+4.51 … +5.91) | +3.20 (+2.69 … +3.69) | +1.22 (+0.22 … +2.05) | +3.32 (+2.46 … +4.18) |
| 50 | +1.13 (+0.57 … +1.90) | +0.89 (+0.31 … +1.36) | +0.66 (−0.68 … +1.49) | +0.67 (−1.15 … +2.48) |
| 100 | −1.31 (−1.72 … −0.81) | −0.49 (−0.82 … −0.08) | +0.34 (−0.40 … +1.12) | −1.62 (−2.11 … −1.13) |
| 150 | −1.95 (−2.34 … −1.56) | −0.76 (−1.03 … −0.49) | +0.42 (−0.02 … +0.86) | −1.95 (−2.02 … −1.89) |

- At 150 mT and low K_eff: s = −2 in the main fit and in the pair (P ∝ A). At high K_eff: s = 0 (no A
  dependence). L_eff 19 (coded ln Q 0, the middle K_eff level of the deck) lies between: −0.76, the interval contains
  −1 (walls) and excludes 0.
- At 9 mT s = +1 … +5: a smaller A raises the loss strongly, most at low K_eff. This is not an energy scale; it is a
  change of the state or of its mobility (not resolved here).
- The A lever on β follows from the two ends: [s(150) − s(9)] · ln(300/212) / ln(150/9) = +0.88 / +0.49 / +0.10 at
  L_eff 30 / 19 / 12 (main fit; the direct A levers of the audit 7.66 at L_eff 30 / 18 / 12: +0.96 / +0.51 / +0.15)
  and +0.64 for the
  pair. Thus the A lever is the sum of two parts: more loss at small B and (at low K_eff) less loss at large B when A
  falls.
- The model and the pair agree at 100 … 150 mT; at 9 … 20 mT the model gives the larger s at L_eff 30 (7.82).
- The reading "exchange energy" for s = −2 is a hypothesis. A check with the non-uniformity g of the loop-mean 𝐦 in the
  snapshot cut (rad per particle diameter, 3 seeds: 4.2 … 4.3 at A 20, 4.4 … 4.6 at A 10, the same at 9 … 150 mT, +3 …
  +7 % per seed): the exchange energy A·g² of these textures falls to × 0.56 at A 10, the loss at 150 mT to × 0.51
  (0.50 … 0.52). Close, but outside the interval; g comes from one cut. The same g at 9 mT gives no explanation for
  × 3.1 there: the cause at small B is open.
- Limits: coarse mesh; quadratic model; the pair is one corner.

**Item 1: what the "other frequencies" of 7.70 are** (A 10 pJ/m only, K1A10_0 … 2, n 3; `eval_pair.py` section 5;
A 20 follows from the C′ runs K1_k_f1). Shares of the cell loss, mean ± SD over the seeds:

| B (mT) | f | 2 f | 3 f | 4 … 7 f | rest (> 7 f or not repeating) |
|---|---|---|---|---|---|
| 9 | 0.89 ± 0.04 | 0.03 ± 0.02 | 0.02 ± 0.01 | 0.04 ± 0.02 | 0.01 ± 0.01 |
| 50 | 0.54 ± 0.11 | 0.07 ± 0.01 | 0.05 ± 0.01 | 0.12 ± 0.04 | 0.22 ± 0.06 |
| 150 | 0.48 ± 0.04 | 0.11 ± 0.00 | 0.05 ± 0.00 | 0.10 ± 0.00 | 0.26 ± 0.04 |

- At ≥ 35 mT about half of the loss outside f is in the harmonics 2 … 7 f and half is in the rest. The reading rule of
  7.72 (rest large or small) does not give a clear class: both parts are large.
- 2 f is the largest harmonic and grows with B (0.03 → 0.11). The 4 … 7 f shares are flat in k (≈ 0.02 each): a
  step-like motion in each cycle and a slow drift within a cycle both give a flat spectrum; not separable here.
- The rest and the harmonics sit in the hot spots: at ≥ 50 mT the 10 % of the cells with the highest loss carry
  ≈ 90 % of the rest, ≈ 84 % of the harmonic loss and ≈ 60 % of the loss at f.
- The rest (0.22 … 0.27 at ≥ 50 mT) is larger than the jump share RF / 2 of the box (0.06 … 0.11): most of it is
  local fast or non-repeating motion that does not show as a jump of the box dissipation.
- Bessel check: 6 … 38 % of the cells have Σ_k p_k > 1.001 p (O2 test: only at 9 mT, 8 %; at 50 mT 0.01 %).
  Energetically small: the negative rest in the cut plane is ≤ 1.5 % of the loss (most stages < 0.5 %). The cause is
  not determined (candidates: drift within a cycle, motion between the samples). The 4 … 7 f shares and the rest carry
  this uncertainty (estimate: up to ≈ 0.1 can move between them).

**Item 2: where the channels sit** (A 10 only, `eval_channel_maps.py` → results/channel_maps_K1A10.md and .png;
descriptive). Concentration = share of the channel in a zone / volume share of the zone, mean ± SD over 3 seeds, 3D.

| zone (volume share) | B (mT) | total | f | 2 … 7 f | rest |
|---|---|---|---|---|---|
| texture top 10 % (0.10) | 50 | 7.1 ± 0.4 | 5.5 ± 0.1 | 8.2 ± 0.1 | 9.4 ± 0.3 |
| texture top 10 % (0.10) | 150 | 5.9 ± 0.1 | 4.7 ± 0.1 | 6.6 ± 0.3 | 7.5 ± 0.4 |
| centre r < 0.5 R (0.12) | 9 / 150 | 2.7 / 1.4 | 2.8 / 1.6 | 2.2 / 1.2 | – / 1.0 |
| rim r ≥ 0.8 R (0.49) | 9 / 150 | 0.54 / 1.02 | 0.49 / 0.79 | 0.88 / 1.15 | – / 1.35 |

- Every channel sits mostly where the loop-mean 𝐦 is most non-uniform; the further from f, the stronger: at ≥ 20 mT
  f ≈ × 5 … 6, harmonics ≈ × 7.5 … 8, rest ≈ × 9 … 11 of an even spread.
- The rim share of the total rises with B (7.69) for two reasons: every channel moves to the rim (f 0.49 → 0.79), and
  the harmonics and the rest, which sit more at the rim than f, grow.
- In the cut (K1A10_2, one plane) the harmonics and the rest form a spot in the middle of the central particle and
  bands at the edges of the corner particles; the loss at f is spread wider. The type of the structure at the spot
  (vortex core or wall) is not determined.
- At 9 mT the rest is ≈ 1 % and its zone shares are not meaningful (negative rest in cells with a Bessel excess).

### 7.84 Design model against the two direct tests of the A effect; the map cut (2026-10-13, context-free reviews of the deck)

**Model against direct tests** (`eval_model_points.py` → results/model_points.md; coarse mesh; A effect = β(A 20) −
β(A 10) at d = 1 µm, fit over 9 / 50 / 100 mT, 90 % intervals):

| point | measured | main fit at the same factors |
|---|---|---|
| corner K1 (K_eff 1.0 kJ/m³, K_σ 1.5 K_eff, φ 0.62), pair, n 3 paired | +0.67 (+0.61 … +0.73) | +0.91 (+0.78 … +1.03) |
| base point (K_eff 6.2 kJ/m³, no stress, φ 0.65), M212 / M300, n 4 + 4 (Welch) | +0.65 (+0.40 … +0.91) | +0.37 (+0.14 … +0.58) |

- The two direct tests give the same A effect at both K_eff ends (difference 0.02, about −0.2 … +0.2). The model gives
  a difference of 0.54. Thus the K_eff dependence of the A effect in the model (A × K_eff interaction, 7.66) is not
  seen in the 2 direct tests. At K1 the intervals do not overlap, at the base point they overlap.
- Limits: 2 points; the base-point runs use one job per amplitude from the relaxed virgin state (not the ascending
  design history), and different seeds at 212 and 300; the pair uses the design protocol. No generalisation.
- The deck now says this on the slides for A, the direct test, the mesh levers, validity and the key messages.

**The map cut**: the slices (acc_slices.npz, slices.npz) are the plane z = 0 (the first cell layer, a box face). The
plane holds only 2 of the 4 particles: particle 3 whole in the middle and particle 0 as 4 corner quarters (periodic
box); checked with the ids of K1A10_2 and S041. The deck said "cut through the 4 particles": corrected on the slides
7, 8, 11 and 15, and the claim "each particle holds a vortex" (slide 7) is removed (the cut shows 2 particles).

### 7.85 B0 stopped: energy balance fails on dx 0.75; precision test (2026-10-14, decision Chris: option B)

- B0 started 2026-10-13 21:31 UTC after C′. The 4 relaxations on dx 0.75 (d/l_ex 96) converged after 11 500 / 14 500 /
  10 500 / 5 500 iterations (≈ 0.98 s per iteration, 1.5 … 3.9 h); one cycle took 1.2 … 1.3 h (as measured in 7.80).
- The cycles were not calm: w_dis / w_loop = 2.0 / 1.9 / 2.4 (s941), 2.3 / 3.0 (s942), 1.44 / 1.05 / 1.12 (s943),
  3.9 / 5.1 / 3.6 / 3.7 (s944). In s944 the loop is periodic (closure 1·10⁻⁵, w_loop changes by 0.3 % from cycle to
  cycle), so the energy balance must hold; a ratio of 3.7 is a numerical error, not physics. On dx 3 and dx 1.5 the
  calm ratio is ≈ 1. Hypotheses: float32 rounding noise (exchange field ∝ 1/dx²) or the error of the time integrator
  (atol 1e-5) on the fine mesh.
- With this, 3 of 4 dx 0.75 runs would reach 12 cycles without a calm tail, and the pre-registered reading would be
  "not determined". Chris chose option B: B0 stopped 04:5x UTC (chain, queue and runs ended by PID, chain removed from
  onstart; the outputs and the 4 relaxed dx 0.75 states are kept).
- Precision test (jobs/prec48.txt, conv_jobs.py --prec): d/l_ex 48 (8 × fewer cells, the same dx / l_ex), base
  microstructure, seed 951, 9 mT, exactly 5 cycles: dx 0.75 single, dx 0.75 double, dx 0.75 single with atol 1e-6,
  dx 1.5 single. Reading (fixed now): the variant whose w_dis / w_loop stays within ± 0.2 of 1 in cycles 2 … 4 while
  dx 0.75 single does not names the cause (float32 or the time integrator); if dx 0.75 single is calm too, the test is
  not informative at this size; if no dx 0.75 variant is calm, the cause is elsewhere. Started 05:0x UTC.

### 7.86 Frequency test C′: result (2026-10-14; analysis of PLAN 4.8, eval_freq.py → results/freq_Cprime.md)

15 runs (corner K1, A 20 and A 10 pJ/m, seeds 2011 … 2013, f/2 / f / 2 f = 15 / 30 / 60 MHz, the same sample interval);
846 MB fetched, md5 identical with the instance (runs/MANIFEST_freq_Cprime.md5). All 7 amplitudes settled (≥ 2 calm
cycles) in all 18 run-frequency combinations (the f runs at A 10 are the pair runs K1A10_k). Checks of the method:
energy balance in the last cycle 0.90 … 1.07 at all stages; K1_0_f1 reproduces the design run K1_0 at f (the same w);
at 9 mT the state markers (offset of the loop, pairs above 60°) are the same at the 3 frequencies, so the frequency
dependence there is not a different state.

| A | B (mT) | w(f/2) / w(f) | w(2 f) / w(f) |
|---|---|---|---|
| 20 | 9 / 50 / 150 | × 0.34 / 0.34 / 0.30 | × 7.6 / 2.0 / 1.16 |
| 10 | 9 / 50 / 150 | × 0.17 / 0.29 / 0.59 | × 0.96 / 1.9 / 2.2 |

(w = loss per cycle, paired by seed, n 3; full table with 90 % CI in results/freq_Cprime.md.)

- The loss per cycle grows faster than ∝ f from f/2 to f at all B except A 10 at 100 … 150 mT (correction 7.88:
  there n = 2.07 and 1.75 (1.57 … 1.93) in P ∝ f^n). Thus w = w_h + c · f does not hold in general. In this model the
  loss at 15 … 60 MHz is a rate-dependent loss.
- β (9 … 150 mT): A 20 2.45 / 2.45 / 1.78, A 10 2.28 / 1.81 / 2.12 at f/2 / f / 2 f. A effect on β β(A 20) − β(A 10),
  paired: +0.16 (+0.00 … +0.32) at 15 MHz, +0.64 (+0.53 … +0.75) at 30 MHz, −0.34 (−0.39 … −0.29) at 60 MHz. The A
  effect on β is a property of the dynamics at 30 MHz, not a static property; it changes its sign between 30 and
  60 MHz. Answer to the question of PLAN 4.8: the A effect sits in the rate-dependent loss (a rate-independent part is
  not measurable here).
- Harmonics (7.72 reading rule): at 9 … 35 mT the shares move with f (e.g. A 20 at 20 mT: R1 0.92 / 0.81 / 0.68,
  2 … 7 × f 0.04 / 0.14 / 0.27): a rate effect. At ≥ 50 mT R1 depends little on f (A 10: 0.67 / 0.54 / 0.58 at 50 mT,
  0.48 / 0.48 / 0.46 at 150 mT): there the split is closer to a static nonlinearity.
- Bessel check (7.83): cells with Σ_k p_k > 1.001 p: 1 / 6 / 30 % (A 20) and 3 / 20 / 31 % (A 10) at f/2 / f / 2 f.
  Neither pre-registered case: not the same at all f (motion between the samples) and not largest at f/2 (drift
  over a long cycle); it grows with f. The rest at 2 f and 9 mT is negative (−0.09 … −0.10). Cause open.
- Limits: one corner (low K_eff, high stress), 3 seeds, coarse mesh dx 3 (the 9 mT point is suspect on this mesh,
  7.78), α 0.1. Above 9 mT the states at the 3 frequencies drift apart (state markers), so part of the difference at
  higher B can be state selection.

### 7.87 Final local evaluation: frequency groups at both A, loss location at both A, frequency exponent (2026-10-14)

Request of Chris (go 2026-10-14): evaluate what the files give, then the conclusion "what sets β". Scripts: eval_pair.py
(section 5 with A 20 from the C′ runs K1_k_f1), eval_channel_maps.py (both A), eval_freq.py (sections 6 … 8).

- Frequency groups at f, A 10 − A 20 paired by seed (n 3, 90 % CI): at 100 / 150 mT the share at f is higher at A 10
  (+0.10 / +0.11, intervals exclude 0) and the rest is lower (−0.11 / −0.15). At 9 … 50 mT no clear difference. Thus at
  large B the higher A loses more in fast or irregular motion. Rest share at 150 mT: 0.40 (A 20), 0.26 (A 10).
- Location (concentration c of the total loss, 3 seeds): in the 10 % with the highest non-uniformity g, c = 4.7 / 6.3
  / 4.9 at A 20 and 6.4 / 7.1 / 5.9 at A 10 (9 / 50 / 150 mT): the loss is somewhat more concentrated at A 10, most at
  small B. Rim and centre: the same trend with B at both A (centre 3.1 → 1.5 and 2.7 → 1.4, rim 0.45 → 1.04 and
  0.54 → 1.02).
- Exponent of P in f (P ∝ f^n; n = 1: loss per cycle independent of f, n = 2: loss per cycle ∝ f), paired, n 3:
  from 15 to 30 MHz n = 2.5 … 3.1 (A 20) and 1.75 … 3.6 (A 10) over B; from 30 to 60 MHz at 9 mT n = 3.9 (A 20) but
  0.94 (A 10); at 150 mT 1.2 (A 20) and 2.1 (A 10). The loss per cycle at 9 mT of A 10 rises × 6 from 15 to 30 MHz and
  stays from 30 to 60 MHz; that of A 20 rises × 2.9 and × 7.6. Loss ratio A 10 / A 20 at 9 mT: × 1.6 / 3.1 / 0.4 at
  15 / 30 / 60 MHz.
- Hypothesis (not tested; wording corrected 2026-10-14, see 7.94): the loss at small B comes from a slow local motion
  in the non-uniform zones, with a characteristic frequency far below the gyromagnetic resonance of the particles
  (GHz). This frequency lies between 15 and 60 MHz at A 10 and above 60 MHz at A 20. Candidates (not separated): a
  motion or relaxation of the texture that the exchange sets (rate ∝ A), or jumps above a threshold (depinning). n > 2
  cannot come from a linear relaxation (n ≤ 2 there); it needs the threshold process. Then A and the particle size set
  β at 30 MHz through the position of 30 MHz relative to this frequency. Test: a frequency sweep at 9 mT (e.g. 7.5 …
  240 MHz) for both A.
- Location against frequency: c in the 10 % with the highest g grows a little with f at small B (A 20 at 9 mT: 4.1 /
  4.7 / 5.7), and the rim share grows with f at all B.
- Energy weight of the Bessel excess: the value of this section (from the cut plane) is not valid. See 7.90.

### 7.88 Precision test result; third context-free review of the deck (2026-10-14)

**Precision test (7.85)**, d/l_ex 48, base microstructure, 9 mT, dx 0.75, cycles 1 … 4 (cycle 0 = transient):

| run | w_dis / w_loop | w_loop (cycles 1 / 2 / 3) | wall per cycle |
|---|---|---|---|
| single, atol 1e-5 | 1.96 … 1.99 | 8.320 / 8.439 / 8.443 ·10⁻⁷ | 730 … 755 s |
| double, atol 1e-5 | 2.49 … 2.51 | 8.319 / 8.438 ·10⁻⁷ | 1240 s |
| single, atol 1e-6 | 1.01 … 1.02 | 8.320 / 8.439 / 8.443 ·10⁻⁷ | 842 … 878 s |
| dx 1.5 single, atol 1e-5 (reference) | 1.25 … 1.26 | – | 219 s |

- Reading (fixed in 7.85): the variant with atol 1e-6 is calm while single at atol 1e-5 is not: the cause is the
  error of the time integrator, not float32 (double precision does not help, it is even worse). The step count is
  almost the same (≈ 59 000 per cycle); the cost of atol 1e-6 is + 14 … 16 % wall time.
- w_loop (the loss value of the study) is the same to 4 digits in all three variants: only the dissipation integral
  w_dis, and thus the calm check, is wrong at atol 1e-5 on fine meshes. The loss values of B0 and of the mesh checks
  are not affected; the steady rule is. On dx 1.5 the excess is already + 26 % at this size (the dx 1.5 mesh checks at
  d/l_ex 212 / 300 passed the rule).
- Consequence: a B0 rerun needs --atol 1e-6 (atol is an INIT_KEY: the dx 0.75 relaxations must run again).

**Third review of the deck** (35 slides; one agent): fixed: "established" holds on the coarse mesh only; the claim
"w grows faster than ∝ f at all B" (see the correction above); the hypothesis "at large B the loss follows the
exchange energy" is rejected, because P(A 10) / P(A 20) at 150 mT is × 0.51 (0.49 … 0.54) at 30 MHz but × 1.03
(0.75 … 1.40) at 15 MHz and × 0.96 (0.81 … 1.13) at 60 MHz (eval_freq.py section 5b; 100 mT: × 0.93 / 0.57 / 1.11):
a static energy scale cannot give this. Further: open points completed (no fine-mesh check at K1; the 9 mT point;
K_eff and stress; α; state drift in the frequency test), "vortex-like" replaced by "non-uniform zones", any linear
relaxation gives n ≤ 2, the prediction f₀(A 20) > f₀(A 10) named as the test of the eigen-motion hypothesis, n,
intervals and fit ranges added, the K_eff sign change marked as a model result, terms and status lines unified.

### 7.89 Frequency sweep at 9 mT, then B0 with atol 1e-6 (2026-10-14, decision Chris: option B; rules fixed now)

- Sweep (jobs/sweep9.txt, pair_jobs.py --sweep, 32 jobs): corner K1, A 20 (K1_k) and A 10 (K1A10_k), seeds 2011 …
  2013, only the 9 mT stage, started from the relaxed state of the same seed (--init_from, the start state of C′),
  harmonics k ≤ 7, samples max(256, 512 · 30 MHz / f), steady rule; f = 21 / 42 / 85 / 120 / 240 MHz (15 / 30 / 60 MHz
  from the first stage of C′). Control: K1_0 and K1A10_0 at 30 MHz must reproduce the C′ value (relative < 1e-3).
- Question: is the loss at 9 mT the low-frequency side of a damped eigen-motion whose frequency is lower at A 10 than
  at A 20? Reading (eval_sweep.py, synthetic test in test_study.py): "supports" if w(f) of A 10 has an interior maximum
  and the maximum of A 20 lies higher (more than one grid step, or at/above 240 MHz while still rising); "against" if
  both rise without an interior maximum with n ≤ 2 above 60 MHz, or both maxima lie within one grid step; else "not
  determined".
- Then B0 again (jobs/conv_b0.txt, names B0b_*): unchanged design and reading (PLAN 4.10), --atol 1e-6 for all arms
  (7.88). Chain instance/chain_sweep_b0.sh: starts after the precision test queue.
- Cost (estimate): sweep ≈ 7 … 12 GPU h, ≈ 2 … 3.5 h wall, ≈ 1 … 1.5 $; B0b ≈ 14 … 28 h, ≈ 6 … 12 $ (atol 1e-6 adds
  ≈ 15 % time).
- Rule added 2026-10-14 (7.90): the frequency groups (acc.json) of a sweep stage are read only if its net
  over-completeness e = max(Σ_{k ≤ 7} R_k − 1, 0) ≤ 0.02. At 60 MHz (256 samples per cycle) 7 of 42 C′ stages failed
  this check, 5 of them at 9 mT. The loss values (w_loop) are not affected.

### 7.90 Error of the frequency split: the cut-plane check was wrong (2026-10-14, question of Chris)

Question of Chris: can noise or jumps of 𝐦 make the per-cell Fourier split wrong, or do they only go into the rest?
Go of Chris for the patch of eval_freq.py (3D check), deck slides 9 and 11 and this record; second go ("Ok") for the
patch of acc_reduce.relative_slice, slide 9 sections 4 and 7, and the check of the cause of ε at 9 mT. Context-free
review of the correction done; its points are in this section.

- Mechanisms (properties of the method): (i) motion above 7f goes into the rest without error, because p̄ comes from
  the damping term at each sample, not from the Fourier series; (ii) motion that does not repeat from loop to loop
  (drift, changing jumps) gives the per-loop Fourier series a step at the end of the loop; the step adds a part of about
  equal size at each k (|m_k| ∝ 1/k, p_k ∝ k² |m_k|²). Thus the border between 2f … 7f and the rest is not sharp, and
  Σ_{k ≤ 7} p_k can exceed p̄ in a cell (negative rest); (iii) an excess of p̄ over the true loss (integrator error,
  or energy that the state gives away when it drifts) goes into the rest.
- Check of (iii) with the energy balance of the kept cycles, ε = w_dis / w_loop − 1, against the rest, 126 stages of
  18 C′ runs (the stages of one run are not independent; no error is given): ε = +0.009 (mean, SD 0.020, −0.016 …
  +0.117), rest 0.19 (mean), 0.26 at 35 … 150 mT.
  - At 20 … 150 mT (108 stages): ε ≤ +0.085, correlation of ε and rest −0.25. The rest there is not an excess of p̄.
  - At 9 mT (18 stages): correlation +0.55. At 15 MHz ε and rest are equal in 5 of 6 runs (ε +0.047 … +0.117, rest
    +0.048 … +0.105; 6th run ε +0.031, rest +0.104). At 30 MHz ε and rest are both ≤ 0.03.
  - Cause check (go of Chris; eval_freq.py section 9, local, no GPU). Rest without ε = R_rest (1 + ε) − ε (the rest
    relative to w_loop when the excess sits in the rest). At 15 MHz and 9 mT: −0.002 … +0.004 in 5 of 6 runs (6th
    +0.076). In these 5 runs the loop closes (|ΔM| / M_peak ≤ 0.0001) and ε is the same in each kept cycle (spread
    ≤ 0.021 over 3 … 5 cycles): the state is periodic. A periodic state with ε > 0 breaks the energy balance (exact
    dynamics give w_dis = w_loop over a closed cycle). Thus ε there is an error of the computed trajectory or of the
    sampled integral, not a loss, and not a drift of the state (a drift would not close the loop and would change ε
    from cycle to cycle). Established for these 5 stages; which numerical part (integrator tolerance atol 1e-5, as in
    7.88, or the sampling of p at 1024 points per cycle) is not separated. Excess power (W_dis − W_loop) f: 0.02 …
    0.04 W/cm³ at 15 MHz; at 30 MHz 0.003 … 0.045 W/cm³ in the 3 runs with a steady ε (A 20). A source with about
    the same power per time (time step 5.9 ps at all f) gives a large ε where the loss per cycle is small: 9 mT,
    15 MHz.
  - Effect on the reported values: the loss values (w_loop) are not affected. The rest at 9 mT and 15 MHz is ≈ 0 after
    the correction. At 20 … 150 mT the correction changes the mean rest by ≤ 0.026 (A 10, 30 MHz, 35 mT: 0.200 →
    0.174); at 50 … 150 mT and 30 MHz by ≤ 0.009. The deck values (rest at 30 MHz, 50 … 150 mT; the 9 mT rest is not
    used) do not change; no slide change.
  - At 60 MHz and 9 mT the state is less periodic: ε changes from cycle to cycle by up to 0.3 (A 20) and the loop
    closes less well (|ΔM| / M_peak 0.0004 … 0.0067, against ≤ 0.0006 at 15 and 30 MHz). This fits mechanism (ii) as
    the cause of the over-completeness there. Not proven: the run with the largest e (K1A10_0_f2, e 0.173) has a steady
    ε (+0.027, +0.028, 2 cycles) and |ΔM| / M_peak 0.003.
- Finding: the check "energy weight of the Bessel excess in the cut plane ≤ 0.6 %" (eval_freq.py section 8, 7.87,
  deck slides 9 and 11) is wrong. Cause: acc_reduce.relative_slice stores a plane relative to its mean and stores zeros
  when the mean is ≤ 0 (the signed mean itself is stored). A plane with a net negative rest thus shows rest 0. This
  hits the 60 MHz stages with the largest excess and also 3 planes at 30 MHz. Local acc_slices.npz with a zeroed rest
  plane (of 194 files): K1A10_0_f2 amp00/01, K1A10_1 amp00/02, K1A10_1_f2 amp00/05, K1A10_2_f2 amp00, K1_0_f2
  amp00/01, K1_1_f1 amp00, K1_1_f2 amp00, K1_2_f2 amp00. The two runs of the deck maps (K1A10_2, K1_2_f1) have no
  zeroed plane. A second weakness of relative_slice: a very small positive mean gives very large relative values, which
  can overflow float16 (maximum 65504).
- New check (eval_freq.py section 8): (a) 3D from acc.json, net over-completeness e = max(Σ_{k ≤ 7} R_k − 1, 0) per
  stage. e is a lower bound of the per-cell excess x = Σ max(Σ_k p_k − p, 0) / Σ p (cells with a positive rest
  compensate). (b) x in the cut plane z = 0 (2 of the 4 particles) where the plane mean rest is > 0; for a zeroed plane
  the stored signed mean gives a lower bound −rest_mean / p_mean.

| A (pJ/m) | 15 MHz | 30 MHz | 60 MHz |
|---|---|---|---|
| 20, e (3D) | 0 in 21 of 21 stages | 0 in 21 of 21 | up to 0.127 (seed 0, 9 mT); 5 of 21 stages > 0.02 |
| 20, x (cut) | ≤ 0.009 (21 planes) | ≤ 0.004 (20 planes); 1 zeroed, ≥ 0.001 | ≤ 0.059 (17 planes); 4 zeroed, up to ≥ 0.25 |
| 10, e (3D) | 0 in 21 of 21 | 0 in 21 of 21 | up to 0.173 (seed 0, 9 mT); 2 of 21 stages > 0.02 |
| 10, x (cut) | ≤ 0.008 (21 planes) | ≤ 0.015 (19 planes); 2 zeroed, ≥ 0.001 | ≤ 0.031 (16 planes); 5 zeroed, up to ≥ 0.38 |

- Cells with Σ_k p_k > 1.001 p̄ (section 5, mean over seeds and B): 1 / 6 / 30 % (A 20) and 3 / 20 / 31 % (A 10) at
  15 / 30 / 60 MHz. At 60 MHz, 7 of 42 stages have e > 0.02 (5 of them at 9 mT, 2 in K1_0_f2 at 20 and 35 mT).
- Consequences: the loss values, β, the frequency ratios and the A effect are not affected (they come from w_loop).
  At 15 and 30 MHz the excess is small where it can be checked (x ≤ 1.5 % of the loss in 81 of 84 cut planes, e = 0 in
  3D); the frequency groups of the pair and the channel maps (30 MHz) keep their values, with this limit. At 60 MHz the
  groups are not reliable at 9 mT and in K1_0_f2 also at 20 and 35 mT. The rest at 9 mT and 15 MHz is mostly an excess
  of p̄ (see the check of (iii)). Deck: slide 9 sections 5 and 6 and slide 11 corrected; the numbers come from
  freq_Cprime.json.
- Patch of acc_reduce.relative_slice (2026-10-14): a plane is stored relative to its mean |value| (key
  <stage>_<x>_scale), so a signed plane with a mean ≤ 0 is kept; float32 when a relative value does not fit into
  float16. Readers (eval_freq, eval_channel_maps, deck) use slice_values / plane_known and read old and new files; the
  results of the old files are unchanged (freq_Cprime, channel_maps_K1A10 and the deck figures byte-identical). The
  sweep and B0b chain on the instance runs the old commit 371004f: its acc_slices.npz have the old format. The per-cell
  accumulators (acc_*.npz) are deleted after the postprocessing, so old planes cannot be recomputed.
- Open: the cause of the over-completeness at 60 MHz (motion that does not repeat, supported by the cycle-to-cycle
  change of ε, or short events between the 256 samples per cycle); which numerical part gives ε at 9 mT and 15 MHz (a
  run with atol 1e-6 would separate it; GPU).

### 7.91 Frequency sweep at 9 mT: result (2026-10-14; reading fixed in 7.89)

- Run 06:52 … 08:37 UTC on the instance, 32 jobs, no failure, post-processing exit code 0. Data local: runs/*_sw*
  without checkpoint.pt (23 MB, 258 files, md5 equal to the instance, runs/MANIFEST_sweep9.md5); results/sweep9.json
  and results/sweep9_report.md (eval_sweep.py). B0b (atol 1e-6) started at 08:37 UTC (chain 7.89).
- Control: w at 30 MHz equals the C′ value (relative difference 2.9e-6 at A 20, 1.1e-6 at A 10; limit 1e-3).

| f (MHz) | 15 | 21 | 30 | 42 | 60 | 85 | 120 | 240 |
|---|---|---|---|---|---|---|---|---|
| P A 20 (W/cm³, mean ± SD of 3 seeds) | 0.39 ± 0.05 | 0.85 ± 0.13 | 2.29 ± 0.55 | 14.5 ± 5.9 | 33.5 ± 4.0 | 23.5 ± 4.9 | 85.1 ± 14.8 | 219 ± 23 (n 2) |
| P A 10 (W/cm³) | 0.59 ± 0.12 | 2.57 ± 0.31 | 6.92 ± 0.36 | 5.37 ± 0.31 | 13.5 ± 0.7 | 39.6 ± 5.7 | 46.9 ± 3.3 | 222 ± 30 |
| P(A 10) / P(A 20), paired, 90 % CI | 1.53 (0.94 … 2.49) | 3.01 (2.15 … 4.22) | 3.07 (2.24 … 4.22) | 0.39 (0.17 … 0.93) | 0.40 (0.35 … 0.46) | 1.70 (1.50 … 1.93) | 0.56 (0.45 … 0.69) | n 2 |

- The values at 15, 30 and 60 MHz are those of C′ (first stage, 7.89); the sweep ran 21, 42, 85, 120 and 240 MHz and
  the two 30 MHz controls.
- Pre-registered reading (on the loss per cycle): **not determined**. The largest loss per cycle lies at the edge
  (240 MHz) for both A. The rule as coded and
  tested before the runs (eval_sweep.reading) counts a maximum at the edge as "no maximum found": "supports" needs an
  interior maximum of A 10; "against" needs either both curves rising with n ≤ 2 above 60 MHz (here n up to 3.7) or
  two interior maxima within one grid step. Open point for Chris: the wording of 7.89 ("both maxima lie within one grid
  step") can also be read as "against" for two edge maxima at the same f. The coded reading is used here.
- Observation, not pre-registered (3 seeds, no claim): the loss is not monotonic in f. A 10 has a local maximum of
  the loss per cycle at 30 MHz and a dip at 42 MHz, A 20 at 60 and 85 MHz. Both rise again above. In P the dips are
  13 … 38 % per seed (7.97). The local maxima lie lower
  at A 10 than at A 20, in the direction of hypothesis H1 (7.87), but the fixed rule does not support H1. The seeds of
  A 20 at 42 MHz scatter strongly (P 8.3 … 20.2 W/cm³). Causes of the dips (several eigen-motions, or a different steady
  state at each frequency): not checked.
- Next step: decision of Chris.

### 7.92 Sweep at 9 mT: is the steady state the same at all frequencies? (2026-10-14, option 1 of Chris)

Local check, no GPU (eval_sweep_state.py → results/sweep9_state.md / .json). 48 runs: corner K1, A 20 and A 10, seeds
0 … 2, 8 frequencies, all from the same relaxed start state of the seed. Markers per run (mean over the kept cycles):
the mean of M along the field in each of the 4 particles (Mp0 … Mp3), the box mean ⟨M⟩, n60 (pairs above 60°); and the
energy that the state gives away in the start-up events.

- Reading fixed before the markers were looked at (docstring): two runs are in the same state if every marker differs
  by at most 3 × the larger cycle-to-cycle spread + a floor. Result: **same state** at both dips. A 10 at 42 MHz: same
  state as a neighbour frequency in 3 of 3 seeds, A 20 at 85 MHz in 2 of 3. Thus the dips of w(f) do not come from
  another steady state (as far as these markers see). All 3 seeds show the dips (A 10: P 5.06 … 5.69 W/cm³ at 42 MHz
  against 6.52 … 7.23 at 30 MHz; A 20: P 17.9 … 27.3 W/cm³ at 85 MHz against 28.9 … 36.4 at 60 MHz).
- Added after the reading, against the 30 MHz state of the same seed: clear differences (ratio ≥ 2) in 3 of 42 runs,
  all at A 20 (seed 2 at 42 and 85 MHz, seed 0 at 120 MHz). A 20 gives away more energy in the start-up at ≥ 42 MHz
  (up to 5.9 µJ/cm³) than at 30 MHz and below (at most 0.05 µJ/cm³). A 10: no ratio above 1.6.
- Added after the reading (observation, not pre-registered): the loss energy per loop E(f) (µJ/cm³) of A 20 has the
  shape of the A 10 curve shifted by the factor 2 in f. SD of ln E20(s f) − ln E10(f) over the common frequencies
  (mean curves; per seed pair):
  s 1: 0.83 (0.84 ± 0.07); s √2: 0.63 (0.65 ± 0.10); **s 2: 0.17 (0.23 ± 0.01)**; s 2√2: 0.78 (0.80 ± 0.15).
  At s 2, E20(2 f) ≈ 2.4 × E10(f) (6 frequencies, 15 … 120 MHz of A 10), so P20(2 f) ≈ 4.8 × P10(f).
- Hypothesis (not tested; it extends H1 of 7.87; wording corrected 2026-10-14, see 7.94): at 9 mT the loss comes
  from a slow local motion whose characteristic frequencies scale with A, far below the gyromagnetic resonance (GHz).
  An exchange-dominated motion or relaxation of a texture of fixed size scales so (rate ∝ γ A / (M_s L²)); a
  magnetostatic gyration of a vortex core would hardly depend on A. Then A sets the loss at a fixed drive frequency
  through the position of the drive relative to these frequencies.
- Limits: the markers are coarse (a change inside a particle that keeps its mean, such as a core shift or a chirality,
  is not seen); 3 seeds of one corner (K1), 9 mT, coarse mesh, α 0.1; the shift is tested only in steps of √2; the
  observation was found after the reading, so it needs a test fixed in advance.
- Next step: decision of Chris.

### 7.93 Sweep at 9 mT: loss maps of A 10 at f against A 20 at scaled f (2026-10-14, option 1 of Chris, second step)

Local check, no GPU (eval_sweep_maps.py → results/sweep9_maps.md / .json). Question: does A 20 at 2 f lose in the same
zones as A 10 at f (hypothesis of 7.92)? Cut z = 0, 9 mT; the A 10 map (152 cells per edge) as block mean on the A 20
grid (108); Pearson correlation r of ln p over the common particle cells, for A 20 at √2^s · f, s 0 … 3,
f = 15 / 21 / 30 / 42 MHz; R_s = mean over these f.

| seed | R_0 | R_1 | R_2 | R_3 |
|---|---|---|---|---|
| 0 | 0.672 | 0.717 | **0.730** | 0.660 |
| 1 | **0.583** | 0.541 | 0.498 | 0.426 |
| 2 | 0.722 | 0.756 | 0.820 | **0.856** |

- Reading (fixed before the comparison): **not determined**. Largest R_s at s 2 / 0 / 3; "supports" needs s 2 in all 3
  seeds, "against" s 0 in at least 2. The map of the motion at f (|m1|², reported only) gives s 1 / 0 / 2.
- Why the test does not decide (observation): the maps follow the state. Seed 1 of A 20 is in another state at 60 and
  120 MHz (7.92); every pair with A 20 at 60 MHz has r 0.31 … 0.39. Seed 2 of A 20 changes its state at 42 MHz and
  above; its maps at 42 … 240 MHz resemble each other (r up to 0.91 against A 10). Within one seed, r changes by 0.1 …
  0.3 from one f to the next, as much as between the shifts.
- Consequence: the map test cannot separate the shift hypothesis from the state changes of A 20. The w(f) shift of 7.92
  stays an observation found after the reading. A decisive test needs a prediction fixed in advance on new runs.
- Next step: decision of Chris.

### 7.94 Correction: the characteristic frequency of H1 is not a resonance of the particles (2026-10-14, remark of Chris)

Remark of Chris: the resonance frequency is much larger. Check (Snoek's law, Snoek 1948): (μ − 1) f_r ≈ (2/3) (γ/2π)
J_s, γ/2π = 28 GHz/T, J_s = 1.5 T. With the powder μ_r of the simulation (6.5 … 6.9): f_r ≈ 5 GHz; for a sphere at the
demagnetisation limit (χ ≈ 3): f_r ≈ 9 GHz. The data agree: at 30 MHz μ_r is constant and the loss angle δ is about
0.4°; near a resonance δ would be large. Thus the drive lies far below the resonance of the particles, and the features
of w(f) at 30 … 85 MHz are not that resonance. The wording "eigen-motion below a weakly damped resonance" (7.87, 7.92,
deck slide 32) was misleading and is corrected (go of Chris for the patch). New wording: a slow local motion in the
non-uniform zones with a characteristic frequency far below the gyromagnetic resonance. It can carry the loss and
change δ only a little, because 10 % of the volume carry about 74 % of the loss. Candidates (hypothesis, not tested):
(i) gyration of a vortex core (in µm-sized elements typically 0.1 … 1 GHz, literature on discs, e.g. Guslienko 2002;
not checked for spheres; magnetostatic, so hardly dependent on A: it does not fit the shift by the factor 2 of 7.92);
(ii) a motion or relaxation of a texture of fixed size that the exchange sets (rate ∝ A γ / (M_s L²); fits the shift);
(iii) jumps above a threshold (depinning): the only candidate for the local n up to 3.9 … 5.2, which a linear
relaxation cannot give; the start-up events and the state changes of A 20 at ≥ 42 MHz (7.92) point to such processes.
The pre-registered question of 7.89 ("damped eigen-motion") stays as written; its reading (7.91) does not change.

Addition (2026-10-14, remark of Chris "so small differences of f cannot give different physics"; go for the patch):
the GHz resonance holds for the particle as a whole (against its demagnetising field). Inside the particles, motions
that keep the flux closure (rotation in walls and textures) feel only the anisotropy field μ0 H_K = 2 K / M_s. Their
time scale (Kittel, γ/2π = 28 GHz/T, α 0.1; estimate):

| K | μ0 H_K | f_K = γ μ0 H_K / 2π | relaxation time 1 / (α 2π f_K) |
|---|---|---|---|
| 1.0 kJ/m³ (K1, K_eff only) | 1.7 mT | 47 MHz | 34 ns |
| about 2.5 kJ/m³ (K1 with K_σ 1.5 K_eff, rough sum) | 4.2 mT | 117 MHz | 14 ns |
| 6.2 kJ/m³ (high K_eff) | 10.4 mT | 290 MHz | 5.5 ns |

The drive periods are 67 / 33 / 17 ns at 15 / 30 / 60 MHz. Thus in this simulation the drive lies close to the internal
time scale, and differences between 15 and 60 MHz are possible. Still open: n > 2 needs a threshold process (candidate:
pinning at the coarse mesh, the core moves less than one cell at 9 mT; earlier analysis eval_pinning.py), and f_K does
not depend on A at fixed K_eff (it does not explain the shift by the factor 2 of 7.92).
Transfer (estimate): real nanocrystalline powders have a much smaller K_eff (a few J/m³ from the grains, Herzer 1990;
up to some 100 J/m³ with stress). Their internal time scale is then of the order of 0.1 … 5 MHz, and a drive at 30 MHz
lies above it, in another regime than in the simulation. The frequency results of the simulation may not transfer to a real
powder.

### 7.95 Frequency data at 9 mT: what was computed, and how accurate (2026-10-14, question of Chris)

Question of Chris: the time scale of 7.94 does not depend on A; then something does not fit. Statistics or accuracy?
Audit (eval_sweep_audit.py → results/sweep9_audit.md, all 48 runs listed with their values).

- What was computed: corner K1, seeds 2011 … 2013, 9 mT only, from the relaxed start state of each seed; 15 / 30 /
  60 MHz = first stage of the C′ runs, 21 / 42 / 85 / 120 / 240 MHz = sweep; drive correction after cycle 0, 7 … 12
  cycles, steady rule ≥ 2 calm cycles; 1024 … 256 samples per cycle; coarse mesh, atol 1e-5, single precision.
- The run code works in units of l_ex at A 10 pJ/m. The "A 20" and "A 10" runs differ only in the particle size:
  d/l_ex 216 against 304 (0.72 against 1.02 µm at A 10 pJ/m, box 108 against 152 cells). The same in both: K_eff
  995 J/m³, K_p 1492 J/m³, L_eff 100 nm, cell 3 l_ex, α 0.1, J_s 1.5 T, packing, stress axes, drive, frequency. With
  the particle size the number of Herzer cells per particle changes (about 316 against 790), and the grain axes differ.
- Accuracy of each value: reached B_peak 8.89 … 9.10 mT (at most 3 % in the loss); cycle-to-cycle SD of the loss
  median 3.6 %, largest 18 % (A 10, 42 MHz, seed 0); choice of the cycles up to 20 % (7.97). The dips are in all 3
  seeds: P falls by 13 … 30 % (A 10, 30 → 42 MHz) and 23 … 38 % (A 20, 60 → 85 MHz). In the loss per cycle they looked
  like a factor 2, because the frequency rises at the same time. Seed scatter is large where A 20 changes its state
  (42 MHz: P 8.3 … 20.2 W/cm³).
- Consistency: all time scales that do not contain d are equal in both runs: the anisotropy time of 7.94 (the same
  K_eff), the exchange time on the scale of the grain cells or of the mesh, the magnetostatic time. They cannot shift
  the curve. The only time scale that changes is the exchange time on the particle scale, ∝ A / d²; its ratio is
  (304 / 216)² = 1.98, the shift factor 2 found in 7.92. Thus the observation is consistent with an exchange-governed
  motion on the particle scale, not with the anisotropy time. The anisotropy estimate of 7.94 only shows that the
  frequency range is not quasi-static for internal motions.
- Limits: the shift test of 7.92 rests on two features per curve (local maximum and dip); a pure power law fits every
  shift equally, so the test is weak. A mesh artefact cannot give an A-dependent shift through the cell size (the cell
  is 3 l_ex in both runs); through the particle size (101 against 72 cells per diameter) it is not excluded. All values
  are on the coarse mesh.

### 7.96 Where the factor 2 is, and where it can come from (2026-10-14, question of Chris)

- Where: on the frequency axis at 9 mT. The w(f) curve of the smaller particle ("A 20", d/l_ex 216) has its features at
  twice the frequency of the larger one ("A 10", 304): local maximum 60 against 30 MHz, dip 85 against 42 MHz (7.92). At
  these corresponding frequencies the smaller particle loses about 2.4 × more energy per loop (P about 4.8 × more,
  because f is doubled). The factor is not in β.
- At larger B it does not hold (eval_sweep_audit.py, C′ with 15 / 30 / 60 MHz, E the loss energy per loop): for
  E20(2 f) = c · E10(f) the two shifted ratios E20(30)/E10(15) and E20(60)/E10(30) must be equal. At 9 mT 1.92 and 2.43; at 20 … 150 mT 2.72 … 3.67
  against 1.58 … 2.27. Thus the shift is an observation at 9 mT only.
- Where it can come from, when only the particle size differs (estimates):
  (i) exchange on the particle scale: a pattern that scales with d has a restoring field ∝ A / d², ratio
  (304 / 216)² = 1.98. Size: μ0 H_ex = (2 A / M_s) (2π / λ)² gives 0.64 mT (18 MHz) for λ = 1.02 µm and 1.3 mT
  (36 MHz) for λ = 0.72 µm. These fields are smaller than the anisotropy fields, which do not depend on d (K_eff
  1.7 mT, stress K_p 2.5 mT). The total restoring field then changes much less than × 2. A clean factor 2 needs pattern
  features of at most about 0.25 µm (μ0 H_ex ≥ 10 mT). Not consistent so far.
  (ii) number of grain cells per particle (316 against about 790): the random part of the anisotropy averaged over a
  particle is about K_eff / √N = 35 … 56 J/m³, small against K_p = 1492 J/m³. Negligible.
  (iii) mesh: 72 against 101 cells per diameter; pinning at the mesh is not excluded, but gives no simple factor 2.
- Conclusion: the factor 2 at 9 mT has no consistent explanation. It rests on two features of one flux density and can
  be chance. No further interpretation until fine-mesh data exist. The patch proposal for slide 32 (factor 2 as
  (d/l_ex)²) is withdrawn.

### 7.97 Loss density at 9 mT with its scatter; why the seed scatter is small here (2026-10-14, question of Chris)

Table (eval_sweep_audit.py, steady rule of eval_sweep): P = f K_d w in W/cm³ (core volume), mean ± SD between the
seeds (3 independent microstructures), n seeds with a steady value.

| f (MHz) | large particle ("A 10", d/l_ex 304) | small particle ("A 20", d/l_ex 216) |
|---|---|---|
| 15 | 0.59 ± 0.12 (n 3) | 0.39 ± 0.05 (n 3) |
| 21 | 2.57 ± 0.31 (n 3) | 0.85 ± 0.13 (n 3) |
| 30 | 6.92 ± 0.36 (n 3) | 2.29 ± 0.55 (n 3) |
| 42 | 5.37 ± 0.31 (n 3) | 14.5 ± 5.9 (n 3) |
| 60 | 13.5 ± 0.7 (n 3) | 33.5 ± 4.0 (n 3) |
| 85 | 39.6 ± 5.7 (n 3) | 23.5 ± 4.9 (n 3) |
| 120 | 46.9 ± 3.3 (n 3) | 85.1 ± 14.8 (n 3) |
| 240 | 222 ± 30 (n 3) | 219 ± 23 (n 2; seed 2 has no steady value after 12 cycles) |

- Dips, paired by seed (the same microstructure at both f), factor on P, 90 % CI with t(2): large particle 30 → 42 MHz
  × 0.78 (0.64 … 0.94; seeds 0.70 / 0.87 / 0.77); small particle 60 → 85 MHz × 0.69 (0.57 … 0.84; 0.77 / 0.62 /
  0.69). The cycle choice changes single values by up to 20 % (median 0 %); the dips stay with both rules in all seeds.
- Remark of Chris: in his simulations the loss at small B often varies by about 80 % between seeds. Why the scatter is
  smaller here: (i) the seeds differ only in the random axes and offsets of the grain cells (about 790 per particle,
  so the random part averages out strongly) and in the random start state before the relaxation. The particle
  positions (fcc), the packing, the particle size and shape and the stress axes (37° and 79°) are the same for all
  seeds. Thus the seed scatter of the simulation is only a part of the variation of a real powder. (ii) n 3 cannot show a small
  scatter: the true SD can be up to 4.4 × the sample SD (90 % bound, χ² with 2 degrees of freedom). The centre point
  of the design (n 8, the same rules) has an SD of ln P of 0.16 … 0.32 over B, max/min up to 2.25; the corner groups
  (n 3) 0.03 … 0.40. The SD of 5 … 7 % of the large particle at 30 … 60 MHz is likely too small.
- Consequence: comparisons across frequency are paired (same seed), so most of the microstructure scatter cancels.
  Error bars of absolute values and of comparisons between different microstructures must use the larger scatter of
  the centre group, and they still miss the variation that the simulation model does not contain (packing, size, shape, random
  stress axes).
- Verdict (question of Chris "noise or effect?"): (a) the dips and the factor-2 shift: not shown. Their size
  (13 … 38 % per seed) is about the uncertainty of the method (cycle choice up to 20 %); "all 3 seeds fall" has a
  probability of 1/8 = 0.125 without any effect (sign test); the t-interval rests on an SD from 3 values. Treat them as
  noise of the method until the fine mesh shows them. (b) Effects in the simulation, far larger than every uncertainty:
  P rises × 9.7 … 15 (large particle) and × 5.3 … 6.5 (small particle) per seed from 15 to 30 MHz at 9 mT; the loss
  ratio of the two particle sizes turns from × 3.07 (2.24 … 4.22) at 30 MHz to × 0.40 (0.35 … 0.46) at 60 MHz. All on
  the coarse mesh, where 9 mT has the largest mesh error: whether they are physical stays open (B0b, fine mesh).

### 7.98 Why the A effect on β changes its sign with f (2026-10-14, question of Chris)

- Data: C′ = 18 runs (2 particle sizes × 3 frequencies × 3 seeds). Each run passes the 7 amplitudes (9 → 150 mT) one
  after the other, and each stage starts from the end state of the stage before. Thus 126 stage values, but per
  particle size and frequency only 3 independent values (the seeds); β is one value per run. The sweep (7.91) ran
  only 9 mT and gives no β.
- β per frequency and particle size (eval_freq.py sections 3, 10, 11; mean ± SD of 3 seeds; an SD from 3 values
  can be up to 4.4 × too small, 7.97). Large particle = A 10 runs (d/l_ex 304), small = A 20 runs (216).

| | 15 MHz | 30 MHz | 60 MHz |
|---|---|---|---|
| β 9 … 150 mT, large particle | 2.28 ± 0.12 | 1.81 ± 0.05 | 2.12 ± 0.03 |
| β 9 … 150 mT, small particle | 2.45 ± 0.03 | 2.45 ± 0.11 | 1.78 ± 0.01 |
| A effect β(small) − β(large), paired, 90 % CI | +0.16 (+0.00 … +0.32) | +0.64 (+0.53 … +0.75) | −0.34 (−0.39 … −0.29) |
| the same, fit 20 … 150 mT (without 9 mT) | +0.16 (+0.10 … +0.23) | +0.65 (+0.43 … +0.86) | −0.17 (−0.43 … +0.09) |
| the same, fit 9 … 100 mT (without 150 mT) | +0.18 (−0.06 … +0.43) | +0.65 (+0.55 … +0.75) | −0.45 (−0.50 … −0.40) |

Local β in windows of 3 amplitudes (figure results/freq_beta_windows.png):

| window (mT) | large 15 | large 30 | large 60 | small 15 | small 30 | small 60 |
|---|---|---|---|---|---|---|
| 9 … 35 | 2.10 ± 0.22 | 1.92 ± 0.18 | 2.28 ± 0.01 | 2.27 ± 0.04 | 2.48 ± 0.08 | 1.61 ± 0.13 |
| 20 … 50 | 2.36 ± 0.25 | 1.97 ± 0.20 | 2.28 ± 0.10 | 2.66 ± 0.07 | 2.23 ± 0.46 | 1.66 ± 0.20 |
| 35 … 70 | 2.49 ± 0.30 | 1.60 ± 0.23 | 1.96 ± 0.25 | 2.44 ± 0.27 | 2.47 ± 0.40 | 1.83 ± 0.04 |
| 50 … 100 | 2.38 ± 0.03 | 1.68 ± 0.20 | 1.96 ± 0.16 | 2.59 ± 0.12 | 2.81 ± 0.41 | 2.02 ± 0.10 |
| 70 … 150 | 2.32 ± 0.13 | 1.78 ± 0.05 | 1.89 ± 0.30 | 2.54 ± 0.24 | 2.37 ± 0.17 | 2.04 ± 0.06 |

- Reading: at 30 MHz the large particle has the lowest β in every window (1.60 … 1.97), the small one 2.23 … 2.81;
  the A effect on β is the same in every fit range (+0.64 … +0.66): it does not rest on the 9 mT point or on the
  150 mT point. At 60 MHz the small particle has a low β at small B (1.61 … 1.83 up to 70 mT) and the large one a high
  β there (2.28); the sign change rests mostly on the range up to 50 mT (without 9 mT: −0.17, the interval contains 0).
  At 15 MHz both lie at 2.1 … 2.7, the small particle a little higher.
- Correction (same day): a first version of this section split the A effect into a "9 mT end" (+0.40) and a "150 mT
  end" (+0.24) and said that about 60 % of it rests on the 9 mT point. That split referred both ends to an equal loss
  (ratio 1), an arbitrary reference; the test without the end points above shows that the statement was wrong.

### 7.99 What happens at 60 MHz: the motion of the small particle stops repeating at small B (2026-10-14)

Question of Chris: 60 MHz is the outlier; what happens there, and what is the characteristic frequency? Local check
(eval_sweep_regime.py → results/sweep9_regime.md): per run and stage, from the last 3 cycles, RN = share of the
dissipation that changes from cycle to cycle, RF = share above the smooth periodic curve, loop closure, span of the
energy balance ε over the 3 cycles; mean of 3 seeds.

| 9 mT | 15 | 21 | 30 | 42 | 60 | 85 | 120 | 240 MHz |
|---|---|---|---|---|---|---|---|---|
| small particle (A 20 runs): RN | 0.020 | 0.009 | 0.004 | 0.037 | **0.201** | 0.086 | 0.178 | 0.130 |
| small: ε span | 0.002 | 0.002 | 0.001 | 0.034 | **0.523** | 0.165 | 0.564 | 0.377 |
| large particle (A 10 runs): RN | 0.023 | 0.016 | 0.042 | 0.041 | 0.048 | **0.115** | 0.107 | 0.072 |
| large: ε span | 0.007 | 0.026 | 0.031 | 0.127 | 0.107 | **0.330** | 0.323 | 0.143 |

- At 60 MHz the motion of the small particle no longer repeats at small B: 20 % of the dissipation changes from cycle
  to cycle (0.4 % at 30 MHz), the loop does not close (closure 0.0047 against 0.0000), the energy balance of single
  cycles spreads by 0.52. Over B (C′): at 60 MHz RN 0.10 … 0.20, largest at 9 … 35 mT; at 30 MHz RN 0.004 … 0.07 up
  to 35 mT and 0.13 … 0.31 above 50 mT.
- Effect on β, small particle 30 → 60 MHz: P rises × 15 at 9 mT but only × 2.3 at 150 mT (paired, from the loss per
  cycle × 7.6 and × 1.16), so β falls from 2.45 to 1.78. Large particle: P × 1.9 at 9 mT (the loss per cycle stays,
  × 0.96) and × 4.3 at 150 mT, so β rises from 1.81 to 2.12. Both changes together give the sign change of the A
  effect on β.
- Characteristic frequency, measured as the onset of the irregular motion at 9 mT: small particle between 42 and
  60 MHz, large particle between 60 and 85 MHz. Estimate (7.94): the relaxation time of a motion that the anisotropy
  field restores is 14 … 34 ns at K1; the drive period falls below it at about 30 … 74 MHz. Then the magnetisation
  cannot settle in each half cycle. This time does not depend on A or d; it explains the range, not the difference
  between the sizes. The onset is lower for the smaller particle, opposite to the shift by the factor 2 of 7.92: the
  shift hypothesis is not supported by this measure.
- Physical or numerical: irregular motion of a driven nonlinear system can be real, or it can be seeded by the
  integrator tolerance (atol 1e-5). Not separated by these data.

### 7.100 Are the evaluated loops settled? Stop rule checked; long tolerance test at 60 MHz (2026-10-14)

Remark of Chris: when the loop does not close, the stage has not settled. True for the small particle at 60 MHz and
9 mT: P jumps by ±15 % from cycle to cycle to the end (26 … 40 W/cm³), the loop does not close (0.001 … 0.010), and
K1_2_f2 still has events in cycles 6 and 7 (ε +1.4, +1.2); the large particle at 60 MHz drifts (22 → 14 W/cm³). At
30 MHz (K1_0_f1) P is constant from cycle 3 and the loop is closed. Cause (method): the stop rule (≥ 2 calm cycles at
the end, |w_dis / w_loop − 1| ≤ 0.2, at most 12 cycles) checks the energy balance of each cycle, not whether the cycles
repeat. Thus the 60 MHz values at 9 mT (dips, factor 2, β 1.78 of the small particle, the sign change of the A effect)
rest on cycles that have not settled.

Option 1 of Chris (local, eval_periodic.py → results/periodic.md / .json): every stage of the 116 runs of the main
analysis (30 MHz) and of C′, last 3 cycles. Classes fixed before the data were looked at: periodic RN ≤ 0.05 and
closure ≤ 0.001; not periodic RN > 0.10 or closure > 0.002; marginal otherwise.

| B (mT) | 9 | 20 | 35 | 50 | 70 | 100 | 150 |
|---|---|---|---|---|---|---|---|
| main data: periodic / marginal / not periodic | 93 / 7 / 16 | 68 / 16 / 32 | 48 / 9 / 59 | 30 / 18 / 68 | 21 / 22 / 73 | 9 / 28 / 79 | 6 / 48 / 62 |
| main data: w not stationary (added after the classes) | 19 | 30 | 33 | 29 | 27 | 31 | 20 |

- The class "not periodic" comes mostly from RN (the time pattern of the dissipation changes from cycle to cycle, for
  example jumps at other times); the median closure is at most 0.0007. RN does not test whether the loss value has
  settled. Added after the classes were looked at: the loss per cycle is not stationary (|drift| over the 3 cycles or
  cycle SD of w above 5 %) in 19 … 33 of 116 stages per B (16 … 28 %). The drift has no sign: up and down about equal,
  mean drift +0.1 … +0.8 % at most B (SE 0.4 … 1.0 %; 35 mT +3.9 %, SE 5.8 %, from single outliers). Thus the
  unsettled stages add noise to single values, not a bias; the bootstrap intervals of the main fit contain this noise
  through the scatter between the runs. Single values (and n 3 comparisons) carry an extra 5 … 10 % in about a quarter
  of the stages.
- At large B the motion is intermittent (it differs from cycle to cycle) but the loss per cycle is stationary to about
  2 % (median cycle SD 1.6 … 2.3 %): a statistically steady state with jumps, not an unsettled one.

Option 2 of Chris (instead of the 3 short jobs of 7.99; after B0b on the same instance; jobs/atol60.txt from
pair_jobs.py --atol60): the small particle (K1_k, seeds 0 … 2) at 60 MHz, 9 mT only, exactly 40 cycles, with atol 1e-5
and 1e-6 (6 jobs, names K1_k_sw060L5 / L6; otherwise the sweep line). Cost (estimate) about 4 GPU h, 0.5 … 1 $.
Reading (fixed now, eval_atol60.py, synthetic test in test_study.py):
  per run, from cycle 2 on: c_settle = first cycle from which every later window of 3 cycles has RN ≤ 0.05, closure
  ≤ 0.001 and w within ±5 % of the mean of the last 10 cycles; "periodic" if c_settle exists, "stays irregular" if the
  last 10 cycles have median RN > 0.10 or median closure > 0.002, else "marginal".
  tolerance artefact      with 1e-6 periodic with c_settle ≤ 12 in ≥ 2 of 3 seeds, with 1e-5 not periodic by cycle 12
                          in ≥ 2 of 3 seeds;
  settles slowly          both tolerances periodic in ≥ 2 of 3 seeds, and c_settle > 12 in ≥ 2 of 3 seeds for one of them
                          (the limit of 12 cycles was too short);
  stays irregular         both tolerances "stays irregular" in ≥ 2 of 3 seeds;
  otherwise               not determined.
  Loss value: P = mean of the last 10 cycles per run; "the C′ value holds" if it lies within ±10 % of the 9 mT value of
  K1_k_f2 for both tolerances in ≥ 2 of 3 seeds, "the C′ value is off" if outside ±10 % in ≥ 2 of 3 seeds, else not
  determined.

### 7.101 Loop closure of the main run: recurrent jumps about a fixed state; closure independent of atol (2026-10-14)

Question of Chris: are nearly half of the 116 runs of the main analysis not settled, because the amplitudes rise
without a relaxation between them? Why did the rule of 2 calm cycles not prevent this? Answer: no. The class "not
periodic" of 7.100 means that the time pattern of the dissipation changes from cycle to cycle (RN). It does not mean an
open loop. The numbers below come from config.json, summary.json and samples_*.csv of the 116 runs (α 0.1, 30 MHz).
The check ran locally, without a GPU.

Numbers: eval_closure_state.py → results/closure_state.md / .json (all tables of this section; the B0 / B0b sample
files were fetched, md5 of the B0 files in runs/MANIFEST_b0_samples.md5; the B0b files are partial copies of a
running job and are replaced at the final fetch).

Terms of this section. Event: a cycle with |w_dis / w_loop − 1| > 0.2 or w_loop ≤ 0 (7.53). Calm cycle: the opposite
of an event. Steady tail: the calm cycles after the last event of a stage (7.54). Closure: |M(end) − M(start)| / M_peak
of one cycle. ΔM: the same difference with its sign. Δe = w_loop − w_dis: the change of the internal energy over one
cycle from H = 0 to H = 0, in units of K_d (energy theorem). Half loop width: half the opening of the loop at H = 0,
M_peak sin δ. For an elliptic loop sin δ = w_loop / (2π h b), with h = H_peak / M_s and b = B_peak / J_s (estimate).

**Protocol of the main run.** virgin_asc: random start, relaxation at H = 0 until converged, then 7 stages
9 → 20 → 35 → 50 → 70 → 100 → 150 mT. The stages follow directly. There is no relaxation and no pause between them.
The amplitude changes at H = 0. Each stage starts from the end state of the stage before.
- 103 runs: 7 cycles per stage, no stop rule during the run.
- 13 reruns (_mc2): at least 7 cycles, then more until the steady tail has 2 cycles, at most 15 (7 … 12 used).
- Analysis: cycle 0 (drive correction) and cycle 1 are discarded. The other cycles are measured. The steady-tail rule
  is applied afterwards.
- Thus the rule of 2 calm cycles acted during the run in 13 of 116 runs only. In the other runs it is a filter of the
  analysis.

**Calm cycles.** All 812 stages (116 runs × 7 amplitudes) end with a steady tail of at least 2 cycles (minimum 2, at
150 mT 4). All measured cycles are calm in 70 / 98 / 91 / 102 / 112 / 109 / 114 of 116 stages (9 … 150 mT). The rule
checks only the energy balance of each cycle (tolerance 20 %). It does not check the closure and the repetition of the
cycles. Thus a steady tail does not mean "settled" in the sense of 7.100: w is not stationary to 5 % in 16 … 28 % of
the stages (7.100). For C′ the case 60 MHz at 9 mT stays open (P drift of the large particle, 7.100; atol60 test).

**Closure of the main run** (mean of the cycles of the steady tail per stage; number of stages of 116 above a limit):

| B (mT) | 9 | 20 | 35 | 50 | 70 | 100 | 150 | sum (of 812) |
|---|---|---|---|---|---|---|---|---|
| > 2 % (closure_tol of the run) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| > 0.5 % | 1 | 1 | 0 | 0 | 1 | 1 | 0 | 4 |
| > 0.2 % | 9 | 5 | 13 | 7 | 12 | 12 | 12 | 70 |
| > 0.1 % | 23 | 20 | 35 | 35 | 35 | 40 | 47 | 235 |
| median | 0.02 % | 0.03 % | 0.05 % | 0.06 % | 0.05 % | 0.08 % | 0.08 % | |

Against the half loop width, the closure has a median of 1 … 10 %. It is above 25 % of the half width in 9 … 18 % of
the cycles, and above 50 % in about 3 %.

**Monotonic change or recurrent jumps.** Criteria (fixed before the data were looked at), per B over the steady tails:
"recurrent jumps" if |mean Δe / w_loop| ≤ max(2 SE, 0.02) and the share of consecutive cycles with the same sign of ΔM
is ≤ 0.5 + 2 SE; "still relaxing" if mean Δe / w_loop < −max(2 SE, 0.02) or the share is > 0.5 + 2 SE; else not
determined.

| B (mT) | Δe / w_loop per cycle, mean ± SE (n 116 stages) | consecutive cycles with the same sign of ΔM | stages with an event in the measured cycles | events: energy release / uptake |
|---|---|---|---|---|
| 9 | −0.74 ± 0.25 % | 39 % (398 pairs) | 46 | 62 / 2 |
| 20 | −0.68 ± 0.17 % | 35 % (431) | 18 | 29 / 0 |
| 35 | −0.18 ± 0.17 % | 36 % (409) | 25 | 48 / 0 |
| 50 | −0.27 ± 0.09 % | 25 % (439) | 14 | 20 / 0 |
| 70 | −0.31 ± 0.11 % | 30 % (453) | 4 | 6 / 0 |
| 100 | −0.23 ± 0.13 % | 29 % (450) | 7 | 9 / 0 |
| 150 | +0.04 ± 0.06 % | 25 % (462) | 2 | 2 / 0 |

- Reading: "recurrent jumps" at all B.
- Expected same-sign share. The criterion above used 0.5, the value of a wandering state. The correct reference for
  a fixed state is 1/3 (added after the data were looked at):

| picture | M at H = 0 from cycle to cycle | expected share of consecutive cycles with the same sign of ΔM |
|---|---|---|
| fixed state with scatter | independent values about a fixed mean | 1/3 |
| wandering state | independent steps, M moves away freely | 1/2 |
| monotonic change | M moves in one direction | nearly 1 |

  Why 1/3: ΔM of cycle c is M_{c+1} − M_c, with M_c the value at the start of cycle c. Two consecutive ΔM have the
  same sign only if M_c, M_{c+1} and M_{c+2} are in rising or in falling order. For independent values from any
  continuous distribution, the probability of this is 2 / 3! = 1/3.
  Limit of the test (simulation of the second review): a monotonic change of M gives a share near 1 only when the
  change per cycle is large against the cycle-to-cycle scatter of M. A change of 1 × the scatter per cycle gives 0.55,
  and 0.3 × the scatter gives 0.36. Thus the share cannot find a monotonic change smaller than about the scatter of M.
  For the loss, the stationarity of w (7.100) covers this case.
- SE of the share (added after the review): bootstrap over whole stages, 2000 draws. Inside one stage, neighbouring
  pairs share one ΔM and are negatively correlated. This alone makes the binomial SE too large. Differences between
  the stages make it too small. The bootstrap SE is 0 … 20 % above the binomial SE at the measured share.

| B (mT) | 9 | 20 | 35 | 50 | 70 | 100 | 150 | all B |
|---|---|---|---|---|---|---|---|---|
| same-sign share | 38.9 % | 34.8 % | 36.4 % | 24.6 % | 29.6 % | 29.1 % | 24.7 % | 30.9 % |
| SE, stage bootstrap | 2.7 % | 2.5 % | 2.8 % | 2.1 % | 2.3 % | 2.4 % | 2.3 % | 0.9 % |
| distance to 1/2, in SE | −4.1 | −6.0 | −4.9 | −12.0 | −8.9 | −8.7 | −11.2 | −20.3 |
| distance to 1/3, in SE | +2.1 | +0.6 | +1.1 | −4.1 | −1.6 | −1.8 | −3.8 | −2.6 |

- Reading (after the data were looked at): at all B the share is 4 … 12 SE below 1/2. Thus M does not wander, and it
  does not change monotonically by more than its scatter.
- At 20, 35, 70 and 100 mT the share agrees with 1/3 (within 2 SE).
- At 50 and 150 mT the share is about 4 SE below 1/3. Thus a ΔM is followed by a ΔM of the opposite sign more often
  than for independent values of M. An example is a change back and forth from cycle to cycle.
- At 9 mT the share is 2.1 SE above 1/3: a weak persistent part.
- C′ (21 stages per group): 24 … 34 %, within 1.9 SE of 1/3 and 3.3 … 5.7 SE below 1/2.
- The share does not give the size of the scatter. The scatter of M at H = 0 shows in the closure (median
  0.02 … 0.08 % of M_peak, table above). The cycle SD of w is a different quantity (median 0.4 … 2.3 %,
  results/periodic.md).
- The unsettled part is the start of a stage. After the amplitude step, the state still changes irreversibly in the
  measured cycles in 2 … 46 of 116 stages per B. 176 of 178 events release energy. The steady-tail rule removes these
  cycles. At least 2 calm cycles remain after the last event.
- More relaxation time between the amplitudes would not have helped. The 9 mT stage starts from the converged
  relaxation at H = 0, and it has the most events. The larger field switches regions for the first time. A relaxation
  at H = 0 stays in the nearest local minimum and cannot do this. More cycles per stage would give more calm cycles.
- Δe / w_loop is below 0 by 2.8 … 4.0 SE at 9, 20, 50 and 70 mT (−0.3 … −0.7 %), but within the 2 % limit. Candidate
  cause: the integrator error of w_dis at atol 1e-5 (7.88), not an energy release. Estimate for dx 3: K1_k at 9 mT has
  a closure of 0 from cycle 4 … 5 on and still w_dis / w_loop 1.015 … 1.03, that is w_dis 1.5 … 3 % too large. w_loop,
  and thus P, does not depend on atol (7.88, and B0 / B0b below).

**C′ per frequency** (corner K1, 3 seeds × 7 amplitudes = 21 stages per group; the same measures; steady rule during
the run, at most 12 cycles):

| f | particle | stages with a steady tail ≥ 2 | largest number of cycles | events release / uptake | median closure | same sign of ΔM | median cycle SD of w |
|---|---|---|---|---|---|---|---|
| 15 MHz | large | 21 / 21 | 7 | 4 / 0 | 0.006 % | 28 % | 1.2 % |
| 15 MHz | small | 21 / 21 | 7 | 3 / 0 | 0.004 % | 24 % | 1.6 % |
| 30 MHz | large | 21 / 21 | 9 | 9 / 0 | 0.048 % | 28 % | 2.8 % |
| 30 MHz | small | 21 / 21 | 8 | 4 / 0 | 0.053 % | 33 % | 2.3 % |
| 60 MHz | large | 21 / 21 | 7 | 4 / 0 | 0.068 % | 34 % | 2.4 % |
| 60 MHz | small | 21 / 21 | 12 | 18 / 0 | 0.101 % | 26 % | 4.4 % |

- All 126 stages have a steady tail of at least 2 cycles. In no group does M change monotonically by more than its
  scatter (limit of the test above). These are medians over all B. They do not replace 7.100 for 9 mT at 60 MHz.
- 15 MHz is the cleanest group: its closure is 8 … 13 × smaller than at 30 MHz.
- The small particle at 60 MHz has 2 … 6 × more events. It needs up to 12 cycles (the limit) and has twice the cycle
  scatter of w. This is an observation from 3 seeds, not a proven effect.
- The closure grows with f (0.005 % → 0.05 % → 0.07 … 0.10 %). Candidate (estimate, Kittel with H_K, 7.94): a jump
  rings down with τ = 1 / (α 2π f_K) ≈ 5 … 34 ns. The period is 67 / 33 / 17 ns. An integrator error would rather grow
  at 15 MHz, because a cycle has more steps there. This is a hint, not a proof.

**Is the closure an integrator error? B0 against B0b (go of Chris, reading fixed before the comparison).**
- Reading: "physical" if the median closure with atol 1e-6 is at least 0.5 × the median with atol 1e-5, per seed;
  else "numerical".
- Data: B0 (atol 1e-5, stopped in 7.85) has 2 … 4 cycles per seed in samples_00_amp00.csv. Its summary.json does not
  exist, because run_loops writes it at the end of a stage. B0b (atol 1e-6, running) has the same cycles. Same seeds
  941 … 944, d/l_ex 96, dx 0.75, 9 mT.
- Result: in all 12 common cycles of the 4 seeds the closure with 1e-6 is 0.9998 … 1.0014 × the closure with 1e-5.
  This holds also for large values (0.3 … 3 %) and through the jumps (s943: w_loop × 3.8 in 3 cycles). W_loop agrees to
  5·10⁻⁵. The relaxations need the same number of iterations (11 500 / 14 500 / 10 500 / 5 500). w_dis differs (s941
  cycle 2: 1.09·10⁻⁶ against 5.12·10⁻⁷), as in 7.88.
- **Reading: physical** in 4 of 4 seeds (in the measured cycles from cycle 2 on: 3 of 3 seeds; B0 s942 has no cycle 2).
- Limits: dx 0.75, not the mesh of the main run (dx 3); 9 mT only; 1 … 4 cycles per seed, mostly the first cycles of
  the stage.
- The precision test P48 (d/l_ex 48, dx 0.75, 9 mT) agrees: the closure per cycle is the same to 5 digits for single
  1e-5, double 1e-5 and single 1e-6.

**Error in the atol60 jobs, found before the runs.**
- atol is an INIT_KEY (7.88). A job with --init_from K1_k/init.pt (relaxed with atol 1e-5) and --atol 1e-6 stops at
  once ("relaxed with other arguments"). The 3 L6 jobs of jobs/atol60.txt (7.100) had this error.
- Cause: the job builder did not compare the INIT_KEYS with the source run. The dry run tested only the evaluation.
- Fix (patch, go of Chris): a run with atol 1e-6 relaxes again with atol 1e-6 from the same seed, as a production run
  with atol 1e-6 would. The relaxation of K1_0 took 298 s. The L5 jobs keep --init_from.
- New test (test_init_from_jobs_have_the_init_keys_of_their_source): every job line with --init_from <run>/init.pt has
  the INIT_KEYS of the line of <run>. The test fails on the old jobs/atol60.txt (3 lines, atol). It passes on all 194
  such lines now.
- At 60 MHz, L5 and L6 now start from two relaxations. eval_atol60.py reports a start control (E / K_d and w_loop of
  cycle 0, L6 relative to L5). In B0 / B0b the relaxed state was the same for both atol. The readings of 7.100 do not
  change.

**Tolerance twins at 30 MHz** (decision Chris: option 1 with patch; still in the chain, see the remark at the end).
- Jobs jobs/atol30.txt (pair_jobs.py --atol30): K1_k_L6 (k 0 … 2) is the production line of K1_k (7 amplitudes, 7
  cycles each, virgin relaxation) with --atol 1e-6. Only the snapshots are removed (disk). Reference: K1_k (atol 1e-5,
  the same line; not K1_k_mc2).
- Reading (eval_atol30.py, synthetic test in test_study.py, fixed now):
  - per run, c = median closure over the measured cycles (cycles 2 … 6 of all 7 stages, no filter);
  - per seed, r_k = c(K1_k_L6) / c(K1_k); the reference values are c(K1_k) = 0.055 / 0.057 / 0.054 %;
  - "physical" if r_k ≥ 0.5 in at least 2 of 3 seeds, "numerical" if r_k < 0.5 in at least 2 of 3 seeds;
  - "not determined" if a run is missing or incomplete, or if c(K1_k) = 0, so that no 2 seeds agree.
- Controls (reported, not part of the reading): the relaxed start state of K1_k_L6 against K1_k (E / K_d, M / M_s);
  w_loop of cycle 0 at 9 mT ("same start" if within 1e-3 relative, else the comparison is statistical only).
- Also reported, not part of the reading: per B the closure, w_dis / w_loop and P (W/cm³) with both tolerances.
- Output: runs/atol30_report.md and results/atol30.json.
- Chain instance/chain_atol30.sh: pin /root/atol30_commit, read after the wait. The chain starts when the atol60 queue
  has ended (runs/ATOL60_QUEUE_DONE) and the atol60 chain has released its lock.
- The atol60 chain now also reads its pin after the wait, and it deletes an old ATOL60_QUEUE_DONE before its queue. It
  was restarted once with the new pin.
- Failure mode: if the atol60 queue never ends (for example, B0b gives no report), the atol30 chain waits. Its log
  shows the wait. The state is checked at each status report.
- Cost (estimate from the measured cycle time of K1_k, 1.2 h per run, + 15 % for atol 1e-6, + about 5 min relaxation):
  about 4.4 GPU h, about 1.5 h wall time, about 0.6 $.
- Remark: B0 / B0b already give "physical" on dx 0.75 at 9 mT. The twins test the mesh of the main run (dx 3) and
  35 … 150 mT, where the closure of K1 is. Chris can remove them from the chain without cost before they start.

### 7.102 Results of atol60, atol30 and B0b; restart of the instance; state files on the disk (2026-10-15)

Terms used here: P = f K_d w_loop (W/cm³). R1, R2 … R7: share of the loss at k × the drive frequency (per-cell
Fourier sums, 7.72); rest = 1 − Σ R_k (7.90). RN: share of the dissipation that changes from cycle to cycle; RF: share
that is not smooth within the cycle (7.51, 7.99). Steady tail, event, closure, ΔM: 7.101. L5 / L6: atol 1e-5 / 1e-6.
F, I15, I3, T15, T3: the B0b arms (below). Each pre-registered reading comes first; every later result is marked
"exploratory".

**Restart of the instance.** The instance restarted on 2026-10-14 at 18:47 UTC (cause not known). The chains resumed
through the onstart lines.
- B0b_T15_s941 had written its init.pt at 18:41. After the restart the file had 0 bytes: run_loops renamed the file
  without os.fsync, and its data were still in the page cache. Both attempts after the restart stopped with EOFError
  in torch.load.
- F s941 and F s942 were inside their 9 mT stage and repeated it from init.pt (a run with one stage has no checkpoint
  inside the stage). T15 s942 restarted during its relaxation. F s943, F s944 and the other runs had finished before
  the restart.
- Fix (patch, go of Chris), commit 515b260: run_loops.replace_durable flushes the file (os.fsync) before os.replace and
  the directory after it. save_atomic (init.pt, checkpoint.pt), write_atomic (summary.json, config.json, DONE), the
  accumulator files and run_queue.save_queue_state use it. The status file of run_queue (monitoring only) does not.
  New test test_state_files_are_on_the_disk_before_the_rename.
- Rerun of B0b_T15_s941 at commit 515b260 (only file writing differs from 371004f of the other B0b runs). The failed
  folder is kept as runs/B0b_T15_s941.failed_reboot, the first report as runs/conv_b0_report_n19.md. The B0b chain was
  started again with /root/sweep_commit 515b260; it skipped the finished runs. B0b_T15_s941 finished at 05:03 UTC (5
  cycles, periodic, P 4.13 W/cm³), the queue ended with exit code 0, and eval_conv.py wrote the report with n = 4. The
  local run of eval_conv.py gives the same report.
- Check of all files on the instance after the restart (summary.json readable, all stages complete, rows of the
  samples files = cycles × (samples + 1), every .pt loaded with torch.load): all 29 runs of B0b, atol60 and atol30 are
  complete. Local copies without .pt: md5 equal (runs/MANIFEST_b0b.md5, MANIFEST_atol60.md5, MANIFEST_atol30.md5).
- Reproducibility floor (exploratory, results/b0b_check.md): the first attempt of F s941 / s942 relaxed in the same
  process, the second started from init.pt (another start path, 7.59). w_loop of cycles 0 … 3 is the same to the 4
  printed digits, the closure to 2 digits. w_dis differs by up to 3·10⁻³ (s941) and 1.3·10⁻² (s942). The B0 / B0b
  comparison of 7.101 with the final files: closure ratio 0.9998 … 1.0014 and W_loop within 5·10⁻⁵, as stated there
  (results/closure_state.md).

**atol60** (reading fixed in 7.100; L6 relaxes itself, 7.101; runs/atol60_report.md, md5 runs/MANIFEST_atol60.md5).
Small particle (K1_k, A 20 runs), 60 MHz, 9 mT, 40 cycles:

| seed | class 1e-6 / 1e-5 | P of the last 10 cycles (W/cm³, mean ± SD), 1e-6 / 1e-5 | C′ value (W/cm³) |
|---|---|---|---|
| 0 | stays irregular / stays irregular | 35.7 ± 1.8 / 35.7 ± 1.8 | 35.3 |
| 1 | stays irregular / stays irregular | 28.5 ± 1.1 / 28.5 ± 1.1 | 28.9 |
| 2 | marginal / marginal | 38.1 ± 4.4 / 38.1 ± 4.4 | 36.4 |

- Start control: the two relaxations give the same state (E / K_d equal, w_loop of cycle 0 within 1·10⁻⁵).
- **Reading (motion): stays irregular. Reading (loss value): the C′ value holds.**
- The irregular motion is not an error of the time integrator (the mesh is not tested here). Over 40 cycles it stays
  irregular in seeds 0 and 1 and marginal in seed 2. P scatters by 4 … 12 % (SD over the last 10 cycles). Its mean
  agrees with the C′ value within 5 %.
- Exploratory (after the readings, results/atol60.md): w_loop of L6 and L5 agree cycle by cycle to ≤ 1·10⁻⁵ in seeds 1
  and 2 over all 40 cycles. In seed 0 the difference grows from 6·10⁻⁶ to 1.2·10⁻⁴ (cycles 30 … 39). The irregular
  motion is reproducible. The running mean of P from cycle 2 on stays within ±5 % of the mean of the last 10 cycles
  after cycle 2 / 10 / 6 (seeds 0 / 1 / 2), with both tolerances.
- Exploratory: loss channels at 60 MHz, 9 mT. L6 and L5 agree to 0.003: R1 0.915 … 0.923, R2 0.033 … 0.053, R4 … R7
  0.038 … 0.064. C′ K1_k_f2 kept 4 cycles (seed 0: 2, 4, 5, 6; seed 1: 3 … 6; seed 2: 4, 5, 8, 9): R1 0.854 … 0.882 and
  R4 … R7 0.11 … 0.14. Candidate cause: the C′ window holds early cycles with more harmonic content. This is not
  measured (R_k of L5 over the C′ cycles was not computed). The zone shares agree within 0.03.

**atol30** (reading fixed in 7.101; runs/atol30_report.md, md5 runs/MANIFEST_atol30.md5). Main-run line of K1_k
(30 MHz, 7 amplitudes, 7 cycles each) with atol 1e-6 against the main runs K1_k (atol 1e-5):
- c (median closure of the measured cycles of all 7 stages) with 1e-5 / 1e-6: 0.0547 / 0.0545 %, 0.0570 / 0.0570 %,
  0.0539 / 0.0538 %; r = 1.00 in all 3 seeds. The same start in all 3 seeds (E / K_d equal, w_loop of cycle 0 within
  3·10⁻⁵).
- **Reading: physical** (3 of 3 seeds). At K1, the loop closure of the main run is not an error of the time
  integrator.
- P with 1e-6 / 1e-5 (mean of 3 seeds): 1.000 at 9 … 100 mT; 0.992 at 150 mT, where the trajectories of the two
  tolerances separate (median closure 0.0020 with 1e-5, 0.0019 with 1e-6). At K1 the loss values do not depend on the tolerance to 1 %.
- w_dis / w_loop at 9 mT: 1.021 with 1e-5, 1.001 with 1e-6. This confirms the w_dis error of atol 1e-5 at 9 mT (7.88,
  7.101) at K1 (n 3). The error acts only on the calm check.
- Change of w_dis / w_loop from 1e-5 to 1e-6: −0.28 % at 20 mT, at most 0.07 % at 35 … 100 mT, +0.4 % at 150 mT. Thus
  at 20 mT the tolerance gives part (≈ 0.3 %) of the mean Δe / w_loop of −0.68 % of 7.101. At 50 and 70 mT it does not
  give the −0.3 % of 7.101; that cause is open (K1 only, n 3).
- Exploratory (after the reading, results/atol30.md): the event cycles are the same with both tolerances in 21 of 21
  stages, so at K1 the steady-tail rule does not depend on atol. β over 9 … 150 mT per seed (main-fit definition):
  2.311 / 2.306, 2.553 / 2.553, 2.482 / 2.480 (1e-5 / 1e-6). The windows 9 … 35 and 35 … 70 mT are equal to 4 digits;
  70 … 150 mT changes by −0.023 / 0.000 / −0.012. The lag of the fundamental is equal to 0.01°. The loss shares per
  zone change by ≤ 0.003 at 20 … 150 mT and by 0.014 at 9 mT.
- Exploratory: R1 at 9 mT is 0.971 with 1e-5 and 0.990 with 1e-6 (mean of 3 seeds); 0.804 / 0.806 at 20 mT; equal
  to 0.001 at 35 … 150 mT. The denominator of R1 is the accumulated dissipation, which atol 1e-5 makes too large at
  9 mT. Thus at 9 mT and 30 MHz the non-fundamental share 1 − R1 of the 1e-5 runs (0.029) is about 2/3 numerical (0.010
  with 1e-6), at K1 (n 3).

**B0b** (mesh convergence at 9 mT, PLAN 4.10, reading fixed before the runs; runs/conv_b0_report.md, md5
runs/MANIFEST_b0b.md5). Base microstructure (L_eff 12 l_ex, r_p 0, φ 0.65, α 0.1), d/l_ex 96, 9 mT, seeds 941 … 944,
atol 1e-6. F: dx 0.75. I15 / I3: dx 1.5 / dx 3, independent relaxation. T15 / T3: dx 1.5 / dx 3, start from the relaxed
F state, then relaxed again.

| arm | runs | P (W/cm³), geometric mean (a) | SD of ln P over the seeds (a) | ln P per seed 941 / 942 / 943 / 944 (b) | median cycle SD of P (b) | median RN (b) |
|---|---|---|---|---|---|---|
| F, dx 0.75 | 4 | 10.5 | 0.96 | 2.51 / 1.97 / 3.59 / 1.33 | 5.2 % | 0.061 |
| I15, dx 1.5 | 4 | 3.25 | 0.65 | 1.44 / 0.97 / 1.92 / 0.39 | 0.3 % | 0.004 |
| I3, dx 3 | 4 | 1.47 | 0.09 | 0.52 / 0.31 / 0.39 / 0.33 | 0.0 % | 0.001 |
| T15, dx 1.5 | 4 | 4.02 | 0.54 | 1.42 / 2.09 / 0.79 / 1.26 | 0.1 % | 0.001 |
| T3, dx 3 | 4 | 1.51 | 0.07 | 0.48 / 0.47 / 0.33 / 0.39 | 0.0 % | 0.000 |

(a) eval_conv.py: ln of the mean w_loop of the steady tail per run. (b) eval_b0b_check.py: mean of ln P over the tail
cycles per run; it differs from (a) in the third digit (T15: 4.01 / 0.540).

Paired differences of ln P (same seed, n 4, mean ± SE, 90 % CI with t(3)), from eval_conv.py:

| difference | independent relaxation | start from the F state |
|---|---|---|
| dx 1.5 − dx 0.75 | −1.17 ± 0.17 (−1.57 … −0.77) | −0.96 ± 0.67 (−2.53 … +0.62) |
| dx 3 − dx 1.5 | −0.79 ± 0.31 (−1.51 … −0.07) | −0.98 ± 0.24 (−1.55 … −0.41) |
| dx 3 − dx 0.75 | −1.96 (unpaired Welch, −3.09 … −0.84) | −1.93 ± 0.50 (−3.10 … −0.77) |

As factors on P: dx 1.5 gives × 0.31 of dx 0.75 (independent), dx 3 gives × 0.45 of dx 1.5 and × 0.14 of dx 0.75. The
cycle noise is 0.038 (pooled SD of ln w_loop, dof 40). The SD over the seeds (0.54 … 0.96 on dx 1.5 and dx 0.75)
dominates the errors. Transferred against independent start: dx 3 +0.03 ± 0.05, dx 1.5 +0.21 ± 0.51.

**Reading (pre-registered): stop.** SD_T = 1.34 ≥ 0.10 (SD_T: SD of the paired difference dx 1.5 − dx 0.75 of the T
arm). The transferred state does not remove the scatter: stop and re-scope (PLAN 4.10). With n 3 (before the rerun)
the reading was the same (SD_T 1.64, runs/conv_b0_report_n19.md).

**Check of the measuring method** (exploratory, after the reading; go of Chris; eval_b0b_check.py →
results/b0b_check.md / .json):
- All 20 runs: 5 cycles, no event, a steady tail of 3 cycles, b_peak / b_target 0.96 … 1.01.
- dx 3: periodic loops in all 8 runs (closure ≤ 2·10⁻⁴, P equal within 0.1 % over the tail, RN ≤ 0.001, RF ≤ 0.014),
  w_dis / w_loop 1.000 … 1.002.
- dx 1.5: periodic in 7 of 8 runs (RN ≤ 0.009, w_dis / w_loop 0.98 … 1.01). In 4 of them RF is 0.15 … 0.23 with RN
  ≤ 0.009: a jump that repeats in every cycle at the same phase. T15 s942 is not periodic (RN 0.19, closure up to
  2·10⁻², P 7.5 … 9.5 W/cm³).
- dx 0.75: 3 of 4 seeds are not periodic (RN 0.037 … 0.19, closure 2·10⁻³ … 9·10⁻³, cycle SD of P 5 … 7 %). These 3
  seeds have the high values of P (7.1 / 12.3 / 36.1 W/cm³), in the order of their RN. The periodic seed s944 (RN
  0.012, RF 0.025) has 3.8 W/cm³, still 2.6 × its dx 3 values (1.39 / 1.48 W/cm³). n 4: an observation, no effect
  claim.
- **w_dis is not reliable on dx 0.75.** F s944 is periodic (closure 7·10⁻⁵), and still w_dis / w_loop is 1.071 …
  1.073 in all tail cycles. A closed loop does not change the energy of the state, so w_dis carries an error of about
  7 % there, also with atol 1e-6 (7.88 measured 1.01 … 1.02 at d/l_ex 48). Thus the values 0.86 … 1.16 of the other
  F seeds are not shown to be changes of the energy of the state. The cause is open. The reading uses w_loop and stays
  valid; the calm check of the F runs (tolerance 20 %) carries this error.
- Neighbour pairs above 60° in the relaxed state (structures that the grid does not resolve): 18 … 26 on dx 3, 0 … 11
  on dx 1.5, 0 … 23 on dx 0.75. The relaxed energy falls with the cell size: E / K_d 0.00212 … 0.00226 on dx 3,
  0.00182 … 0.00190 on dx 1.5, 0.00146 … 0.00157 on dx 0.75.
- Reading of the check: on dx 3 and dx 1.5 the measuring method holds (all runs calm and on target, energy balance
  within 0.2 % and 2 %). On dx 0.75, w_loop and thus P hold, w_dis does not. Jumps occur at 9 mT in 3 of 4 seeds on
  dx 0.75, in some runs on dx 1.5 (mostly repeating ones) and in none on dx 3. The jumps do not explain the whole
  mesh factor: the periodic seed s944 also has 2.6 × the loss of dx 3. The F values carry an extra uncertainty: the
  stop rule ended them after 5 cycles, with 3 tail cycles of a motion that does not repeat (cycle SD 5 … 7 %, small
  against the seed SD of 0.96 in ln P).
- Earlier mesh checks, for comparison (hints only):
  - d/l_ex 96, the same microstructure, old pairs with a relaxation that did not converge (7.64): at 50 mT dx 1 gave
    × 1.64 and dx 3 × 0.33 of dx 1.5 (n 3, the 90 % CIs include 0).
  - K3 at d/l_ex 212, another microstructure (7.63): the dx 3 grid cannot hold the states of dx 1.5; dx 3 gives × 0.61
    of dx 1.5 at 9 mT.
  - d/l_ex 300, the same microstructure (7.64): dx 3 gives × 1.48 of dx 1.5 at 9 mT, but the 90 % CI touches 0.
- The B0b factors hold for this point only (base microstructure, d/l_ex 96, 9 mT). They are no rule for the design.

**Statements that this section replaces or limits:**
- 7.100, "the 60 MHz values at 9 mT rest on cycles that have not settled": for the small particle at 9 mT the C′ value
  is a valid mean of a statistically steady irregular state (scatter of P 4 … 12 % SD over cycles). Its C′ channel
  shares at 9 mT possibly depend on the cycle window (candidate). Still open: the small particle at 20 and 35 mT (7.99
  shows irregular motion up to 35 mT) and the large particle at 60 MHz (P drift, 7.100). Thus the factor 2 and the
  sign change of the A effect at 60 MHz (7.98, 7.100) stay open.
- 7.101, candidate cause "integrator error of w_dis" for Δe / w_loop below 0: confirmed at 9 mT and in part at 20 mT;
  not supported at 50 and 70 mT; all at K1 only (n 3).

**Open decision for Chris:** the re-scope after the B0b reading "stop" (PLAN 4.10).

**Decision of Chris (2026-10-14):** no further evaluation of the wall area from VTI snapshots. The amount of wall does
not set β (7.19, 7.46); the texture now comes from the per-cell accumulators.

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

(Pilot error budget; for the Sobol design superseded by PLAN 4.0: Gamma GLM
with bootstrap, decision rules with Holm and TOST ±0.15, expected SE ≈ 0.06.)

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

Sobol design v2 (PLAN 4.0, prepared 7.30): 96 + 8 points at α 0.1 over
L_eff, r_p, φ and d/l_ex at dx 3, 4 corner groups × 3 seeds, 8 centre
replicates at α 0.01 (124 runs); outputs β_core (3-point windows, all) and
ln E[w] at 10 / 50 / 100 mT, the local-response measures from the snapshots
(9 / 50 / 150 mT) and the per-cell accumulators (all amplitudes); main fit
Gamma GLM on w, OLS and Huber as cross-checks; decision rules fixed before
the start (PLAN 4.0, PROTOKOLL 7.33). The earlier
plans (one-factor, factorial, 128 + 8) are replaced.

## 10. Known limitations

- 4 particles per box, identical spheres, FCC order, no contacts.
- T = 0 (no thermal activation). No eddy currents: estimate for a sphere
  (π f B̂ d)²/(20 ρ) at 30 MHz, d = 1 µm, ρ ≈ 1.2 µΩm: 3.7 W/cm³ of particle
  volume for B̂ = 0.10 T in the particle, ≈ 8 W/cm³ for B̂ ≈ 0.15 T in the
  particle (≈ B_core/φ at 100 mT); per core volume (× φ) ≈ 5 W/cm³, against
  576 W/cm³ in the model (core volume, 100 mT): small (audit 7.32, 7.33).
- Herzer length only up to a factor of order 1.
- K_eff is set, not computed: the Herzer averaging from the grains to the
  cubes is assumed, not simulated (section 2.1).
- No anisotropy on intermediate scales (≈ 300 nm), by choice.
- Vortex cores and Bloch points not resolved (flagged). The loss at small B
  sits at the vortex-core region; the absolute loss at dx 3 was higher
  than at dx 1.5 in 3 of 4 seeds (+21 … +43 % in the mean, not
  significant), the direction to the converged value is not known
  (7.25, 7.28).
- Staircase surface (the volume-fraction model is not usable, 7.29).
- One base point for the mechanism studies (7.17 … 7.29); f fixed at 30 MHz
  (the loss per cycle depends on f, 7.26).

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

**Status 2026-10-06:** all checks before the Sobol design are done
(7.25 … 7.29). Open: (1) Chris's go for PLAN 4.0 (audit 3 handled, 7.33);
(2) credit top-up to ≈ 45 $ (now ≈ 11 $). Decided: α 0.1 main + 8 centre
replicates at α 0.01 (item 3 below extended); mesh dx 3 with the limits
of 7.28 and the mesh-check rule of PLAN 4.0 (item 4 superseded); staircase
surface; paired α comparison done (7.31).

**Kept in mind (Chris, 2026-10-04): how does the model reach the high β of
measured cores?** In the model μ = 8 is fixed by the stray field of the
spheres (the particles are intrinsically very soft) and the response is
linear. (Corrected after the audit, 7.32: this does not imply β ≤ 2. With
β = 2 + d ln sin δ/d ln B, β > 2 is possible at constant μ; the demag
shearing does not change ∮H dM at a given ΔM, so Rayleigh hysteresis
inside the particles survives as β → 3, but its size scales with
(ΔM/χ_i)³ and is small for very soft particles.) A larger β needs a loss
angle that grows with B: particles whose own susceptibility is not ≫ 1/N, e.g. a particle
anisotropy K_p ≈ 90 kJ/m³ (r_p ≈ 14, estimate); the Sobol range r_p ≤ 3 is
probably still linear. Options: (A) extend r_p to 10 … 15 (local wall width
≈ 1 cell at dx 3: needs dx 1.5, ≈ 30× cost); (B) pre-test r_p 3 and 6 at
1 µm, dx 3, 3 seeds each (≈ 20 V100-h ≈ 2 $); (C) compare with Chris's
measured β (frequency, B range, material, particle size, core μ).
Recommended: B and C. Not started.


**Status 2026-10-04:** the target is the total loss at 30 MHz with the
model α (Chris 2026-10-03); the hysteresis-method question below is
superseded. Step 2 is decided (PLAN 4.2, option 3) but **on hold** until
the checks of 7.24 (mesh300, 9 mT frequency test, phase-resolved snapshots
and ring-down test) are done (Chris 2026-10-04).

Open 2026-10-03 (superseded): the method for the hysteresis loss (7.9). Options: (A)
quasi-static loops by energy minimization (no damping loss by
construction; validation against LLG at 9 … 20 mT); (B) LLG at a much lower
f (cost ∝ 1/f, ≈ 100× for 100 mT); (C) smaller α and f. Step 2 of PLAN.md
waits for this decision.

(Superseded by PLAN 4.0 / 7.33:) Decisions 2026-10-02 for the production (PLAN.md): D1 128 + 8 Sobol points
first, extension to 256 if necessary; D2 7 amplitudes (9 … 150 mT); D3
d/l_ex 212 … 300 with d/L_eff ≥ 7 and a residual check vs d/L_eff; D4
frequency test at the base point and at the corner L_eff 30, r_p 3.


Older open points (superseded by PLAN 4.0 / 7.33; kept as record):

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
