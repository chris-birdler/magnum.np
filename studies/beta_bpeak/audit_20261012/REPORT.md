# Audit of the β(B_peak) powder study (2026-10-12)

Auditor: context-free agent. Scope: all local run data in `runs/`, the analysis code, the result files, the deck `slide_decks/beta_powder_sobol`. Study state: PROTOKOLL 7.65, commit 080bd8f. I changed no existing file. Scripts and outputs: `audit_20261012/` (list in F.3).

(Saved by the main session from the agent's hand-back text; the agent could not write this file itself. Two findings were re-checked by the main session with its own fit: the A × K_eff interaction and the α 0.01 loops at 9 mT. Both hold.)

Terms:
- **Run**: one simulation folder. **Stage**: one amplitude of a run. **Measured cycles**: cycles 2 … end of a stage.
- **Event cycle**: w_loop ≤ 0 or |w_dis/w_loop − 1| > 0.20. **Calm tail**: the measured cycles after the last event cycle.
- **Main rule**: a calm tail of ≥ 2 cycles, and a finished rerun X_mc2 replaces X (7.54, 7.65).
- **P**: loss density at 30 MHz, P = f·K_d·w, in W/cm³ per core volume.
- **Range effect**: prediction at the factor maximum minus at its minimum, the other factors averaged over the 116 design points.
- **Signs**: physical. "A 10 → 20 pJ/m" = d/l_ex 300 → 212 at d = 1 µm.
- **dx 3 / dx 1.5**: cell size in units of l_ex. **Mesh error**: value(dx 3) − value(dx 1.5) for the same seed.
- **Errors**: SE = standard error over seeds; CI = confidence interval; n = number of independent runs; n.s. = not significant.

## (A) Executive summary

**A.1 Data currency.** 275 current folders: 116 design runs (incl. 13 reruns), 8 α 0.01, 18 conditioning, 27 mesh300 jobs, 75 K3 jobs, 18 same-state jobs.
- All have DONE, a converged relaxation and a job line equal to the config; stdout.log equals summary.json.
- Every stage has a calm tail ≥ 2 (the 8 α 0.01 runs are exempt by design).
- The 19 reruns are bit-identical to their originals up to the first extended stage.
- The M212 jobs are pending (0 local).
**The deck is stale.** 41 items checked, 36 need a change; three statements are wrong in substance (C.2).

**A.2 Levers for β (9 … 150 mT).** dx 3, 116 runs, Gamma GLM with 15 terms, bootstrap 1000. All point estimates of results/sobol_analysis_steady2_mc2.json reproduce to 4 digits.

| Lever (low → high) | Δβ (95 % CI) | Class | Robust? |
|---|---|---|---|
| A 10 → 20 pJ/m at d = 1 µm (= d 1.0 → 0.71 µm at A 10 pJ/m) | **+0.57 (+0.48 … +0.66)** | lever | yes: +0.53 … +0.73 in all 9 variants |
| K_eff 1.0 → 6.2 kJ/m³ | **+0.15 (+0.07 … +0.23)** | lever | yes in all steady variants (p ≤ 0.002); not in the pre-registered all-cycle fit (+0.10, p 0.03) |
| K_σ = r_p·K_eff, r_p 0 → 3 | +0.14 (+0.03 … +0.26) | at the decision limit | p 0.011 … 0.019 over 5 bootstrap seeds vs Holm limit 0.010; "lever" in 3 of 9 variants |
| φ 0.55 → 0.69 | −0.06 (−0.18 … +0.02) | not determined | 90 % CI edge −0.150 … −0.161 over 5 seeds; never inside ±0.15 |
| α 0.1 → 0.01 (centre, n 8) | +0.17 ± 0.11 (all cycles) / −0.21 ± 0.15 (calm cycles ≥ 1) | not determined | sign depends on the rule |

**New: strong A × K_eff interaction.**
- The A lever is +0.93 (+0.80 … +1.05) at K_eff 1.0 kJ/m³ and +0.16 (−0.03 … +0.35) at K_eff 6.2 kJ/m³; it falls by 0.78 (0.52 … 1.01).
- The K_eff lever is +0.53 (+0.38 … +0.66) at A 10 pJ/m and **−0.24 (−0.40 … −0.09) at A 20 pJ/m**.
- Lowest β of the design (≈ 1.6): low A and low K_eff. Highest (≈ 2.5): high A and low K_eff (figures/map_beta_A_Keff.png).
- With the Herzer coupling K_eff ∝ A⁻³, A × 2 gives Δβ +0.40 (+0.28 … +0.53) (slight extrapolation: the design covers K_eff × 6.25, A × 2 needs × 8).

**A.3 Partial ranges** (Δβ, 95 % CI):

| Lever | 9 … 35 mT | 35 … 70 mT | 70 … 150 mT |
|---|---|---|---|
| A 10 → 20 pJ/m | +0.31 (+0.12 … +0.50) | +1.00 (+0.72 … +1.28) | +0.35 (+0.13 … +0.59) |
| K_eff 1.0 → 6.2 kJ/m³ | +0.14 (−0.02 … +0.29) | +0.12 (−0.13 … +0.33) | **+0.31 (+0.10 … +0.49)** |
| r_p 0 → 3 | n.s. | +0.29 (+0.02 … +0.55) | n.s. |

Both main levers act where dx 3 has the largest β error (A.5).

**A.4 β from power-law fits per run** (7 amplitudes, realised B):
- β = 2.12 (sd 0.25, SE 0.023, n 116), range 1.59 … 2.60; β > 2 in 77 (nominal-B fit, as the design analysis) / 78 (realised-B fit) of 116.
- rms residual of ln P 0.13 (median) against a loop error of ≈ 1 %: a single power law describes the trend but is not exact.
- Local β ≈ 2.15 from 9 to 70 mT, then 2.07 (50 … 100 mT) and 1.98 (70 … 150 mT); paired β(70 … 150) − β(9 … 35) = −0.17 ± 0.04 (p 8·10⁻⁵).
- Mean curvature dβ/d ln B = −0.08 ± 0.02.
- μ_r is constant over B to 0.1 %, so **β − 2 = d ln(sin δ)/d ln B exactly**: every lever acts through the loss angle δ.

**A.5 Mesh: dx 3 vs dx 1.5.** Base point d/l_ex 300 (n 4); K3 corner d/l_ex 212 (n 6).

| Output | Needs dx 1.5? | Error to assume with dx 3 |
|---|---|---|
| P level | **yes** | Sign depends on the microstructure. Base ln P +0.39 ± 0.20 / +0.34 ± 0.16 / +0.22 ± 0.07 (× 1.48 / 1.40 / 1.25 at 9 / 50 / 100 mT). K3 −0.47 ± 0.09 / −0.36 ± 0.13 / −0.14 ± 0.09 (× 0.63 / 0.70 / 0.87). K3 − base: −0.86 ± 0.22, −0.70 ± 0.21, −0.36 ± 0.11, each p ≈ 0.01. Assume ±0.5 in ln P (× 0.6 … 1.6) at 9 … 70 mT and ±0.25 at 100 mT, sign unknown (estimate). |
| β full range | no, with an error bar | Base −0.06 ± 0.05 (95 % CI −0.23 … +0.12). K3 +0.07 ± 0.03 (9 … 100 mT), +0.12 ± 0.03 (9/50/100). Assume ±0.2 (estimate). |
| β below 50 mT | no | K3 9 … 35 mT −0.01 ± 0.11. 9 → 50 mT: base −0.02 ± 0.03, K3 +0.05 ± 0.05. |
| β above 50 mT | **yes** | K3 50 … 100 mT **+0.35 ± 0.11** (CI +0.06 … +0.63). Base 50 → 100 mT −0.17 ± 0.15 (other sign, n.s.). |
| A lever on β (K3 − base pair) | partly | +0.54 ± 0.11 (dx 3) vs +0.36 ± 0.07 (dx 1.5): mesh share +0.18 ± 0.13, ratio 0.67 ± 0.18. Design lever on the fine mesh ≈ **+0.38 ± 0.11 (estimate)**. Above 50 mT the pair effect is mostly mesh-made (+0.76 → +0.25). |
| A lever on P | **yes** | At 50 mT the sign reverses: × 0.79 (dx 3) vs **× 1.59 (dx 1.5, p 0.034)**. The dx 3 crossing at ≈ 70 mT moves below 50 mT. |
| K_eff and r_p levers | unknown | K2 not run. Their effects come from > 35 mT, where the β mesh error reaches ±0.35. |
| Jump share RF | yes | The 35 … 70 mT RF window is mesh-made (6 of 6 seeds at 70 and 100 mT). |

Caveats: the K3 − base pair also changes L_eff, r_p and φ. The dx 3 model predicts Δβ +0.43 (+0.18 … +0.71) for this pair; 63 % of it comes from d/l_ex. Nobody has shown that dx 1.5 is converged.

## (B) Inventory and data currency (task 1)

Script inventory.py → inventory.csv, inventory.md.

**B.1 Status of the 441 folders** (with the 22 files and logs, 463 entries).

| Group (PROTOKOLL) | Folders | Status |
|---|---|---|
| Sobol main S000 … S095 (7.44 / 7.54 / 7.65) | 96 | 86 current; 10 replaced by reruns |
| Centre C0 … C7 | 8 | current |
| Corners K1_0 … K4_2 | 12 | 9 current; 3 replaced (K1_0, K4_0, K4_1) |
| Design reruns X_mc2 | 13 | current (in-process relaxation) |
| α 0.01 C0_a001 … C7_a001 (7.48 / 7.62) | 8 | current; 23 stages with calm tail < 2 (all cycles used, by decision) |
| Conditioning cond_C0 … C5 × b9/b20/b35 | 18 | current |
| mesh300 + reruns | 24 + 3 | 21 + 3 current; 3 replaced |
| K3 9/50/100 mT (7.52) | 36 | 35 current; MK3_dx15_s2034_b50 replaced |
| K3 fill-in 20/35/70 mT (7.57) | 36 | 34 current; MK3_dx15_s2034_b35, MK3_dx3_s2035_b35 replaced |
| K3 reruns | 3 | current |
| Same-state MK3cg (7.63) | 18 | current |
| Earlier converged studies: B300 8, PH 2, F9 6, V 2, A010c 3, X010/X100 6, R300f 12, R300v01 12 | 51 | valid for their sections (all-cycle rule) |
| R300f01_s921_b9 | 1 | stopped, no summary, not used |
| Withdrawn, relaxation not converged (7.13): A010/A030, F_base/F_corner, P1/P2/P3/R_dx15, R212, S212/S300_dx3f, M_L12/M_floor, S_L12/S_floor, F_do/F_si, I_*, D_*, C_dx* | 112 | withdrawn; old mesh pairs are hints only |
| bench_base | 1 | test, unfinished |
| _archive_sobol_v2_oldaxes | 1 (20 + 4 runs) | archived, old stress axes |
| X_mc2init (retired) | 0 local | renamed on the instance only; local S020_mc2, S026_mc2, M300_dx15_s922_b9_mc2 are the in-process variants (init_from None, bit-identical to the originals) |
| M212_* (PLAN 4.6) | 0 local, 24 jobs | **pending on the GPU instance** |

**B.2 Integrity checks**
- DONE and a converged relaxation in 275 of 275 current folders. Every job line equals its config.
- stdout.log equals summary.json; every cycle is printed.
- W_loop and W_dis recomputed from samples_*.csv match summary.json to 5·10⁻⁸ and 9·10⁻⁸ (relative).
- Every stage of the design, mesh300, K3, fill-in, same-state and conditioning runs has ≥ 2 calm cycles (confirms the 7.65 acceptance).
- Determinism: 19 of 19 reruns bit-identical up to the first extended stage; conditioning runs repeat C0 … C5 bit for bit (12 of 12 checks). Commits uniform per group; one GPU type.
- md5: all 4297 entries of the three manifests verify. **Not in any local manifest:** 55 new folders (13 + 6 reruns, 18 conditioning, 18 same-state) and 7 old test folders (bench_base, M_floor × 3, M_L12 × 3). (Main session, 2026-10-12: runs/MANIFEST_mc2b.md5 now covers these 63 folders, 506 files, verified.)

**B.3 Which result files are current**
- Current (12 Oct): sobol_analysis_steady2_mc2.*, k3_mesh.md, k3_beta_fits.md, k3_fine.md, mesh212.md (no M212 data yet), jumps_time.md, conditioning.md.
- Stale: deck_data.json (11 Oct, no reruns), k3_same_state.md (before the rerun MK3_dx15_s2034_b50_mc2).
- Superseded, not used by the deck: crossover.md, mesh300.md and older files.

## (C) Deck items to update (task 2)

The deck (11 Oct 16:45) reads deck_data.json (11 Oct 16:42), built from sobol_analysis_steady2.json (no reruns, 13 stages missing), interim K3 data and 3 fine-mesh repeats at the base point. A rebuild takes the _mc2 file and current mesh data automatically, except for the parts in C.3. Sources: deck_check.md, sobol_audit.md, mesh_audit.md. 41 items checked; items 3, 4, 9, 13 and 23 need no change.

**C.1 Numbers (deck → current)**

| # | Slide | Item | Deck | Current |
|---|---|---|---|---|
| 1 | 1 | date / status | 11 Oct, interim | 12 Oct, final (PLAN 4.5 done) |
| 2 | 5–13 | footnote "Interim: reruns … running" | present | remove |
| 3 | 5 | "every value rests on ≥ 2 steady loops" | false at build time | now true |
| 4 | 5 | loop scatter | 2.5 % | 2.5 % (n 783) |
| 5 | 6 | β range | 1.56 … 2.60 | 1.59 … 2.60 |
| 6 | 6 | β mean (sd) | 2.11 (0.25) | 2.12 (0.25) |
| 7 | 6 | β > 2 | 78 of 116 | 77 of 116 |
| 8 | 6, 15 | repeat scatter (20 repeats) | 0.14 | 0.13 (dof 15) |
| 9 | 6 | spread from parameters | 0.21 | 0.21 |
| 10 | 7, 16, 17 | A lever on β | +0.58 (+0.48 … +0.68) | +0.57 (+0.48 … +0.67) |
| 11 | 7, 9, 16, 17 | K_eff lever | +0.15 (+0.07 … +0.23) | +0.15 (+0.07 … +0.22) |
| 12 | 7, 16, 17 | stress lever | +0.16 (+0.03 … +0.27) | +0.14 (+0.02 … +0.26) |
| 13 | 7, 16, 17 | packing | −0.06 (−0.16 … +0.03) | unchanged |
| 14 | 7, 17 | P at 50 mT: A / K_eff / stress / packing | −26 / −6 / +4 / +6 % | −26 % (−40 … −9) / −6 / +4 / +7 % |
| 15 | 7 | fine-mesh A lever (callout, diamond) | +0.31 | +0.38 ± 0.11 (factor 0.67 ± 0.18) |
| 16 | 7 | diamond P at 50 mT | +48 % | +59 % |
| 17 | 7, 8, 15 | repeats of the mesh pair | 5 and 3 | 6 and 4 |
| 18 | 8 | Δβ of the pair, coarse / fine | interim | +0.54 ± 0.11 / +0.36 ± 0.07 |
| 19 | 8, 15 | "about half the size" | half | about two thirds |
| 20 | 9 | K_eff effect per window | interim | +0.14 (−0.01 … +0.28) / +0.12 (−0.14 … +0.33) / +0.31 (+0.10 … +0.52) |
| 21 | 10 | fine mesh at 9 mT | × 0.81 (p 0.36, 6/3) | × 0.82 (p 0.31, 6/4) |
| 22 | 10 | fine mesh at 50 mT | × 1.48 (p 0.09, 5/4) | × 1.59 (p 0.034, 6/4) |
| 23 | 10 | fine mesh at 100 mT | × 1.91 (p < 0.001) | unchanged (6/4) |
| 24 | 10 | coarse crossing, pair / all runs | ≈ 75 / ≈ 69 mT | ≈ 68 / ≈ 70 mT |
| 25 | 10, 16 | "probably at 50 mT (hint)" | p 0.09 | p 0.034 (one uncorrected test, confounded pair) |
| 26 | 11 | K3 β 9 … 100 mT, coarse / fine | 2.38 / 2.27 | 2.37 / 2.30 |
| 27 | 11 | paired difference | +0.11 ± 0.02 (n 4) | +0.07 ± 0.03 (n 6, p 0.057) |
| 28 | 11 | "above, coarse gives too high β" | — | only 50 … 100 mT resolved (+0.35 ± 0.11); 20 … 50 mT −0.16 ± 0.10; base 50 → 100 mT −0.17 ± 0.15 (other sign) |
| 29 | 13 | K3 repeats | 4 … 6 | 5 … 6 |
| 30 | 15 | steady loops "running" | running | done |
| 31 | 15 | long-run test "running" | running | done: drift ≤ 3 % (90 % CI), Δβ −0.003 ± 0.018 / −0.004 ± 0.021 (n 6) |
| 32 | 15 | damping row | +0.17 ± 0.11, one pair +0.85 | +0.17 ± 0.11, max +0.83 (all cycles); −0.21 ± 0.15 (calm): not determined |
| 33 | 16 | next step (1) | open | done; add: M212 running |
| 34 | 17 | stress p footnote | 0.0069 vs 0.015 | 0.0145 (main), 0.0069 (no reruns), 0.015 (≥ 1) |
| 35 | 17 | packing 90 % edge | −0.144 vs −0.152 | −0.150 (main), −0.144, −0.152 |
| 36 | 18 | "× 0.60 … 0.85, 20 … 150 mT" | — | × 0.60 … 0.84 (all cycles); calm: × 0.63 … 0.82 at 35 … 150 mT (n 7) |

**C.2 Statements that are wrong or not supported**

| # | Slide | Statement | Finding | Correct statement |
|---|---|---|---|---|
| 37 | 9, 15, 16 | "Higher K_eff raises β" | holds only at low A; at A 20 pJ/m: −0.24 (−0.40 … −0.09) | "Higher K_eff raises β at low A (+0.53 at 10 pJ/m) and lowers it at A 20 pJ/m (−0.24)." |
| 38 | 7, 8, 16 | A lever as one number (+0.57) | +0.93 at K_eff 1.0, +0.16 (n.s.) at 6.2 kJ/m³ | show the K_eff dependence (map figure) |
| 39 | 18 | "× 0.34 at 9 mT" | 6 of 8 α 0.01 seeds have no calm loop at 9 mT; 8 of 8 have fewer than 2; 5 of 8 at 20 mT. Loop areas negative in 5 of 8 (balance −180 … +194). | "At 9 mT the α 0.01 loops do not settle in 9 loops: no value." |
| 40 | 12 | "61 % of P from 10 % of the volume" | that is the texture-top 10 %; the p-top 10 % holds 69 % | "The 10 % most textured volume (walls, cores) carries 61 % of P at 9 mT." |
| 41 | 7, 16, 17 | r_p / φ "at the decision limit" | correct per 7.65, but a rebuild shows "not determined" (C.3) | keep, with the reason in E.7 |

Also check on rebuild: the slide 11 title (holds only at K3, one point); slide 13 R1 values (acc.json, events included: hint only).

**C.3 Code points of the deck build**
- make_deck.disp_class gives "at the decision limit" only when the ≥ 1 and ≥ 2 analyses disagree. With _mc2 both say "not determined". Fix: compare with steady2 (no reruns), or set the class explicitly per 7.65.
- INTERIM, DATE and the slide-15 "running" rows are hard-coded.
- "about half" is fixed text although a_fine_factor computes the ratio.
- jumps_block takes R1 from acc.json of the originals, events included: keep, but mark.

## (D) Detailed results

### D.1 Task 3: levers, independent re-analysis

Scripts: design_data.py (loader, 6 rule variants), sobol_audit.py (own IRLS Gamma GLM, scatter model, bootstrap), sobol_extra.py (interactions, bootstrap noise, flag-column check, mesh-point predictions).
Reproduction: per-run outputs equal the study CSV files to 10⁻¹⁴. Without reruns the 9 … 35 mT window differs by ≤ 0.003 (E.9). With the study bootstrap seed, CIs and p values are identical.

**Variants (β all):**

| Variant | A | K_eff | r_p | φ |
|---|---|---|---|---|
| main (≥ 2 calm, reruns, weights, 15 terms) | +0.572 (+0.479 … +0.663) lever | +0.150 (+0.065 … +0.227) lever | +0.145 (+0.027 … +0.258) n.d., p 0.0135 | −0.064 (−0.182 … +0.024) n.d. |
| no weights | +0.590 lever | +0.136 lever | +0.133 n.d. | −0.072 n.d. |
| 5 terms, weights | +0.727 lever | +0.156 lever | +0.188 lever (p 0.008) | −0.105 n.d. |
| 5 terms, no weights | +0.577 lever | +0.145 lever | +0.184 lever (p 0.009) | −0.103 n.d. |
| ≥ 2, no reruns | +0.577 lever | +0.152 lever | +0.157 lever (p 0.006) | −0.058 n.d. |
| ≥ 1, no reruns | +0.572 lever | +0.151 lever | +0.145 n.d. | −0.065 n.d. |
| all cycles (pre-registered, 115 runs) | +0.531 lever | +0.103 n.d. (p 0.032) | +0.135 n.d. | −0.035 n.d. |
| last 2 calm cycles | +0.597 lever | +0.138 lever | +0.140 n.d. | −0.069 n.d. |
| w_dis | +0.581 lever | +0.133 lever | +0.128 n.d. | −0.073 n.d. |

(n.d. = not determined.) The 5-term weighted fit gives +0.73 for A: the scatter weights down-weight high-K_eff runs (slope +0.82, p 0.008) and the A lever is larger at low K_eff. The range effect is a design average, not a material property.

**Effects on P** (dx 3, factors):

| Lever | 10 mT | 50 mT | 100 mT | 150 mT |
|---|---|---|---|---|
| A 10 → 20 pJ/m | × 0.33 (0.27 … 0.40) | × 0.73 (0.60 … 0.88) | × 1.21 (1.01 … 1.38) | × 1.32 (1.17 … 1.50) |
| K_eff 1.0 → 6.2 kJ/m³ | × 0.73 (0.62 … 0.87) | × 0.94 (0.78 … 1.12) | — | × 1.17 (1.06 … 1.28) |
| φ 0.55 → 0.69 | × 1.34 (1.12 … 1.68); × 1.07 without dilution | — | — | — |

The A effect on P changes sign at ≈ 70 mT.

**Conditional effects** (95 % bootstrap CI):

| Effect | At | β all | β 9 … 35 | β 35 … 70 | β 70 … 150 | ln P 10 mT | ln P 100 mT |
|---|---|---|---|---|---|---|---|
| A 10 → 20 | K_eff 1.0 | +0.93 (+0.80 … +1.05) | +0.69 | +1.37 | +0.71 | −1.73 | +0.41 |
| A 10 → 20 | K_eff 2.8 | +0.55 (+0.45 … +0.64) | +0.28 | +0.97 | +0.32 | −1.06 | +0.17 |
| A 10 → 20 | K_eff 6.2 | +0.16 (−0.03 … +0.35) | −0.13 | +0.57 | −0.07 | −0.40 | −0.07 |
| K_eff 1.0 → 6.2 | A 20 | −0.24 (−0.40 … −0.09) | −0.28 | −0.29 | −0.09 | +0.36 | −0.26 |
| K_eff 1.0 → 6.2 | A 10 | +0.53 (+0.38 … +0.66) | +0.55 | +0.51 | +0.69 | −0.97 | +0.22 |
| r_p 0 → 3 | A 20 / A 10 | +0.01 (−0.23 … +0.24) / +0.28 (+0.11 … +0.47) | | | | | |

- Model-free support: β slope over the d/l_ex range within each L_eff level: −0.18 ± 0.10 (L 12, n 27), −0.49 ± 0.16 (L 18, n 38), −0.77 ± 0.11 (L 24, n 24), −0.95 ± 0.07 (L 30, n 27) (OLS, HC3).
- The lnd × lnQ term raises R² from 0.53 to 0.66 (p < 10⁻⁴). The 15-term model leaves a residual sd of 0.14, equal to the repeat scatter (0.13).

**Bootstrap noise** (5 seeds × 1000): p(r_p) = 0.011 … 0.019 (normal), 0.014 … 0.022 (percentile); φ 90 % edge −0.150 … −0.161. Both decisions lie inside the Monte-Carlo noise.

**Real alloy** (K_eff ∝ A⁻³), A × 2: Δβ +0.40 (+0.28 … +0.53); ln P −0.75 at 10 mT, −0.23 (n.s.) at 50 mT, +0.21 at 100 mT.

**Damping.** With calm loops (≥ 1) the P ratio α 0.01 / α 0.1 is × 0.78 / 0.82 / 0.73 / 0.78 / 0.63 at 35 / 50 / 70 / 100 / 150 mT (n 7). At 9 and 20 mT there is no valid value (6 of 8 and 5 of 8 seeds have no calm loop). Δβ not determined.

### D.2 Task 4: power-law fits per run

Script beta_fits.py → beta_fits.md, beta_fits.csv, figures/beta_fits.png, figures/beta_A_by_Leff.png.

| Window (mT) | Mean ± SE (n 116) | sd | Share > 2 |
|---|---|---|---|
| 9 … 150 | 2.116 ± 0.023 | 0.250 | 0.67 |
| 9 … 35 | 2.150 ± 0.027 | 0.288 | 0.72 |
| 20 … 50 | 2.151 ± 0.041 | 0.446 | 0.62 |
| 35 … 70 | 2.163 ± 0.046 | 0.493 | 0.62 |
| 50 … 100 | 2.067 ± 0.041 | 0.440 | 0.52 |
| 70 … 150 | 1.980 ± 0.034 | 0.366 | 0.49 |

- rms residual of ln P median 0.134 (quartiles 0.087 … 0.175). χ²/dof against the loop error 125 … 25 000 (quartiles): the misfit is real structure of each run, not noise.
- Mean curvature −0.077 ± 0.023; 24 of 116 runs have a significant curvature (18 negative, 6 positive).
- A broken power law improves 35 of 116 runs; the best knot scatters over 20 … 100 mT (33 / 22 / 17 / 31 / 13 runs): no common break point.
- β(70 … 150) − β(9 … 35) = −0.17 ± 0.04; β(20 … 50) − β(9 … 35) = 0.00 ± 0.04.
- Factors (per-run OLS with lnd × lnQ, β 9 … 150 mT, R² 0.66): lnd −0.56 ± 0.05, lnQ +0.13 ± 0.04, r_p +0.14 ± 0.06, φ −0.06 ± 0.05, lnd × lnQ +0.41 ± 0.05.
- Groups: A ≤ 12.5 pJ/m 1.94 ± 0.04; A ≥ 15.8 pJ/m 2.34 ± 0.02; L_eff 30 high A 2.44 ± 0.03 (n 12); L_eff 30 low A 1.74 ± 0.05 (n 7). High-A runs peak at 35 … 70 mT (2.57 ± 0.07).
- r_p: β 9 … 35 mT residual −0.27 ± 0.10 at r_p < 0.3 (n 9), ≈ 0 above. Hint only (n 9; no run has r_p = 0 with the new axes).
- Code comparison: design β_all differs from my realised-B fit by ≤ 0.009. analyze.features β_pow weights stages by the number of kept cycles; differs by up to 0.078; not used in the design analysis.

### D.3 Task 5: mesh

Script mesh_audit.py → mesh_audit.md, figures/mesh_level.png.

| Point | B (mT) | P dx 1.5 (W/cm³) | P dx 3 | Δ ln P (n) | Factor |
|---|---|---|---|---|---|
| base | 9 / 50 / 100 | 5.2 / 134 / 496 | 7.7 / 188 / 619 | +0.39 ± 0.20 / +0.34 ± 0.16 / +0.22 ± 0.07 (4) | × 1.48 / 1.40 / 1.25 |
| K3 | 9 / 20 / 35 | 4.3 / 22.0 / 84.4 | 2.7 / 17.0 / 51.0 | −0.47 ± 0.09 / −0.26 ± 0.13 / −0.50 ± 0.17 (6) | × 0.63 / 0.77 / 0.60 |
| K3 | 50 / 70 / 100 | 213 / 503 / 948 | 148 / 335 / 826 | −0.36 ± 0.13 / −0.41 ± 0.05 / −0.14 ± 0.09 (6) | × 0.70 / 0.67 / 0.87 |

Old pairs (relaxation not converged, hints) at d/l_ex 212 with the base microstructure: −0.36 ± 0.34 (9 mT), −0.06 ± 0.09 (50 mT) (n 4); at d/l_ex 96 −1.0 … −1.1. The sign may follow the particle size; M212 tests this.

| Point | Range (mT) | β dx 3 | β dx 1.5 | Δβ (n, 95 % CI) |
|---|---|---|---|---|
| base | 9/50/100 | 1.83 ± 0.10 | 1.89 ± 0.05 | −0.06 ± 0.05 (4, −0.23 … +0.12) |
| base | 9 → 50 / 50 → 100 | 1.86 / 1.72 | 1.89 / 1.89 | −0.02 ± 0.03 / −0.17 ± 0.15 |
| K3 | 9 … 100 | 2.37 ± 0.02 | 2.30 ± 0.04 | +0.07 ± 0.03 (6, −0.00 … +0.14) |
| K3 | 9/50/100 | 2.37 ± 0.03 | 2.25 ± 0.04 | +0.12 ± 0.03 (+0.04 … +0.19) |
| K3 | 9 … 35 / 20 … 50 | 2.18 / 2.32 | 2.19 / 2.48 | −0.01 ± 0.11 / −0.16 ± 0.10 |
| K3 | 35 … 70 / 50 … 100 | 2.72 / 2.49 | 2.57 / 2.14 | +0.15 ± 0.23 / **+0.35 ± 0.11** |

2-point segments are noisy (K3 20 → 35 mT −0.45 ± 0.11; 70 → 100 mT +0.75 ± 0.14): use windows of ≥ 3 amplitudes.

**Mesh share of the levers**
- A on β: pair +0.54 ± 0.11 (dx 3) vs +0.36 ± 0.07 (dx 1.5). Design model predicts +0.43 (+0.18 … +0.71) for this pair (d/l_ex +0.27, r_p +0.10, L_eff +0.06, φ 0). Fine-mesh design lever ≈ +0.38 ± 0.11 (estimate). The pair lies at L_eff 12 … 18, where the lever is small; no fine-mesh pair at L_eff 24 … 30.
- A on P: × 0.35 / 0.79 / 1.34 (dx 3) vs × 0.82 (n.s.) / 1.59 / 1.91 (dx 1.5) at 9 / 50 / 100 mT. On the fine mesh high A has more loss from ≈ 50 mT on; any crossing lies below 50 mT and is not resolved.
- K_eff and r_p: no data. The β mesh error above 35 mT (+0.35 at K3, −0.17 at base) is the same size as these levers.
- Same-state test (n 6): ln P(cg / dx 3 own) = −0.03 ± 0.09 / −0.22 ± 0.16 / −0.12 ± 0.10; ln P(cg / dx 1.5) = −0.49 ± 0.16 / −0.58 ± 0.14 / −0.26 ± 0.08. First loop at 9 mT releases 3300 … 6800 × the loop loss. Fine states are not stable on the dx 3 grid. (PROTOKOLL 7.63 gives −0.47 at 50 mT: pre-rerun.)
- Relaxed energies E/K_d: base dx 3 −0.00078 … −0.00088, dx 1.5 −0.00097 … −0.00098; K3 dx 3 −0.00086 … −0.00101, dx 1.5 −0.00114 … −0.00118.

### D.4 Task 6: time curves and field data

Script timecurves.py → timecurves.md, figures/tc_loopshape.png.

1. **Linear core.** μ_r = 6.9 (median; quartiles 5.9 … 8.2); B slope −0.0008 ± 0.0015. μ_r depends only on φ (ln μ_r +0.74 over the φ range; ≤ 0.03 for the other factors; R² 0.98). |M3|/|M1| = 3 … 6·10⁻⁴, |M2|/|M1| = 5 … 8·10⁻⁴.
2. **The loss angle carries β.** sin δ = 0.0059 (9 mT) → 0.0083 (70 mT) → 0.0079 (150 mT); δ 0.33 … 0.48°. π·μ0·|H1|·|M1|·sin δ equals W_loop to 0.25 %. β − 2 = d ln sin δ/d ln B (r = 1.000; μ_r term ≤ 0.003). Low A gives a large loss angle at small B that does not grow (lnd effect on ln sin δ: +1.09 ± 0.11 at 9 mT, −0.26 ± 0.07 at 150 mT).
3. **Where the loss is** (acc.json; core = 12.5 % of the volume, shell = 49 %):

| Region | 9 mT | 50 mT | 150 mT |
|---|---|---|---|
| texture-top 10 % of the volume | 0.61 | 0.69 | 0.55 |
| p-top 10 % of the volume | 0.69 | 0.76 | 0.73 |
| core r/R < 0.5 | 0.38 | 0.27 | 0.18 |
| shell r/R ≥ 0.8 | 0.22 | 0.39 | 0.50 |

   Supports 7.23 and 7.46.
4. **R1 (share of the local fundamental).** Snapshots of the last 2 calm cycles: 0.92 / 0.51 / 0.40 at 9 / 50 / 150 mT; acc.json (events included) 0.86 / 0.50 / 0.40. Gradual decline. The R1 < 0.6 onset moves × 0.48 (0.37 … 0.60) from A 20 to A 10 (≈ 40 → 26 mT, p < 0.001); × 1.27 with K_eff (p 0.053). On the fine mesh at K3, R1 differs by ≤ 0.08.
5. **Switched volume.** V_sw/V = 0.65 / 3.3 / 9.8 % per half cycle; exponent 0.98 (0.68 … 1.13); correlates with β (r = 0.91). In the factor model it takes the A lever from −0.56 to −0.13 and the interaction from +0.41 to +0.06 (R² 0.66 → 0.90). Partly by construction. Source: domains.json of the originals, as in the study.
6. **New: per-particle loss.** For r_p > 0.5, particles with the stress axis at 37° dissipate × 1.45 (9 mT) … × 1.84 (150 mT) more per volume than those at 79° (n 100), at equal fundamental susceptibility (× 1.02). No further growth with r_p above ≈ 0.3; at r_p 0.02 … 0.26, ln Q is already +0.24 … +0.34 (n 9). Hypothesis: K_σ exceeds the particle mean of the Herzer anisotropy (K_eff/√N, N = 185 … 8200 Herzer cubes per particle) already at r_p ≈ 0.1. Hint; no r_p = 0 run with the new axes exists.
7. **Events** (cycles 2 … 6, 116 originals): 173 release energy, 2 take it up.
8. **Jump measures:** RF/BF code reviewed; numbers taken from k3_fine.md and jumps_time.md, not recomputed. They support 7.57 and 7.60.

Not re-verified: the event-precursor claim (7.55 T5), the ring-down and frequency tests, the cut images.

**Check of the earlier interpretations**

| Claim | Verdict |
|---|---|
| 7.22 linear response | holds for the whole design |
| 7.22 / 7.24 identity for β | holds |
| 7.23 / 7.27 / 7.46 loss at cores and walls | holds |
| 7.46 switched volume | holds; extends to the interaction |
| 7.47 / 7.62 crossing ≈ 69 mT | holds at dx 3 (≈ 70 mT), not on the fine mesh |
| 7.50 R1 onset follows A | holds at dx 3 |
| 7.51 / 7.57 RF window | holds (mesh-made) |
| 7.53 / 7.55 events transient | holds |
| 7.48 × 0.34 at 9 mT and R(9 mT) | **not supported** |
| 7.63 same-state | conclusion holds; numbers pre-rerun |
| 7.65 decision limits | holds |
| Deck "higher K_eff raises β" | holds only at low A |

### D.5 Task 7: conditioning and convergence

Script cond_audit.py → cond_audit.md.
- eval_cond.py follows PLAN 4.5 D. My recomputation gives the same drifts: +0.001 ± 0.013, −0.001 ± 0.005, −0.004 ± 0.009 (n 6): the design value stands. 6 events, all in cycles 2 … 4.
- Trend of ln P over cycles 5 … 29 per 100 cycles: +0.018 ± 0.009 at 9 mT (95 % CI −0.006 … +0.042), −0.001 ± 0.009 at 20 mT (n.s.), +0.014 ± 0.021 at 35 mT (n.s.). One single run (C0, 9 mT) rises +5.7 ± 1.3 % per 100 cycles.
- Closure falls to 10⁻⁸ … 10⁻⁶ in 12 of 18 runs; stays at 3·10⁻⁴ … 3·10⁻³ in 5 runs (C1_b20, C2_b20, C0_b35, C1_b35, C4_b35) and at 9·10⁻⁵ in C3_b35. Those runs reach no periodic orbit, but their P does not drift.
- Design runs: slope of ln w over the calm tails −0.003 … +0.002 per cycle at all amplitudes; none significant.
- Limit: only the centre point tested. Per PLAN 4.5, 11 of the 13 rerun runs have A ≤ 15 pJ/m and 6 of 13 have L_eff 30.
- Determinism confirmed (B.2).

## (E) Errors and weak points

| # | Where | Finding | Effect |
|---|---|---|---|
| E.1 | deck | pre-7.65 data, "interim" text | 36 items (C.1) |
| E.2 | deck 9 / 15 / 16 | ignores the A × K_eff interaction | "K_eff raises β" wrong at high A |
| E.3 | 7.48, deck 18 | Owner rule: all cycles for the damping test. Finding: at 9 mT the α 0.01 loop areas are negative or erratic, so the mean is undefined. | drop the 9 mT ratio and R(9 mT) |
| E.4 | make_deck.disp_class | cannot show "decision limit" with _mc2 | rebuild mislabels r_p and φ |
| E.5 | analyze_design, steady mode | "flagged" computed on kept cycles → 0 flagged; "without flagged" column = unweighted main fit (checked identical) | label wrong |
| E.6 | alpha_test | nan CI printed as "mixed" | label wrong |
| E.7 | classify | normal p with bootstrap SE; r_p and φ inside the Monte-Carlo noise | correctly "decision limit" in 7.65 |
| E.8 | report | "227 of 4060" assumes 5 cycles per stage; actual 4081 | label |
| E.9 | to_nominal | 20 mT value clamped when 9 mT missing | ≤ 0.003 (variant without reruns only) |
| E.10 | deck 12 | "10 % of the volume" wording | wording |
| E.11 | PROTOKOLL §0 | still the 2026-10-06 state | update |
| E.12 | manifests | 55 new + 7 old test folders not in a manifest | owner rule not met (main session: fixed, MANIFEST_mc2b.md5) |
| E.13 | 7.57 / 7.63 / k3_same_state.md | pre-rerun numbers (e.g. P 197 → 213 W/cm³) | conclusions hold |
| E.14 | "the" lever | design average depends on weights under the interaction | report conditional effects |
| E.15 | acc.json | sums all cycles incl. events; rerun snapshots stay from originals | ≤ 0.05 bias of R1 |
| E.16 | features β_pow | weights ∝ number of kept cycles | ≤ 0.08; not used |
| E.17 | K3 − base "A effect" | also changes L_eff, r_p, φ | state it at every use |

No error found in: the steady-tail code, skip logic and rerun mapping; the GLM, bootstrap, Holm and TOST; eval_k3, eval_k3beta, eval_k3fill, eval_cond; the unit chain (P = f·K_d·w; K_d 895.2 kJ/m³; A at d = 1 µm); the loop integrals.

## (F) Open questions and next checks

1. **M212** (running; ≈ 39 … 87 V100-h): is the sign of the P mesh error set by size or microstructure?
2. **Fine-mesh pair at L_eff 30** (where the A lever is largest): K1-like point d/l_ex 212 and a point at 300, 6 seeds each, 9/50/100 mT. ≈ 13 + 25 V100-h per seed pair ≈ 230 V100-h ≈ 23 $ (estimate). Also gives a pure A contrast.
3. **K2 mesh check** (K_eff, r_p): ≈ 250 V100-h ≈ 25 $ (PLAN estimate).
4. **Report the conditional levers.** No cost.
5. **r_p step test:** r_p 0 / 0.1 / 0.3 / 1, 3 seeds, dx 3; ≈ 30 V100-h ≈ 3 $ (estimate).
6. **α 0.01 at 9 … 20 mT:** ≥ 30 cycles needed; ≈ 5 V100-h (estimate), or accept "not determined".
7. **dx 1.5 convergence:** one dx 1 seed at d/l_ex 212, ≈ 100 V100-h (estimate). Only after items 1 … 3.
8. **md5 manifest** for the folders of E.12. Local, no cost. (Done by the main session.)

**F.3 Files** (in audit_20261012/):
- audit_lib.py, design_data.py, inventory.py → inventory.csv, inventory.md
- sobol_audit.py → sobol_audit.md/.json, sobol_audit_boot.npz; sobol_extra.py → sobol_extra.md/.json
- beta_fits.py → beta_fits.md/.csv; mesh_audit.py → mesh_audit.md/.json; cond_audit.py → cond_audit.md
- timecurves.py → timecurves.md/.json; deck_check.py → deck_check.md; misc_checks.txt; fig_map.py
- figures/: map_beta_A_Keff.png, beta_fits.png, beta_A_by_Leff.png, mesh_level.png, tc_loopshape.png
All scripts run on the CPU in < 5 min in total.

Left out of this report: the earlier base-point studies except as references (base point under the main rule: β 1.755 ± 0.053, n 8, same as 7.17); the ring-down and frequency tests; the archived old-axes runs.
