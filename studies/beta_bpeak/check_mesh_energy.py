#!/usr/bin/env python
"""
Code check of the mesh dependence (question of Chris 2026-10-15 after B0b, PROTOKOLL 7.104): does each energy term of
the production setup converge with the cell size for a given smooth magnetisation? The material and geometry are set
up as in run_loops.py (fcc_box, frac, Ms / A / Ku per cell, Herzer cubes, voids with A = 0 and Ms_void). The fields
are given analytically on the cell centres (lengths in l_ex):
  uniform   m = (1, 1, 1) / √3 in all cells;
  helix     m = (cos k z, sin k z, 0), k = 2π / a (one period over the box; no volume charges);
  vortex    one vortex per particle around the axis z through the particle centre, core radius 6 l_ex:
            m = (−y, x, 0) · s / ρ + (0, 0, c) with s = sin θ, c = cos θ, θ = 2 atan(ρ / r_c)... (smooth, resolved
            at all three meshes only for r_c ≥ 2 dx).
Energies per box volume in units of K_d. A correct discretisation gives the same energies on all meshes up to an
error that shrinks with dx (O(dx²) in the volume, O(dx) at the staircase surface).

    python check_mesh_energy.py [d_lex [phi]] > results/mesh_energy.md

Mode "states" (on the instance, needs runs/B0b_*/init.pt): the energy terms of the relaxed B0b states on their own mesh,
and of the relaxed F state (dx 0.75) after the block average onto dx 1.5 and dx 3 (run_loops.transfer_m, before any
relaxation). The same physical state on three meshes shows which term carries the mesh difference.

    python check_mesh_energy.py states > runs/mesh_energy_states.md

Mode "m300" (go of Chris 2026-10-15; 1 µm particles at A 10 pJ/m): the relaxed dx 1.5 states of M300_dx15_s92k_b9
(base microstructure, d/l_ex 300, seeds 921 … 924) on their own mesh and after the block average onto dx 3, and the
relaxed dx 3 states M300_dx3_s92k_b9 if their init.pt exists. Single precision as in production (memory). Reading fixed
before the run: the B0b law (staircase error ≈ 0.22 (dx/d)^1.23 K_d) predicts a demag difference dx 3 − dx 1.5 of
0.00044 K_d for the same state; "the law holds at 1 µm" if the mean over the seeds lies within ±30 % (0.00031 …
0.00057); below: the error at 1 µm is smaller than predicted; above: larger.

    python check_mesh_energy.py m300 > runs/mesh_energy_m300.md

Mode "corrugation" (pinning by the surface): the vortex of each particle is shifted rigidly along x by δ = 0 … 2 dx in
16 steps; the total energy E(δ) minus its quadratic fit (the smooth part of a sphere) gives the corrugation, peak to
peak, for both surfaces and three meshes.

    python check_mesh_energy.py corrugation [d_lex [phi]] > results/mesh_corrugation.md

Mode "solver" (go of Chris 2026-10-15; on the instance, needs runs/B0b_*/init.pt): the production demag DemagFieldPBC
(finite differences for the potential, Bruckner et al. 2021) against the Newell demag DemagField with image shells in a
cube arrangement plus the compensation of the outer surface charges (newell_pbc). Code test before the run: for
uniform m a cubic-symmetric body arrangement has the exact energy φ (1 − φ) / 3 K_d per box volume for every direction
of m; on the B0b geometry the Newell form gives it within 1e-4 with one shell, DemagFieldPBC deviates by +9.5 / 4.7 /
2.4 % on dx 3 / 1.5 / 0.75 for m along (111) and is exact only for m along an axis. States as in mode "states" (seeds
941 … 944; F, F→1.5, F→3, I15, T15, I3, T3), shells 1, 2, 4, double precision. Readings fixed before the run:
  convergence   the Newell energies are converged if |E(2 shells) − E(4 shells)| ≤ 1e-6 K_d for all states;
  solver share  Δ = E_FD − E_Newell(4 shells) per state; s(dx) = [Δ(dx) − Δ(0.75)] / [E_FD(dx) − E_FD(0.75)] for the
                frozen F state is the part of the 7.104 mesh difference that comes from the solver; mean over the
                seeds: s(3) ≥ 0.5 "the solver dominates", s(3) ≤ 0.2 "the solver is minor", else "both are relevant";
  10 % rule     the solver error matters for the 10 % criterion of 7.104 if the mean |Δ| / |E_K| of the relaxed
                dx 3 states (I3, T3) exceeds 0.1 (E_K = anisotropy energy of the state).

    python check_mesh_energy.py solver > runs/demag_solver.md
    python check_mesh_energy.py solver N0b > runs/demag_solver_N0b.md      # states of the reference-solver runs (4.11)

Mode "shells" (decision Chris 2026-10-16, option A of PROTOKOLL 7.108; on the instance, CPU): convergence of the
reference solver in the number k of image shells. In check 2 (results/demag_solver_N0b.md) E_N(4) − E_N(2) reached
13.6 % of E_N(4) on the N0b F states (B0b states: at most 1.5 % of E_N(4)), with the same ratio (E_N(2) − E_N(1)) /
(E_N(4) − E_N(2)) = 2.65 … 2.85 in all 28 states. States: the relaxed F state of <prefix> (seeds 941 … 944) on dx 0.75
and its block averages on dx 1.5 and dx 3; k = 1, 2, 4, 6, 8. Each state also as its zero-mean variant: the
magnetisation M = Ms m minus its box mean <M> (zero_mean). Hypothesis: after the compensation with N_block = 1/3 the
outer surface charges of the image block leave a non-uniform field (N_block(r) − 1/3) <M> in the central box; it
couples the box mean to the non-uniform part of M, falls with the block size as 1 / (2k + 1)² and has either sign.
One process per mesh (mode "shells <prefix> <dx>", JSON runs/shells_<prefix>_dx<dx>.json after each k), then the
report (mode "shells_report <prefix>"). Readings fixed before the run:
  law          q = (E(6) − E(4)) / (E(8) − E(6)) over the states with |E(8) − E(4)| ≥ 1e-6 K_d: the mean is nearer to
               2.62 (law 1 / (2k + 1)²) or to 2.32 (power law k^−1.42, fitted to the ratio 2.67 of check 2);
  cause        "coupling to the box mean" if max |E0(8) − E0(4)| ≤ 0.1 max |E(8) − E(4)| over the states (E0: the
               zero-mean variant), else "not the box mean";
  4 shells     E∞ from the least-squares fit E(k) = c0 + c1 x + c2 x², x = 1 / (2k + 1)², over k = 2, 4, 6, 8;
               r4 = max |E∞ − E(4)| / |E∞| over the states: "4 shells adequate" if r4 ≤ 1 %, "not adequate" if
               r4 ≥ 5 %, else "marginal". Information: r4 against the FD error Δ / E_N(4) of check 2.

    CUDA_DEVICE=-1 python check_mesh_energy.py shells N0b 3 > runs/shells_N0b_dx3.log      # one process per mesh
    python check_mesh_energy.py shells_report N0b > runs/shells_N0b.md

Static check of the extrapolated kernel (option 1' of Chris 2026-10-16, PROTOKOLL 7.108): mode "shells_ext <prefix>
<dx>" adds E(4x), the energy with DemagField(true_pbc=True, extrapolate=True) and 4 shells, to
runs/shells_<prefix>_dx<dx>.json; the report adds the column E(4x). Reading fixed before the run: "extrapolated
kernel valid" if max |E(4x) − E∞| / |E∞| ≤ 0.1 % over the 12 states, else "not valid".

    CUDA_DEVICE=-1 python check_mesh_energy.py shells_ext N0b 1.5 > runs/shells_ext_N0b_dx1.5.log
"""
import math
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent))

from anisotropy_model import box_edge, cube_axes, _cell_centres, FCC_SITES   # noqa: E402
from geometry import fcc_box                                                   # noqa: E402
from units import units                                                        # noqa: E402

LEFF, SEED, JS, A = 12.0, 941, 1.5, 1e-11
PHI = 0.65                                            # the default mode sets it from argv[2] (main)
SHELLS = (1, 2, 4)


def newell_pbc(n, dx, ms, shells, extrapolate=False):
    """Newell demag of magnum.np (DemagField) on an n³ box with `shells` image shells per side (a cube of (2 shells + 1)³
    boxes). The cube of images adds the field −<M>/3 of its outer surface charges at the box; energy_newell_pbc removes
    it, which gives the zero mean field of the true periodic medium (as DemagFieldPBC). ms: Ms per cell [n, n, n, 1].
    In the cube arrangement the image couplings that the dipole branch of DemagField drops at the local origin sum to
    zero for each cell (cubic symmetry). extrapolate: DemagField(true_pbc=True, extrapolate=True), the kernel
    extrapolated in the number of shells (PROTOKOLL 7.108); its field already has the compensation.
    Returns (state, term); the kernel is computed once at the first call."""
    import torch
    from magnumnp import Mesh, State, DemagField
    st = State(Mesh((n, n, n), (dx,) * 3, pbc=(shells,) * 3))
    st.material = {"Ms": ms if torch.is_tensor(ms) else torch.tensor(ms, dtype=torch.get_default_dtype())}
    return st, (DemagField(true_pbc=True, extrapolate=True) if extrapolate else DemagField())


def energy_newell_pbc(st, term, m, Ms, Kd):
    """Demag energy per box volume in K_d of m [n, n, n, 3] with the outer surface charges compensated:
    E = E_images − μ0 V |<M>|² / 6, i.e. E / (K_d V) − |<M> / Ms|² / 3 (a term with true_pbc has it in its field)."""
    import torch
    st.m = m if torch.is_tensor(m) else torch.tensor(m, dtype=torch.get_default_dtype())
    V = float(np.prod(st.mesh.dx)) * float(np.prod(st.mesh.n))
    if getattr(term, "_true_pbc", False):
        return float(term.E(st)) / (Kd * V)
    mean = (st.material["Ms"] * st.m).mean(dim=(0, 1, 2)) / Ms
    return float(term.E(st)) / (Kd * V) - float((mean ** 2).sum()) / 3.0


def setup(d_lex, dx_lex, seed=SEED, surface="voxel"):
    from magnumnp import Mesh, State, DemagFieldPBC, ExchangeField, UniaxialAnisotropyField
    import torch
    U = units(JS, A)
    Ms, Kd, lex = U["Ms"], U["Kd"], U["l_ex"]
    a_lex, d_used = box_edge(d_lex, PHI, 12.0)
    N = int(round(a_lex / dx_lex))
    geo = fcc_box(N, d_used * lex, dx_lex * lex, sub=8 if surface == "fraction" else 0)   # as run_loops
    ids = geo["ids"]
    mag4 = (ids >= 0)[..., None]
    frac4 = geo["frac"].astype(np.float64)[..., None]
    mesh = Mesh((N, N, N), (dx_lex * lex,) * 3, pbc=(1, 1, 1))
    st = State(mesh)
    T = lambda x: torch.tensor(x, dtype=st.dtype, device=st.device)
    st.material = {"Ms": T(np.where(mag4, np.maximum(frac4 * Ms, 1e-8 * Ms), 1e-8 * Ms)),
                   "A": T(np.where(mag4, frac4 * A, 0.0))}
    ax, _ = cube_axes(ids, dx_lex, a_lex, d_used, LEFF, seed)
    st.material["Ku"] = T(np.where(mag4, frac4 * (1.0 / LEFF ** 2) * Kd, 0.0))
    st.material["Ku_axis"] = T(ax)
    terms = {"exchange": ExchangeField(), "anisotropy": UniaxialAnisotropyField(), "demag": DemagFieldPBC()}
    X, Y, Z = _cell_centres(N, dx_lex)
    return st, terms, T, ids, (X, Y, Z), a_lex, d_used, Kd, geo["phi_vox"]


def vortex(ids, XYZ, a, shift=0.0, rc=6.0):
    """One vortex per particle around the z axis through the particle centre, shifted along x by shift (l_ex)."""
    X, Y, Z = XYZ
    m = np.zeros(X.shape + (3,))
    m[..., 2] = 1.0
    for p, s in enumerate(FCC_SITES * a):
        sel = ids == p
        x, y = X[sel] - s[0], Y[sel] - s[1]
        x -= a * np.round(x / a)
        y -= a * np.round(y / a)
        x = x - shift
        rho = np.hypot(x, y)
        th = np.minimum(2 * np.arctan(rho / rc), math.pi / 2)
        phi = np.arctan2(y, x)
        m[sel] = np.stack([-np.sin(phi) * np.sin(th), np.cos(phi) * np.sin(th), np.cos(th)], -1)
    return m


def corrugation(d_lex):
    L = ["# Corrugation of the vortex energy by the surface (check_mesh_energy.py corrugation, d/l_ex %g, φ %g)" % (d_lex, PHI),
         "", "E(δ): total energy per box volume in K_d with the vortex axes shifted by δ = 0 … 2 dx (16 steps). Corrugation "
         "= peak to peak of E(δ) minus its quadratic fit. Reference: anisotropy energy of the same field.", "",
         "| surface | dx (l_ex) | E(0) | E(2 dx) − E(0) | corrugation | corrugation / abs. anisotropy energy |",
         "|---|---|---|---|---|---|"]
    for surface in ("voxel", "fraction"):
        for dx in (3.0, 1.5, 0.75):
            st, terms, T, ids, XYZ, a, d, Kd, phi = setup(d_lex, dx, surface=surface)
            V = float(st.mesh.dx[0] * st.mesh.dx[1] * st.mesh.dx[2]) * st.mesh.n[0] * st.mesh.n[1] * st.mesh.n[2]
            sh = np.linspace(0.0, 2 * dx, 17)
            E, Ean = [], []
            for x in sh:
                st.m = T(vortex(ids, XYZ, a, x))
                e = {t: float(term.E(st)) / (Kd * V) for t, term in terms.items()}
                E.append(sum(e.values()))
                Ean.append(e["anisotropy"])
            E = np.array(E)
            res = E - np.polyval(np.polyfit(sh, E, 2), sh)
            corr = float(res.max() - res.min())
            L.append("| %s | %g | %.6f | %+.2e | %.2e | %.3f |" % (surface, dx, E[0], E[-1] - E[0], corr,
                                                                corr / abs(np.mean(Ean))))
            print(L[-1], file=sys.stderr, flush=True)
    print("\n".join(L))


def fields(ids, XYZ, a, d):
    X, Y, Z = XYZ
    k = 2 * math.pi / a
    out = {"uniform": np.broadcast_to(np.array([1.0, 1.0, 1.0]) / math.sqrt(3), X.shape + (3,)).copy(),
           "helix": np.stack([np.cos(k * Z), np.sin(k * Z), np.zeros_like(Z)], -1)}
    m = np.zeros(X.shape + (3,))
    m[..., 2] = 1.0
    rc = 6.0
    for p, s in enumerate(FCC_SITES * a):
        sel = ids == p
        x, y = X[sel] - s[0], Y[sel] - s[1]
        x -= a * np.round(x / a)
        y -= a * np.round(y / a)
        rho = np.hypot(x, y)
        th = 2 * np.arctan(rho / rc)                 # 0 at the axis, → π at large rho: a smooth core
        th = np.minimum(th, math.pi / 2)             # in-plane outside the core
        phi = np.arctan2(y, x)
        m[sel] = np.stack([-np.sin(phi) * np.sin(th), np.cos(phi) * np.sin(th), np.cos(th)], -1)
    out["vortex"] = m
    return out


def states():
    """Energy terms of the relaxed B0b states and of the transferred F state (d/l_ex 96, φ 0.65)."""
    import torch
    import run_loops
    L = ["# Energy terms of the relaxed B0b states (check_mesh_energy.py states)", "",
         "Per box volume in units of K_d. 'F→dx' = the relaxed F state, block-averaged onto the coarse mesh, not relaxed.", "",
         "| seed | state | dx | exchange | anisotropy | demag | total | init record E / K_d |", "|---|---|---|---|---|---|---|---|"]
    for seed in (941, 942, 943, 944):
        S = {}
        for dx in (0.75, 1.5, 3.0):
            S[dx] = setup(96.0, dx, seed)
        mF = torch.load(HERE / "runs" / ("B0b_F_s%d" % seed) / "init.pt", map_location="cpu", weights_only=False)
        cases = [("F", 0.75, mF["m"], mF["init"]["E_over_Kd"])]
        for dx in (1.5, 3.0):
            mc, _ = run_loops.transfer_m(mF["m"].double().numpy(), S[0.75][3], S[dx][3])
            cases.append(("F→dx", dx, mc, float("nan")))
        for arm, dx in (("I15", 1.5), ("T15", 1.5), ("I3", 3.0), ("T3", 3.0)):
            c = torch.load(HERE / "runs" / ("B0b_%s_s%d" % (arm, seed)) / "init.pt", map_location="cpu", weights_only=False)
            cases.append((arm, dx, c["m"], c["init"]["E_over_Kd"]))
        for name, dx, m, e0 in cases:
            st, terms, T, ids, XYZ, a, d, Kd, phi = S[dx]
            V = float(np.prod(st.mesh.dx)) * float(np.prod(st.mesh.n))
            m = m.detach().cpu().double().numpy() if torch.is_tensor(m) else np.asarray(m, dtype=np.float64)
            st.m = T(m)
            e = {t: float(term.E(st)) / (Kd * V) for t, term in terms.items()}
            L.append("| %d | %s | %g | %.6f | %.6f | %.6f | %.6f | %.6f |" % (seed, name, dx, e["exchange"], e["anisotropy"],
                                                                         e["demag"], sum(e.values()), e0))
            print(L[-1], file=sys.stderr, flush=True)
    print("\n".join(L))


def solver(prefix="B0b"):
    """DemagFieldPBC against Newell with compensated outer surface charges on the states of the runs <prefix>_* (B0b:
    readings in the module doc; N0b, the runs with the reference solver, PLAN 4.11: information only)."""
    import json
    import torch
    import run_loops
    U = units(JS, A)
    Ms = U["Ms"]
    S, NW = {}, {}
    for dx in (0.75, 1.5, 3.0):
        S[dx] = setup(96.0, dx, SEED)
        st = S[dx][0]
        for sh in SHELLS:
            NW[(dx, sh)] = newell_pbc(st.mesh.n[0], float(st.mesh.dx[0]), st.material["Ms"], sh)
    rows = []
    for seed in (941, 942, 943, 944):
        for dx in (0.75, 1.5, 3.0):                   # the cube axes of this seed (anisotropy energy E_K)
            st, terms, T, ids, XYZ, a, d, Kd, phi = S[dx]
            ax, _ = cube_axes(ids, dx, a, d, LEFF, seed)
            st.material["Ku_axis"] = T(ax)
        mF = torch.load(HERE / "runs" / ("%s_F_s%d" % (prefix, seed)) / "init.pt", map_location="cpu", weights_only=False)
        cases = [("F", 0.75, mF["m"])]
        for dx in (1.5, 3.0):
            mc, _ = run_loops.transfer_m(mF["m"].double().cpu().numpy(), S[0.75][3], S[dx][3])
            cases.append(("F→%g" % dx, dx, mc))
        for arm, dx in (("I15", 1.5), ("T15", 1.5), ("I3", 3.0), ("T3", 3.0)):
            c = torch.load(HERE / "runs" / ("%s_%s_s%d" % (prefix, arm, seed)) / "init.pt", map_location="cpu", weights_only=False)
            cases.append((arm, dx, c["m"]))
        for name, dx, m in cases:
            st, terms, T, ids, XYZ, a, d, Kd, phi = S[dx]
            V = float(np.prod(st.mesh.dx)) * float(np.prod(st.mesh.n))
            m = T(m.detach().cpu().double().numpy() if torch.is_tensor(m) else np.asarray(m, dtype=np.float64))
            st.m = m
            r = {"seed": seed, "state": name, "dx": dx, "E_K": float(terms["anisotropy"].E(st)) / (Kd * V),
                 "E_FD": float(terms["demag"].E(st)) / (Kd * V)}
            for sh in SHELLS:
                r["E_N%d" % sh] = energy_newell_pbc(*NW[(dx, sh)], m, Ms, Kd)
            r["delta"] = r["E_FD"] - r["E_N%d" % SHELLS[-1]]
            rows.append(r)
            print("%d %s %g %.7f %.7f %+.2e" % (seed, name, dx, r["E_FD"], r["E_N%d" % SHELLS[-1]], r["delta"]),
                  file=sys.stderr, flush=True)
    sh = SHELLS[-1]
    L = ["# Demag solver check on the %s states (check_mesh_energy.py solver %s)%s" % (
        prefix, prefix, "" if prefix == "B0b" else "; the readings below were fixed for B0b, here information only"), "",
         "Per box volume in units of K_d, double precision. E_FD = DemagFieldPBC (production); E_N(k) = Newell with k "
         "image shells and compensated outer surface charges; Δ = E_FD − E_N(%d); E_K = anisotropy energy. "
         "'F→dx' = the relaxed F state, block-averaged onto dx, not relaxed." % sh, "",
         "| seed | state | dx | E_K | E_FD | " + " | ".join("E_N(%d)" % k for k in SHELLS) + " | Δ | Δ / abs(E_K) |",
         "|---|---|---|---|---|" + "---|" * len(SHELLS) + "---|---|"]
    for r in rows:
        L.append("| %d | %s | %g | %.6f | %.6f | " % (r["seed"], r["state"], r["dx"], r["E_K"], r["E_FD"])
                 + " | ".join("%.7f" % r["E_N%d" % k] for k in SHELLS)
                 + " | %+.6f | %.2f |" % (r["delta"], abs(r["delta"] / r["E_K"])))
    conv = max(abs(r["E_N%d" % SHELLS[-2]] - r["E_N%d" % sh]) for r in rows)
    g = lambda seed, state: next(r for r in rows if r["seed"] == seed and r["state"] == state)
    seeds = (941, 942, 943, 944)
    out = {"convergence_max_abs": conv, "rows": rows}
    L += ["", "Convergence: max |E_N(%d) − E_N(%d)| = %.1e K_d over all states (limit 1e-6): **%s**."
          % (SHELLS[-2], sh, conv, "converged" if conv <= 1e-6 else "not converged")]
    for dx in (3.0, 1.5):
        dFD = np.array([g(s, "F→%g" % dx)["E_FD"] - g(s, "F")["E_FD"] for s in seeds])
        dN = np.array([g(s, "F→%g" % dx)["E_N%d" % sh] - g(s, "F")["E_N%d" % sh] for s in seeds])
        sh_ = (dFD - dN) / dFD
        out["frozen_dx%g" % dx] = {"dE_FD": dFD.tolist(), "dE_N": dN.tolist(), "share": sh_.tolist()}
        line = ("Frozen F state, dx %g − dx 0.75: E_FD difference %.5f ± %.5f, Newell difference %.5f ± %.5f K_d "
                "(mean ± SE, n 4); solver share s = %.2f ± %.2f" % (dx, dFD.mean(), dFD.std(ddof=1) / 2, dN.mean(),
                                                                 dN.std(ddof=1) / 2, sh_.mean(), sh_.std(ddof=1) / 2))
        if dx == 3.0:
            s3 = sh_.mean()
            rd = "the solver dominates" if s3 >= 0.5 else ("the solver is minor" if s3 <= 0.2 else "both are relevant")
            out["reading_share"] = rd
            line += ". **Reading: %s.**" % rd
        L.append(line + ("." if dx != 3.0 else ""))
    q = np.array([abs(r["delta"] / r["E_K"]) for r in rows if r["state"] in ("I3", "T3")])
    rd = "relevant" if q.mean() > 0.1 else "not relevant"
    out["reading_10pc"] = rd
    out["delta_over_EK_dx3_relaxed"] = q.tolist()
    L.append("Relaxed dx 3 states (I3, T3): |Δ| / |E_K| = %.2f ± %.2f (mean ± SE, n %d). **Reading: the solver error "
             "is %s for the 10 %% criterion.**" % (q.mean(), q.std(ddof=1) / math.sqrt(len(q)), len(q), rd))
    (HERE / "runs" / ("demag_solver.json" if prefix == "B0b" else "demag_solver_%s.json" % prefix)).write_text(
        json.dumps(out, indent=1))
    print("\n".join(L))


def m300():
    """Energy terms of the M300 states (d/l_ex 300) on dx 1.5 and dx 3; one mesh and demag kernel per dx, the cube axes
    per seed."""
    import torch
    import run_loops
    from anisotropy_model import cube_axes
    L = ["# Energy terms of the M300 states, d/l_ex 300 (check_mesh_energy.py m300)", "",
         "Per box volume in units of K_d, single precision. 'dx15→3' = the relaxed dx 1.5 state, block-averaged onto dx 3, "
         "not relaxed.", "", "| seed | state | dx | exchange | anisotropy | demag | total | init record E / K_d |",
         "|---|---|---|---|---|---|---|---|"]
    S = {dx: setup(300.0, dx, 921) for dx in (1.5, 3.0)}
    out = {}
    for seed in (921, 922, 923, 924):
        for dx in (1.5, 3.0):
            st, terms, T, ids, XYZ, a, d, Kd, phi = S[dx]
            ax, _ = cube_axes(ids, dx, a, d, LEFF, seed)
            st.material["Ku_axis"] = T(ax)
            del ax
        R = HERE / "runs"
        c15 = torch.load(R / ("M300_dx15_s%d_b9" % seed) / "init.pt", map_location="cpu", weights_only=False)
        m15 = c15["m"].detach().cpu().double().numpy()
        mc, _ = run_loops.transfer_m(m15, S[1.5][3], S[3.0][3])
        cases = [("dx15", 1.5, m15, c15["init"]["E_over_Kd"]), ("dx15→3", 3.0, mc, float("nan"))]
        f3 = R / ("M300_dx3_s%d_b9" % seed) / "init.pt"
        if f3.exists():
            c3 = torch.load(f3, map_location="cpu", weights_only=False)
            cases.append(("dx3", 3.0, c3["m"], c3["init"]["E_over_Kd"]))
        for name, dx, m, e0 in cases:
            st, terms, T, ids, XYZ, a, d, Kd, phi = S[dx]
            V = float(st.mesh.dx[0] * st.mesh.dx[1] * st.mesh.dx[2]) * st.mesh.n[0] * st.mesh.n[1] * st.mesh.n[2]
            m = m.detach().cpu().double().numpy() if torch.is_tensor(m) else np.asarray(m, dtype=np.float64)
            st.m = T(m)
            e = {t: float(term.E(st)) / (Kd * V) for t, term in terms.items()}
            out[(seed, name)] = e
            L.append("| %d | %s | %g | %.6f | %.6f | %.6f | %.6f | %.6f |" % (seed, name, dx, e["exchange"], e["anisotropy"],
                                                                         e["demag"], sum(e.values()), e0))
            print(L[-1], file=sys.stderr, flush=True)
    dd = np.array([out[(s, "dx15→3")]["demag"] - out[(s, "dx15")]["demag"] for s in (921, 922, 923, 924)])
    an = np.array([abs(out[(s, "dx15")]["anisotropy"]) for s in (921, 922, 923, 924)])
    m, se = dd.mean(), dd.std(ddof=1) / 2.0
    rd = "the law holds at 1 µm" if 0.00031 <= m <= 0.00057 else ("smaller than predicted" if m < 0.00031 else "larger than predicted")
    e15 = dd / (2 ** 1.23 - 1)                        # staircase error at dx 1.5 with the order of B0b (estimate)
    L += ["", "Demag difference dx 3 − dx 1.5 for the same state: %s K_d per seed, mean %.5f ± %.5f (SE, n 4). Prediction "
          "0.00044. **Reading: %s.**" % (" / ".join("%.5f" % x for x in dd), m, se, rd),
          "Estimate with the order 1.23 of B0b: staircase error / |anisotropy| = %.2f at dx 1.5 and %.2f at dx 3 (mean)."
          % (np.mean(e15 / an), np.mean(e15 * 2 ** 1.23 / an))]
    print("\n".join(L))


SHELLS_EXT = (1, 2, 4, 6, 8)
Q_POWER, Q_INV2 = 2.32, 2.62          # (E(6) − E(4)) / (E(8) − E(6)) for k^−1.42 and for 1 / (2k + 1)²
#                                       (2.32: corrected 2026-10-16 07:40 UTC from 2.30, a computing error; the
#                                       boundary between the two laws moves from 2.46 to 2.47)


def zero_mean(m, ms):
    """(m0, ms0): the magnetisation Ms m without its box mean, as a field m0 with the constant Ms0 = max Ms in every
    cell (the void cells carry the negative mean): ms0 · m0 = ms · m − <ms · m>. m [n, n, n, 3], ms [n, n, n, 1]."""
    import torch
    M = ms * m
    M0 = M - M.mean(dim=(0, 1, 2), keepdim=True)
    s = float(ms.max())
    return M0 / s, torch.full_like(ms, s)


def shells(prefix, dx, ks=SHELLS_EXT):
    """E_N(k) of the relaxed F states of <prefix> on mesh dx (0.75: the state itself; 1.5, 3: its block average) and of
    their zero-mean variants, for the shells ks; the JSON is written after each k."""
    import json
    import torch
    import run_loops
    U = units(JS, A)
    Ms = U["Ms"]
    S = {d: setup(96.0, d, SEED) for d in sorted({0.75, dx})}
    st, terms, T, ids, XYZ, a, d, Kd, phi = S[dx]
    ms = st.material["Ms"]
    cases = []
    for seed in (941, 942, 943, 944):
        mF = torch.load(HERE / "runs" / ("%s_F_s%d" % (prefix, seed)) / "init.pt", map_location="cpu", weights_only=False)
        if dx == 0.75:
            m = mF["m"].double().cpu().numpy()
        else:
            m, _ = run_loops.transfer_m(mF["m"].double().cpu().numpy(), S[0.75][3], S[dx][3])
        cases.append((seed, T(np.asarray(m, dtype=np.float64))))
    f = HERE / "runs" / ("shells_%s_dx%g.json" % (prefix, dx))
    out = json.loads(f.read_text()) if f.exists() else {}
    for k in ks:
        if str(k) in out:
            continue
        kk, ext = int(str(k).rstrip("x")), str(k).endswith("x")       # "4x": 4 shells, extrapolated kernel
        stk, term = newell_pbc(st.mesh.n[0], float(st.mesh.dx[0]), ms, kk, extrapolate=ext)
        res = {}
        for seed, m in cases:
            e = energy_newell_pbc(stk, term, m, Ms, Kd)
            m0, ms0 = zero_mean(m, ms)
            stk.material["Ms"] = ms0
            e0 = energy_newell_pbc(stk, term, m0, Ms, Kd)
            stk.material["Ms"] = ms
            res[str(seed)] = {"E": e, "E0": e0}
            print("dx %g k %s seed %d E %.8f E0 %.8f" % (dx, k, seed, e, e0), file=sys.stderr, flush=True)
        out[str(k)] = res
        f.write_text(json.dumps(out, indent=1))
        del term, stk


def shell_fit(ks, E):
    """E∞ of the least-squares fit E(k) = c0 + c1 x + c2 x², x = 1 / (2k + 1)², over ks (at least 3 values)."""
    x = 1.0 / (2.0 * np.asarray(ks, dtype=float) + 1.0) ** 2
    return float(np.linalg.lstsq(np.stack([np.ones_like(x), x, x * x], 1), np.asarray(E, dtype=float), rcond=None)[0][0])


def shell_readings(rows):
    """Readings of mode "shells" (module doc). rows: [{"E": {k: E}, "E0": {k: E}}] with k = 1, 2, 4, 6, 8."""
    big = [r for r in rows if abs(r["E"][8] - r["E"][4]) >= 1e-6]
    q = [(r["E"][6] - r["E"][4]) / (r["E"][8] - r["E"][6]) for r in big]
    qm = float(np.mean(q)) if q else float("nan")
    law = ("not determined" if not q else "1 / (2k + 1)²" if abs(qm - Q_INV2) < abs(qm - Q_POWER) else "power law k^−1.42")
    d = max(abs(r["E"][8] - r["E"][4]) for r in rows)
    d0 = max(abs(r["E0"][8] - r["E0"][4]) for r in rows)
    cause = "coupling to the box mean" if d0 <= 0.1 * d else "not the box mean"
    r4 = max(abs(shell_fit((2, 4, 6, 8), [r["E"][k] for k in (2, 4, 6, 8)]) - r["E"][4]) /
             abs(shell_fit((2, 4, 6, 8), [r["E"][k] for k in (2, 4, 6, 8)])) for r in rows)
    ad = "4 shells adequate" if r4 <= 0.01 else ("not adequate" if r4 >= 0.05 else "marginal")
    return {"q_mean": qm, "q": q, "law": law, "d": d, "d0": d0, "cause": cause, "r4": r4, "adequacy": ad}


def shells_report(prefix):
    import json
    ext = []
    rows, L = [], ["# Shell convergence of the reference solver on the %s F states (check_mesh_energy.py shells, "
                   "PROTOKOLL 7.108)" % prefix, "",
                   "Per box volume in units of K_d, double precision. E(k) = Newell with k image shells and compensated "
                   "outer surface charges; E0(k) = the same for the zero-mean variant (box mean of M removed); E∞ = fit "
                   "c0 + c1 x + c2 x², x = 1 / (2k + 1)², over k = 2, 4, 6, 8.", "",
                   "| seed | state | " + " | ".join("E(%d)" % k for k in SHELLS_EXT) + " | E∞ | (E∞ − E(4)) / E∞ | "
                   + " | ".join("E0(%d)" % k for k in SHELLS_EXT) + " |",
                   "|---|---|" + "---|" * (2 * len(SHELLS_EXT) + 2)]
    for dx, name in ((0.75, "F"), (1.5, "F→1.5"), (3.0, "F→3")):
        f = HERE / "runs" / ("shells_%s_dx%g.json" % (prefix, dx))
        if not f.exists():
            L.append("| – | %s | not computed |" % name)
            continue
        out = json.loads(f.read_text())
        if not all(str(k) in out for k in SHELLS_EXT):
            L.append("| – | %s | only k = %s |" % (name, ", ".join(sorted(out, key=int))))
            continue
        for seed in (941, 942, 943, 944):
            r = {"E": {k: out[str(k)][str(seed)]["E"] for k in SHELLS_EXT},
                 "E0": {k: out[str(k)][str(seed)]["E0"] for k in SHELLS_EXT}}
            einf = shell_fit((2, 4, 6, 8), [r["E"][k] for k in (2, 4, 6, 8)])
            if "4x" in out:
                ext.append((seed, name, out["4x"][str(seed)]["E"], einf))
            rows.append(r)
            L.append("| %d | %s | " % (seed, name) + " | ".join("%.7f" % r["E"][k] for k in SHELLS_EXT)
                     + " | %.7f | %+.2f %% | " % (einf, 100 * (einf - r["E"][4]) / einf)
                     + " | ".join("%.7f" % r["E0"][k] for k in SHELLS_EXT) + " |")
    if rows:
        R = shell_readings(rows)
        L += ["", "Law: q = (E(6) − E(4)) / (E(8) − E(6)) = %.2f (mean over %d states; 1 / (2k + 1)² gives %.2f, k^−1.42 "
              "gives %.2f). **Reading: %s.**" % (R["q_mean"], len(R["q"]), Q_INV2, Q_POWER, R["law"]),
              "Cause: max |E(8) − E(4)| = %.1e K_d, zero-mean variant max |E0(8) − E0(4)| = %.1e K_d. **Reading: %s.**"
              % (R["d"], R["d0"], R["cause"]),
              "4 shells: r4 = max |E∞ − E(4)| / |E∞| = %.2f %% over %d states. **Reading: %s.**"
              % (100 * R["r4"], len(rows), R["adequacy"])]
        if ext:
            L += ["", "Extrapolated kernel (DemagField(extrapolate=True), 4 shells): " + "; ".join(
                "%d %s E(4x) %.7f, (E(4x) − E∞) / E∞ %+.3f %%" % (s, nm, e, 100 * (e - ei) / ei) for s, nm, e, ei in ext)]
            mx = max(abs(e - ei) / abs(ei) for _, _, e, ei in ext)
            R["ext_max_rel"] = mx
            R["ext_reading"] = ("extrapolated kernel valid" if mx <= 0.001 and len(ext) == len(rows)
                                else "not valid" if mx > 0.001 else "not determined (not all states)")
            L.append("max |E(4x) − E∞| / |E∞| = %.4f %% over %d states. **Reading: %s.**" % (100 * mx, len(ext),
                                                                                            R["ext_reading"]))
        (HERE / "runs" / ("shells_%s.json" % prefix)).write_text(json.dumps(R, indent=1))
    print("\n".join(L))


def main():
    from magnumnp import set_precision
    if len(sys.argv) > 1 and sys.argv[1] == "m300":
        set_precision("single")                      # memory: 26 M cells at dx 1.5, as in production
        return m300()
    set_precision("double")                          # once, before the first mesh
    if len(sys.argv) > 1 and sys.argv[1] == "states":
        return states()
    if len(sys.argv) > 1 and sys.argv[1] == "solver":
        return solver(sys.argv[2] if len(sys.argv) > 2 else "B0b")
    if len(sys.argv) > 3 and sys.argv[1] == "shells":
        return shells(sys.argv[2], float(sys.argv[3]))
    if len(sys.argv) > 3 and sys.argv[1] == "shells_ext":
        return shells(sys.argv[2], float(sys.argv[3]), ks=("4x",))
    if len(sys.argv) > 2 and sys.argv[1] == "shells_report":
        return shells_report(sys.argv[2])
    global PHI
    if len(sys.argv) > 1 and sys.argv[1] == "corrugation":
        PHI = float(sys.argv[3]) if len(sys.argv) > 3 else PHI
        return corrugation(float(sys.argv[2]) if len(sys.argv) > 2 else 48.0)
    PHI = float(sys.argv[2]) if len(sys.argv) > 2 else PHI
    d_lex = float(sys.argv[1]) if len(sys.argv) > 1 else 48.0
    rows = {}
    for surface in ("voxel", "fraction"):
        for dx in (3.0, 1.5, 0.75):
            st, terms, T, ids, XYZ, a, d, Kd, phi = setup(d_lex, dx, surface=surface)
            V = float(st.mesh.dx[0] * st.mesh.dx[1] * st.mesh.dx[2]) * st.mesh.n[0] * st.mesh.n[1] * st.mesh.n[2]
            for name, m in fields(ids, XYZ, a, d).items():
                st.m = T(m)
                rows[(surface, name, dx)] = {t: float(term.E(st)) / (Kd * V) for t, term in terms.items()}
                rows[(surface, name, dx)]["phi_vox"] = phi
    L = ["# Energy terms of given smooth fields on three meshes (check_mesh_energy.py, d/l_ex %g, φ %g)" % (d_lex, PHI), "",
         "Energy per box volume in units of K_d; production material and geometry; double precision. Surface: voxel = "
         "staircase (production); fraction = surface cells with Ms, A, K scaled by their volume fraction (--surface "
         "fraction).", "",
         "| surface | field | dx (l_ex) | φ_vox | exchange | anisotropy | demag |", "|---|---|---|---|---|---|---|"]
    for surface in ("voxel", "fraction"):
        for name in ("uniform", "helix", "vortex"):
            for dx in (3.0, 1.5, 0.75):
                r = rows[(surface, name, dx)]
                L.append("| %s | %s | %g | %.4f | %.6f | %.6f | %.6f |" % (surface, name, dx, r["phi_vox"], r["exchange"],
                                                                          r["anisotropy"], r["demag"]))
    k = 2 * math.pi / box_edge(d_lex, PHI, 12.0)[0]
    L += ["", "Analytic exchange of the helix (continuum): φ (l_ex k)² = %.6f at φ %g." % (PHI * k * k, PHI)]
    print("\n".join(L))


if __name__ == "__main__":
    main()
