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
PHI = float(sys.argv[2]) if len(sys.argv) > 2 else 0.65


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


def main():
    from magnumnp import set_precision
    if len(sys.argv) > 1 and sys.argv[1] == "m300":
        set_precision("single")                      # memory: 26 M cells at dx 1.5, as in production
        return m300()
    set_precision("double")                          # once, before the first mesh
    if len(sys.argv) > 1 and sys.argv[1] == "states":
        return states()
    if len(sys.argv) > 1 and sys.argv[1] == "corrugation":
        global PHI
        PHI = float(sys.argv[3]) if len(sys.argv) > 3 else PHI
        return corrugation(float(sys.argv[2]) if len(sys.argv) > 2 else 48.0)
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
