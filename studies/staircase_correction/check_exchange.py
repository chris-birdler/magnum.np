"""
Static exchange test (PROTOKOLL 3.13): error of the exchange energy at the particle surface, voxel model against the
fraction surface of V1 (A v, M_s v; the exchange of magnum.np with the harmonic mean of A between cells).

Helix m = (cos(k·r), sin(k·r), 0) with a fixed wave vector k: in the continuum the exchange energy on the smooth spheres
is A |k|² V_body. The finite-difference stencil gives in an infinite bulk the density
    A Σ_i 2 (1 − cos(k_i Δx)) / Δx²,
so the reference E_ref = that density times V_body isolates the surface error (the stencil error of the bulk drops
out). Control: a box without surfaces (all cells magnetic) must give E_ref exactly. Vortex field (information): a
vortex around an axis through each particle centre, against the continuum integral (stencil and surface error).

    python check_exchange.py --out results/exchange_static.json
"""
import argparse
import json
import math
import os
import pathlib
import sys

import numpy as np
import torch

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "beta_bpeak"))
os.environ.setdefault("CUDA_DEVICE", "-1")
from magnumnp import set_precision, Mesh, State, ExchangeField  # noqa: E402
set_precision("double")
torch.set_num_threads(int(os.environ.get("STAIRCASE_THREADS", os.cpu_count())))
from anisotropy_model import box_edge  # noqa: E402
from geometry import fcc_box, FCC_SITES  # noqa: E402

A0, MS0 = 1e-11, 1.0e6
GEOMS = {"K2": (300.0, 0.62), "K4": (300.0, 0.55), "d96": (96.0, 0.55)}


def exchange_energy(N, dx, A, Ms, m):
    st = State(Mesh((N, N, N), (dx, dx, dx), pbc=(1, 1, 1)))
    st.material["A"] = torch.as_tensor(A).reshape(N, N, N, 1)
    st.material["Ms"] = torch.as_tensor(np.maximum(Ms, 1e-30)).reshape(N, N, N, 1)
    st.m = torch.as_tensor(m)
    with torch.no_grad():
        return float(ExchangeField().E(st))


def helix(N, dx, kvec):
    x = (np.arange(N) + 0.5) * dx
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    ph = kvec[0] * X + kvec[1] * Y + kvec[2] * Z
    return np.stack([np.cos(ph), np.sin(ph), np.zeros_like(ph)], -1)


def stencil_density(kvec, dx):
    return A0 * sum(2.0 * (1.0 - math.cos(k * dx)) / dx ** 2 for k in kvec)


def vortex(N, dx, a, R, rc):
    """Vortex around the z axis through each particle centre (minimum image), core profile θ(ρ) = (π/2) min(1, ρ/rc)."""
    x = (np.arange(N) + 0.5) * dx
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    m = np.zeros((N, N, N, 3))
    best = np.full((N, N, N), np.inf)
    for s in FCC_SITES * a:
        dX = X - s[0]; dX -= a * np.round(dX / a)
        dY = Y - s[1]; dY -= a * np.round(dY / a)
        dZ = Z - s[2]; dZ -= a * np.round(dZ / a)
        r2 = dX ** 2 + dY ** 2 + dZ ** 2
        sel = r2 < best
        rho = np.sqrt(dX ** 2 + dY ** 2) + 1e-300
        th = 0.5 * np.pi * np.minimum(1.0, rho / rc)
        mm = np.stack([-dY / rho * np.sin(th), dX / rho * np.sin(th), np.cos(th)], -1)
        m[sel] = mm[sel]
        best = np.where(sel, r2, best)
    return m


def vortex_continuum(R, rc, n=4000):
    """A ∫ over one sphere of |∇m|² = θ'(ρ)² + sin²θ / ρ², by quadrature in ρ (exact in z)."""
    rho = (np.arange(n) + 0.5) / n * R
    th = 0.5 * np.pi * np.minimum(1.0, rho / rc)
    dth = np.where(rho < rc, 0.5 * np.pi / rc, 0.0)
    g = dth ** 2 + np.sin(th) ** 2 / rho ** 2
    height = 2.0 * np.sqrt(np.maximum(R ** 2 - rho ** 2, 0.0))
    return A0 * float(np.sum(g * height * 2.0 * np.pi * rho) * (R / n))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE / "results" / "exchange_static.json"))
    ap.add_argument("--dx_lex", type=float, default=3.0)
    a = ap.parse_args()
    res = {"control": {}, "helix": {}, "vortex": {}}
    lex = 1.0                                   # lengths in l_ex (the energy ratios do not depend on the unit)
    # control: box without surfaces
    N = 32
    dx = a.dx_lex
    for name, kd in (("x", (1, 0, 0)), ("111", (1, 1, 1))):
        kv = np.array(kd, float) / np.linalg.norm(kd) * 2 * np.pi / (N * dx / 2)   # two periods? one per half box
        kv = np.round(kv * N * dx / (2 * np.pi)) * 2 * np.pi / (N * dx)          # commensurate with the box
        E = exchange_energy(N, dx, np.full((N, N, N), A0), np.full((N, N, N), MS0), helix(N, dx, kv))
        Eref = stencil_density(kv, dx) * (N * dx) ** 3
        res["control"][name] = E / Eref - 1
    print("control (box without surfaces):", res["control"], flush=True)
    for g, (d_nom, phi) in GEOMS.items():
        a_lex, d_lex = box_edge(d_nom, phi, 12.0)
        N = int(round(a_lex / dx))
        R = d_lex / 2
        V_body = 4 * (4.0 / 3.0) * math.pi * R ** 3
        models = {}
        gv = fcc_box(N, d_lex, dx, sub=0)
        models["voxel"] = (gv["frac"].astype(float), gv["frac"].astype(float))
        for q in (2, 4):
            gq = fcc_box(N, d_lex, dx, sub=q)
            models["V1_q%d" % q] = (gq["frac"].astype(float), gq["frac"].astype(float))
        res["helix"][g], res["vortex"][g] = {}, {}
        for kname, kd in (("100", (1, 0, 0)), ("110", (1, 1, 0)), ("111", (1, 1, 1))):
            for lam in (d_lex / 2, d_lex / 4):
                kv = np.array(kd, float) / np.linalg.norm(kd) * 2 * np.pi / lam
                m = helix(N, dx, kv)
                Eref = stencil_density(kv, dx) * V_body
                Econt = A0 * float(kv @ kv) * V_body
                row = {"lambda_lex": lam, "stencil_vs_continuum": Eref / Econt - 1}
                for mn, (f, fm) in models.items():
                    row[mn] = exchange_energy(N, dx, A0 * f, MS0 * fm, m) / Eref - 1
                res["helix"][g]["%s_lam%.0f" % (kname, lam)] = row
                print(g, kname, "lam %.1f" % lam, {k: round(v, 4) for k, v in row.items()}, flush=True)
        for rcf in (0.1, 0.3):
            rc = rcf * R
            m = vortex(N, dx, a_lex, R, rc)
            Ec = 4 * vortex_continuum(R, rc)
            row = {"rc_over_R": rcf}
            for mn, (f, fm) in models.items():
                row[mn] = exchange_energy(N, dx, A0 * f, MS0 * fm, m) / Ec - 1
            res["vortex"][g]["rc%.1f" % rcf] = row
            print(g, "vortex", row, flush=True)
    pathlib.Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    pathlib.Path(a.out).write_text(json.dumps(res, indent=1))
    print("done", a.out)


if __name__ == "__main__":
    main()
