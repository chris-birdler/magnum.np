#!/usr/bin/env python
"""
Demag solver tests with uniform m and with the upstream code (PROTOKOLL 7.105; local CPU, a few minutes):

  python check_demag_uniform.py > results/demag_uniform.md

(1) B0b geometry (d/l_ex 97.5, nominal 96, φ 0.65, voxel spheres on cell vertices, FCC sites): uniform m along (100),
    (110), (111). Exact reference: φ_vox (1 − φ_vox) / 3 K_d per box volume for every direction (cubic symmetry, zero
    mean field). DemagFieldPBC on dx 3 / 1.5 / 0.75; Newell with one image shell and compensated outer surface charges
    (check_mesh_energy.newell_pbc) on dx 3 / 1.5.
(2) Aligned cube of edge a/2 (volume fraction 1/8, exact 7/192) with 4 / 8 / 16 / 32 cells per edge, m along (100),
    (110), (111): DemagFieldPBC; Newell (1 shell) for 4 and 8 cells.
(3) DemagFieldPBC of our code against the upstream code (git fde0c62, the last upstream change, before our commit
    8654eb7): largest relative difference of h for a random m on the B0b geometry at dx 3.
"""
import importlib.util
import math
import pathlib
import subprocess
import sys

import numpy as np
import torch

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent))
from magnumnp import Mesh, State, DemagFieldPBC, set_precision      # noqa: E402
set_precision("double")
import check_mesh_energy as C                                        # noqa: E402
from anisotropy_model import box_edge                                # noqa: E402
from geometry import fcc_box                                         # noqa: E402
from units import units                                              # noqa: E402

U = units(1.5, 1e-11)
Ms, Kd, lex = U["Ms"], U["Kd"], U["l_ex"]
DIRS = {"(100)": (1.0, 0.0, 0.0), "(110)": (1.0, 1.0, 0.0), "(111)": (1.0, 1.0, 1.0)}


def uni(n, u):
    u = np.asarray(u) / np.linalg.norm(u)
    return torch.tensor(np.broadcast_to(u, (n, n, n, 3)).copy())


def fd_energy(ms, m, dx):
    n = ms.shape[0]
    st = State(Mesh((n, n, n), (dx,) * 3, pbc=(1, 1, 1)))
    st.material = {"Ms": torch.tensor(ms[..., None])}
    st.m = m
    return float(DemagFieldPBC().E(st)) / (Kd * dx ** 3 * n ** 3)


def main():
    L = ["# Demag solver tests with uniform m and against the upstream code (check_demag_uniform.py)", "",
         "Energy per box volume in units of K_d, double precision, local CPU. Relative error = E / E_exact − 1. "
         "Newell: one image shell, compensated outer surface charges (check_mesh_energy.newell_pbc).", "",
         "## 1. B0b geometry, exact φ_vox (1 − φ_vox) / 3", "",
         "| dx (l_ex) | φ_vox | exact | m along | DemagFieldPBC | rel. error | Newell | rel. error |", "|---|---|---|---|---|---|---|---|"]
    a_lex, d_used = box_edge(96.0, 0.65, 12.0)
    for dxl in (3.0, 1.5, 0.75):
        N = int(round(a_lex / dxl))
        dx = dxl * lex
        chi = fcc_box(N, d_used * lex, dx)["ids"] >= 0
        ms = np.where(chi, Ms, 1e-8 * Ms)
        f = chi.mean()
        ex = f * (1 - f) / 3
        nw = C.newell_pbc(N, dx, torch.tensor(ms[..., None]), 1) if dxl >= 1.5 else None
        for name, u in DIRS.items():
            e = fd_energy(ms, uni(N, u), dx)
            en = C.energy_newell_pbc(*nw, uni(N, u), Ms, Kd) if nw else float("nan")
            L.append("| %g | %.5f | %.6f | %s | %.6f | %+.2e | %s | %s |" % (
                dxl, f, ex, name, e, e / ex - 1, "%.6f" % en if nw else "–", "%+.1e" % (en / ex - 1) if nw else "–"))
            print(L[-1], file=sys.stderr, flush=True)
    L += ["", "## 2. Aligned cube, edge a/2, exact 7/192 = %.6f" % (7 / 192), "",
          "| cells per cube edge | m along | DemagFieldPBC | rel. error | Newell | rel. error |", "|---|---|---|---|---|---|"]
    for nc in (4, 8, 16, 32):
        n = 2 * nc
        ms = np.zeros((n, n, n))
        ms[:nc, :nc, :nc] = Ms
        nw = C.newell_pbc(n, lex, torch.tensor(ms[..., None]), 1) if nc <= 8 else None
        for name, u in DIRS.items():
            e = fd_energy(ms, uni(n, u), lex)
            en = C.energy_newell_pbc(*nw, uni(n, u), Ms, Kd) if nw else float("nan")
            L.append("| %d | %s | %.6f | %+.2e | %s | %s |" % (nc, name, e, e / (7 / 192) - 1,
                     "%.6f" % en if nw else "–", "%+.1e" % (en / (7 / 192) - 1) if nw else "–"))
            print(L[-1], file=sys.stderr, flush=True)
    src = subprocess.run(["git", "show", "fde0c62:magnumnp/field_terms/demagPBC.py"], cwd=HERE, capture_output=True,
                         text=True, check=True).stdout.replace("from .field_terms", "from magnumnp.field_terms.field_terms")
    up = {}
    exec(compile(src, "upstream_demagPBC", "exec"), up)
    N, dx = 48, 3 * lex
    chi = fcc_box(N, d_used * lex, dx)["ids"] >= 0
    st = State(Mesh((N, N, N), (dx,) * 3, pbc=(1, 1, 1)))
    st.material = {"Ms": torch.tensor(np.where(chi, Ms, 1e-8 * Ms)[..., None])}
    r = torch.randn((N, N, N, 3), generator=torch.Generator().manual_seed(1), dtype=torch.float64)
    st.m = r / r.norm(dim=-1, keepdim=True)
    h1, h2 = DemagFieldPBC().h(st), up["DemagFieldPBC"]().h(st)
    L += ["", "## 3. Our DemagFieldPBC against the upstream code (git fde0c62)", "",
          "Random m on the B0b geometry, dx 3: max |h_ours − h_upstream| / max |h_upstream| = %.1e."
          % (float((h1 - h2).abs().max()) / float(h2.abs().max()))]
    print("\n".join(L))


if __name__ == "__main__":
    main()
