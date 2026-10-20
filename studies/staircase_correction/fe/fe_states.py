#!/usr/bin/env python
"""
FE comparison, preparation (PROTOKOLL 3.24): export a coarse state for the FE code and compute its FD reference
energies with magnum.np (local venv; the FE side runs in the magnum.pi container, which has no torch).

The quantity under test is the stray field of the coarse m (constant in each cell of the FD mesh) in the smooth
spheres. The coarse m must then be defined on every cell that touches a sphere. The state defines m on its own
magnetic cells (for X3 / T3Pc: the cells with body sub-cells at q = 2); each other touching cell takes m of the
nearest such cell of the same particle (periodic centre distance).

    python fe_states.py --run RUN_DIR --out OUT_DIR [--q 2 4] [--cache KERNEL_CACHE]

Writes OUT_DIR/m.npy [N,N,N,3] (float64, zero outside the touching cells), ids.npy (voxel particle index at dx, -1:
void), touch.npy (particle index of the touching cells, -1: void), occ_q<q>.npy (body of the q-staircase) and meta.json (a, dx, d, sites, Ms, K_d and the FD
energies in units of K_d V_box):
    E_vox            voxel body (cell centre in a sphere), exact for full cubes up to the Newell PBC sum
    E_vox_100/111    the same body with uniform m; target φ_vox (1 − φ_vox) / 3 (reading R1)
    E_full           the whole box magnetic, m of the nearest touching cell in the void (no surface; m_full.npy)
    E_X3_q<q>        exact field of the q-staircase for the coarse m (DemagStaircaseFine)
    E_T3Pc_q<q>      production model: Newell with M_s v plus the pairwise correction (p 3, centroid term)
"""
import argparse
import json
import os
import pathlib
import sys
import time

import numpy as np
import torch

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent.parent / "beta_bpeak"))
os.environ.setdefault("CUDA_DEVICE", "-1")

from magnumnp import set_precision  # noqa: E402
set_precision("double")
torch.set_num_threads(int(os.environ.get("STAIRCASE_THREADS", os.cpu_count())))

from magnumnp import Mesh, State, DemagField, DemagStaircaseFine, DemagStaircaseCorrection  # noqa: E402
from magnumnp.field_terms.demag_staircase import staircase_cells  # noqa: E402
from anisotropy_model import box_edge  # noqa: E402
from geometry import fcc_box, FCC_SITES  # noqa: E402
from units import units  # noqa: E402

MU0 = 4e-7 * np.pi


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def touching(N, dx, d, a):
    """Particle index of each cell that has a part inside a sphere (exact cube-ball test, minimum image), -1: void."""
    x0 = np.arange(N) * dx
    out = -np.ones((N, N, N), dtype=np.int64)
    R = 0.5 * d
    for p, s in enumerate(FCC_SITES * a):
        dist2 = 0.0
        for ax in range(3):
            lo = x0 - s[ax]
            lo -= a * np.round((lo + 0.5 * dx) / a)          # minimum image of the cell centre
            hi = lo + dx
            g = np.where(lo > 0, lo, np.where(hi < 0, -hi, 0.0))   # distance of the interval to 0
            shp = [1, 1, 1]
            shp[ax] = N
            dist2 = dist2 + (g ** 2).reshape(shp)
        hit = dist2 < R * R
        if np.any(hit & (out >= 0)):
            raise ValueError("a cell touches two spheres")
        out[hit] = p
    return out


def fill_nearest(m, have, touch):
    """m on the cells of touch that are not in have: m of the nearest cell of have of the same particle."""
    from scipy.spatial import cKDTree
    N = touch.shape[0]
    xyz = np.stack(np.unravel_index(np.arange(N ** 3), (N,) * 3), -1).astype(np.float64)
    mf = m.reshape(-1, 3).copy()
    hv, tc = have.reshape(-1), touch.reshape(-1)
    mf[tc < 0] = 0.0
    n_new = 0
    for p in range(4):
        src = np.flatnonzero((hv == p))
        dst = np.flatnonzero((tc == p) & (hv != p))
        if len(dst) == 0:
            continue
        _, k = cKDTree(xyz[src], boxsize=N).query(xyz[dst], k=1)
        mf[dst] = mf[src[k]]
        n_new += len(dst)
    return mf.reshape(N, N, N, 3), n_new


def energy(term, N, dx, Ms_cell, m):
    state = State(Mesh((N, N, N), (dx, dx, dx), pbc=(4, 4, 4)))
    state.material = {"Ms": torch.tensor(Ms_cell[..., None], dtype=state.dtype)}
    state.m = torch.tensor(m, dtype=state.dtype)
    return float(term.E(state))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="run folder with init.pt and config.json (an X3 or T3Pc run)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--q", type=int, nargs="*", default=[2])
    ap.add_argument("--cache", default=os.environ.get("STAIRCASE_KERNEL_CACHE"))
    ap.add_argument("--state", default="init.pt", help="init.pt (prepared state) or checkpoint.pt")
    a = ap.parse_args()
    run, out = pathlib.Path(a.run), pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((run / "config.json").read_text())
    U = units(cfg["Js"], cfg["A"])
    Ms, Kd, lex = U["Ms"], U["Kd"], U["l_ex"]
    a_lex, d_lex = box_edge(cfg["d_lex"], cfg["phi"], cfg["a_quant"])
    N = int(round(a_lex / cfg["dx_lex"]))
    dx, d, a_box = cfg["dx_lex"] * lex, d_lex * lex, a_lex * lex
    m = torch.load(str(run / a.state), map_location="cpu", weights_only=False)["m"].double().numpy()
    assert m.shape == (N, N, N, 3)
    ids = fcc_box(N, d, dx)["ids"].astype(np.int64)
    _, pid2, _ = staircase_cells(fcc_box(2 * N, d, dx / 2)["ids"], 2)
    n2 = staircase_cells(fcc_box(2 * N, d, dx / 2)["ids"], 2)[0]
    have = np.where(n2 > 0, pid2, -1)                       # the magnetic cells of a state at q = 2
    touch = touching(N, dx, d, a_box)
    assert np.all(touch[have >= 0] == have[have >= 0]) and np.all(touch[ids >= 0] == ids[ids >= 0])
    m_ext, n_new = fill_nearest(m, have, touch)
    log("N %d, voxel cells %d, q2 cells %d, touching cells %d (%d filled from the nearest q2 cell)"
        % (N, (ids >= 0).sum(), (have >= 0).sum(), (touch >= 0).sum(), n_new))
    np.save(out / "m.npy", m_ext)
    np.save(out / "ids.npy", ids)
    np.save(out / "touch.npy", touch)
    # the same m on every cell of the box (each void cell takes m of the nearest touching cell): a body without a
    # surface, so that only the charge sheets on the inner cell faces remain (check of the FE error source)
    from scipy.spatial import cKDTree
    xyz = np.stack(np.unravel_index(np.arange(N ** 3), (N,) * 3), -1).astype(np.float64)
    tf = touch.reshape(-1)
    src_c, dst_c = np.flatnonzero(tf >= 0), np.flatnonzero(tf < 0)
    _, k = cKDTree(xyz[src_c], boxsize=N).query(xyz[dst_c], k=1)
    m_full = m_ext.reshape(-1, 3).copy()
    m_full[dst_c] = m_full[src_c[k]]
    m_full = m_full.reshape(N, N, N, 3)
    np.save(out / "m_full.npy", m_full)

    norm = Kd * a_box ** 3
    E = {}
    demag = DemagField(true_pbc=True, extrapolate=True, cache_dir=a.cache)
    vox = np.where(ids >= 0, Ms, 0.0)
    t = time.time()
    E["E_vox"] = energy(demag, N, dx, vox, m_ext) / norm
    for name, v in (("100", (1.0, 0.0, 0.0)), ("111", (1.0, 1.0, 1.0))):
        u = np.broadcast_to(np.array(v) / np.linalg.norm(v), (N, N, N, 3)).copy()
        E["E_vox_" + name] = energy(demag, N, dx, vox, u) / norm
    phi_vox = float((ids >= 0).mean())
    E["E_full"] = energy(demag, N, dx, np.full((N, N, N), Ms), m_full) / norm
    E["phi_vox"] = phi_vox
    E["R1_target"] = phi_vox * (1 - phi_vox) / 3
    log("voxel: E %.9f, uniform (100) %.9f (111) %.9f, target %.9f (%.0f s)"
        % (E["E_vox"], E["E_vox_100"], E["E_vox_111"], E["R1_target"], time.time() - t))
    for q in a.q:
        t = time.time()
        occ = fcc_box(q * N, d, dx / q)["ids"] >= 0
        nq = occ.reshape(N, q, N, q, N, q).sum(axis=(1, 3, 5))
        if np.any((nq > 0) & (touch < 0)):
            raise ValueError("a body sub-cell lies in a cell that does not touch a sphere")
        np.save(out / ("occ_q%d.npy" % q), occ)          # body of the q-staircase, for the FE side
        fine = DemagStaircaseFine(occ, q, Ms, true_pbc=True, extrapolate=True, cache_dir=a.cache)
        E["E_X3_q%d" % q] = energy(fine, N, dx, Ms * nq / q ** 3, m_ext) / norm
        E["phi_q%d" % q] = float(occ.mean())
        log("X3 q %d: E %.9f, φ %.6f (%.0f s)" % (q, E["E_X3_q%d" % q], E["phi_q%d" % q], time.time() - t))
        # T3Pc: Newell with M_s v plus the pairwise near-field correction with the centroid term (p = 3)
        t = time.time()
        n_s, pid_s, occ_s = staircase_cells(fcc_box(q * N, d, dx / q)["ids"], q)
        corr = DemagStaircaseCorrection(n_s, occ_s, q, 3, Ms, form="pair", centroid=True, demag=demag)
        Mv = Ms * n_s / q ** 3
        E["E_T3Pc_q%d" % q] = (energy(demag, N, dx, Mv, m_ext) + energy(corr, N, dx, Mv, m_ext)) / norm
        log("T3Pc q %d: E %.9f (%.0f s)" % (q, E["E_T3Pc_q%d" % q], time.time() - t))
    meta = {"run": str(run), "state": a.state, "N": N, "dx": dx, "a": a_box, "d": d, "Ms": Ms, "Kd": Kd,
            "sites": (FCC_SITES * a_box).tolist(), "n_filled": int(n_new), "energies_KdVbox": E}
    (out / "meta.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta["energies_KdVbox"], indent=1))


if __name__ == "__main__":
    main()
