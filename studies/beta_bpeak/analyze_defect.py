#!/usr/bin/env python
"""
Unresolved structures in a magnetization state (corner test, K3_1_nx, Chris 2026-10-18). A pair of neighbour cells of
the same particle with an angle above a threshold is "high-angle". The script finds these pairs, groups their cells in
clusters (6-connected), and gives per cluster: the particle, the size and shape (PCA of the cell positions), the
distance to the particle surface (cells), the share of pairs across a cube-grain boundary (different cube axes), and the
topological charge Q on the surface of a cube of half-edge r around the cluster centre (Q = ±1: a Bloch point inside;
Q = 0: no point singularity, e.g. a vortex line that enters and leaves, or a sharp wall). Q needs the whole surface in
the same particle; the share of surface cells in the particle is given.

    python analyze_defect.py --config runs/K3_1_nx/config.json --state K3_1_nx/ckpt_03_amp03.pt [--state ...]
"""
import argparse
import json
import math
import pathlib

import numpy as np
from scipy import ndimage

import run_loops as RL
from anisotropy_model import box_edge
from units import units

THRESHOLDS = (90.0, 120.0, 150.0, 170.0)


def geometry(cfg):
    """(ids, cube axes, N) of the run of config cfg, as run_loops builds them."""
    lex = units(cfg["Js"], cfg["A"])["l_ex"]
    a_lex, d_lex = box_edge(cfg["d_lex"], cfg["phi"], cfg.get("a_quant"))
    N = int(round(a_lex / cfg["dx_lex"]))
    a_lex = N * cfg["dx_lex"]
    geo = RL.locked_geometry(N, d_lex * lex, cfg["dx_lex"] * lex, 1, sub=8 if cfg.get("surface") == "fraction" else 0)
    ids = geo["ids"]
    ax, _ = RL.locked_cube_axes(ids, 1, cfg["dx_lex"], a_lex, d_lex, cfg["Leff_lex"], cfg["seed"])
    return ids, ax, N


def pair_angles(m, ids):
    """For each axis: (angle in degrees [N,N,N] between cell i and its neighbour i − 1 along the axis, mask of the
    pairs inside one particle); periodic, as the diagnostics of run_loops."""
    out = []
    for a in range(3):
        c = np.clip((m * np.roll(m, 1, axis=a)).sum(-1), -1.0, 1.0)
        same = (ids >= 0) & (ids == np.roll(ids, 1, axis=a))
        out.append((np.degrees(np.arccos(c)), same))
    return out


def solid_angle(m1, m2, m3):
    num = np.einsum("...i,...i", m1, np.cross(m2, m3))
    den = 1.0 + np.einsum("...i,...i", m1, m2) + np.einsum("...i,...i", m2, m3) + np.einsum("...i,...i", m3, m1)
    return 2.0 * np.arctan2(num, den)


def charge(m, ids, c, r, pid):
    """Topological charge on the surface of the cube [c − r, c + r]^3 (periodic indices) and the share of its surface
    cells in particle pid."""
    N = m.shape[0]
    q, n_in, n_all = 0.0, 0, 0
    rng = np.arange(-r, r + 1)
    for ax in range(3):
        for side in (-1, 1):
            u, v = [k for k in range(3) if k != ax]
            idx = np.zeros((len(rng), len(rng), 3), int)
            idx[..., ax] = c[ax] + side * r
            idx[..., u] = c[u] + rng[:, None]
            idx[..., v] = c[v] + rng[None, :]
            idx %= N
            M = m[idx[..., 0], idx[..., 1], idx[..., 2]]
            I = ids[idx[..., 0], idx[..., 1], idx[..., 2]]
            n_in += int((I == pid).sum())
            n_all += I.size
            # outward orientation: (u, v, ax) right-handed for side +1
            sgn = side * (1 if (u, v, ax) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)) else -1)
            a, b, cc, d = M[:-1, :-1], M[1:, :-1], M[1:, 1:], M[:-1, 1:]
            q += sgn * (solid_angle(a, b, cc) + solid_angle(a, cc, d)).sum()
    return q / (4.0 * math.pi), n_in / n_all


def analyze(m, ids, ax_cube, thr=150.0):
    pa = pair_angles(m, ids)
    counts = {t: int(sum(((ang > t) & same).sum() for ang, same in pa)) for t in THRESHOLDS}
    cells = np.zeros(ids.shape, bool)
    gb, n_hi = 0, 0
    for a, (ang, same) in enumerate(pa):
        hi = (ang > thr) & same
        cells |= hi
        cells |= np.roll(hi, -1, axis=a)
        n_hi += int(hi.sum())
        diff = np.abs((ax_cube * np.roll(ax_cube, 1, axis=a)).sum(-1)) < 0.999
        gb += int((hi & diff).sum())
    lab, n = ndimage.label(cells)
    out = {"counts": counts, "n_pairs_thr": n_hi, "share_grain_boundary": gb / max(n_hi, 1), "clusters": []}
    for k in range(1, n + 1):
        pos = np.argwhere(lab == k)
        pid = int(np.bincount(ids[tuple(pos.T)][ids[tuple(pos.T)] >= 0]).argmax())
        inside = ids == pid
        dist = ndimage.distance_transform_edt(np.pad(inside, 8, mode="wrap"))[8:-8, 8:-8, 8:-8]
        cen = np.round(pos.mean(0)).astype(int)
        ev = np.sort(np.linalg.eigvalsh(np.cov(pos.T))) if len(pos) > 3 else np.zeros(3)
        qs = [charge(m, ids, cen, r, pid) for r in (2, 3, 4)]
        out["clusters"].append({"cells": len(pos), "particle": pid, "centre": cen.tolist(),
                                "dist_to_surface_cells": float(dist[tuple(pos.T)].mean()),
                                "pca_sd_cells": [round(float(math.sqrt(max(e, 0))), 2) for e in ev],
                                "Q_r2_r3_r4": [round(q, 3) for q, _ in qs], "in_particle_r2_r3_r4": [round(f, 2) for _, f in qs]})
    out["clusters"].sort(key=lambda c: -c["cells"])
    return out


def main():
    import torch
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--state", action="append", required=True)
    ap.add_argument("--thr", type=float, default=150.0)
    a = ap.parse_args()
    cfg = json.loads(pathlib.Path(a.config).read_text())
    ids, ax_cube, N = geometry(cfg)
    for s in a.state:
        m = torch.load(s, map_location="cpu")["m"].double().numpy()
        m /= np.maximum(np.linalg.norm(m, axis=-1, keepdims=True), 1e-30)
        r = analyze(m, ids, ax_cube, a.thr)
        print("== %s  (N %d)" % (s, N))
        print("   pairs in one particle above " + ", ".join("%g°: %d" % (t, n) for t, n in r["counts"].items()))
        print("   pairs above %g°: %d, share across a cube-grain boundary %.2f" % (a.thr, r["n_pairs_thr"],
                                                                             r["share_grain_boundary"]))
        for c in r["clusters"][:6]:
            print("   cluster: %(cells)d cells, particle %(particle)d, centre %(centre)s, mean distance to the surface "
                  "%(dist_to_surface_cells).1f cells, PCA SD %(pca_sd_cells)s cells, Q (r 2/3/4) %(Q_r2_r3_r4)s, "
                  "surface in the particle %(in_particle_r2_r3_r4)s" % c)


if __name__ == "__main__":
    main()
