#!/usr/bin/env python
"""
Counts of the preparation of PLAN 4.21 (PROTOKOLL 7.146), written to results/plan421_scan.json:
  contacts  for each production point of PLAN 4.17 (jobs/sobol_main.txt: 96 Sobol points and 8 centre points; the 4
            corners of jobs/sobol_corners.txt) with a gap of at most 3.5 cells: the cell pairs of two particles with
            f > 0 that touch at a face, an edge or a corner, at q 8 and q 16 (the sub-cell sampling of geometry.fcc_box,
            without its contact rule), and the cells with f > 0 for two particles. Above 3.5 cells no contact is
            possible: a cut-cell centre lies within R + 0.866 dx, and neighbour centres are at most 1.732 dx apart.
  exchange  V2 against the exchange of García-Cervera et al. (2003) at the 4.16 geometry (148^3; q 8 and q 16): the
            outer faces of the inner cells (centre inside) to a cell with f > 0 and the centre outside (the V2 ghost)
            or to an empty cell (no ghost); the weight of the ghost pair in the f-weighted exchange, 2 f_i f_j /
            (f_i + f_j) of the full A (median, P10); the inner cut cells and the weight of their pairs.
    python scan421.py
"""
import json
import math
import pathlib
import shlex

import numpy as np

import run_loops as RL
from anisotropy_model import box_edge
from geometry import FCC_SITES

HERE = pathlib.Path(__file__).resolve().parent


def cells_with_f(N, d, dx, sub, chunk=4000):
    """owner [4, N, N, N]: the cells with f > 0 of each particle (the arithmetic of geometry.fcc_box)."""
    a = N * dx
    R = 0.5 * d
    x = (np.arange(N) + 0.5) * dx
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    q = (np.arange(sub) + 0.5) / sub - 0.5
    QX, QY, QZ = [g.reshape(-1) * dx for g in np.meshgrid(q, q, q, indexing="ij")]
    owner = np.zeros((len(FCC_SITES), N, N, N), bool)
    for p, s in enumerate(FCC_SITES * a):
        c = []
        for C, sc in zip((X, Y, Z), s):
            dd = C - sc
            dd -= a * np.round(dd / a)
            c.append(dd)
        r = np.sqrt(c[0] ** 2 + c[1] ** 2 + c[2] ** 2)
        half = 0.5 * math.sqrt(3.0) * dx
        edge = np.abs(r - R) < half
        e0 = [cc[edge] for cc in c]
        fe = np.zeros(len(e0[0]), bool)
        for i in range(0, len(fe), chunk):
            s2 = sum((e0[k][i:i + chunk, None] + Q) ** 2 for k, Q in enumerate((QX, QY, QZ)))
            fe[i:i + chunk] = (s2 <= R * R).mean(axis=1) > 0.0
        own = r <= R - half
        own[edge] = fe
        owner[p] = own
    return owner


def contacts(owner):
    ids = np.where(owner.any(0), owner.argmax(0), -1)
    n = {1: 0, 2: 0, 3: 0}
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            for sz in (-1, 0, 1):
                k = abs(sx) + abs(sy) + abs(sz)
                if k:
                    nb = np.roll(ids, (sx, sy, sz), axis=(0, 1, 2))
                    n[k] += int(((ids >= 0) & (nb >= 0) & (ids != nb)).sum())
    return {"face": n[1] // 2, "edge": n[2] // 2, "corner": n[3] // 2, "double": int((owner.sum(0) > 1).sum())}


def mesh(rest):
    a = RL.parse_args(shlex.split(rest) + ["--out", "unused"])
    lex = RL.units(a.Js, a.A)["l_ex"]
    a_lex, d_lex = box_edge(a.d_lex, a.phi, a.a_quant)
    N = int(round(a_lex / a.dx_lex))
    return N, d_lex * lex, a.dx_lex * lex, d_lex, a.phi


def scan_contacts():
    lines = [l.split(None, 1) for l in (HERE / "jobs" / "sobol_main.txt").read_text().splitlines() if l.strip()]
    corners = {}
    for l in (HERE / "jobs" / "sobol_corners.txt").read_text().splitlines():
        if l.strip():
            n, rest = l.split(None, 1)
            corners.setdefault(n.split("_")[0], rest)
    lines += [[c, r] for c, r in sorted(corners.items())]
    out = []
    for name, rest in lines:
        N, d, dx, d_lex, phi = mesh(rest)
        rec = {"name": name, "N": N, "d_lex": round(d_lex, 2), "phi": phi, "gap_cells": round((N * dx / math.sqrt(2) - d) / dx, 3)}
        if rec["gap_cells"] <= 3.5:
            for sub in (8, 16):
                rec["q%d" % sub] = contacts(cells_with_f(N, d, dx, sub))
        out.append(rec)
    touch = lambda r, q: ("q%d" % q in r) and sum(r["q%d" % q].values()) > 0
    return {"points": out, "n_points": len(out), "n_gap_le_3.5": sum(r["gap_cells"] <= 3.5 for r in out),
            "n_contact_q16": sum(touch(r, 16) for r in out), "n_contact_q8": sum(touch(r, 8) for r in out),
            "n_face_contact": sum(r.get("q%d" % q, {}).get("face", 0) > 0 for r in out for q in (8, 16))}


def scan_exchange():
    line = [l for l in (HERE / "jobs" / "hf416.txt").read_text().splitlines() if l.startswith("HFD_")][0]
    N, d, dx, _, _ = mesh(line.split(None, 1)[1])
    from geometry import fcc_box
    res = {"N": N}
    for sub in (8, 16):
        g = fcc_box(N, d, dx, sub=sub)
        ids, frac = g["ids"], g["frac"].astype(np.float64)
        inside = (fcc_box(N, d, dx)["ids"] >= 0) & (ids >= 0)
        ghost = (ids >= 0) & ~inside
        n_g = n_v = 0
        w = []
        for ax in range(3):
            for sh in (1, -1):
                nb_ids, nb_f = np.roll(ids, -sh, ax), np.roll(frac, -sh, ax)
                ext = inside & ~(np.roll(inside, -sh, ax) & (nb_ids == ids))
                gf = ext & np.roll(ghost, -sh, ax) & (nb_ids == ids)
                n_g += int(gf.sum())
                n_v += int((ext & ~gf).sum())
                w.append(2 * frac[gf] * nb_f[gf] / (frac[gf] + nb_f[gf]))
        w = np.concatenate(w)
        icut = inside & (frac < 1)
        wi = []
        for ax in range(3):
            nb_in, nb_f, nb_ids = np.roll(inside, -1, ax), np.roll(frac, -1, ax), np.roll(ids, -1, ax)
            sel = inside & nb_in & (nb_ids == ids) & ((frac < 1) | (nb_f < 1))
            wi.append(2 * frac[sel] * nb_f[sel] / (frac[sel] + nb_f[sel]))
        wi = np.concatenate(wi)
        res["q%d" % sub] = {"inner_cells": int(inside.sum()), "ghost_cells": int(ghost.sum()),
                            "outer_faces": n_g + n_v, "faces_to_ghost": n_g, "faces_without_ghost": n_v,
                            "ghost_weight_median": float(np.median(w)), "ghost_weight_P10": float(np.percentile(w, 10)),
                            "ghost_weight_P90": float(np.percentile(w, 90)), "inner_cut_cells": int(icut.sum()),
                            "inner_cut_pair_weight_median": float(np.median(wi))}
    return res


if __name__ == "__main__":
    res = {"contacts": scan_contacts(), "exchange": scan_exchange()}
    (HERE / "results" / "plan421_scan.json").write_text(json.dumps(res, indent=1))
    c, e = res["contacts"], res["exchange"]
    print("contacts: %d points, %d with gap <= 3.5 cells, contacts at q 16: %d, at q 8: %d, face contacts: %d" % (
        c["n_points"], c["n_gap_le_3.5"], c["n_contact_q16"], c["n_contact_q8"], c["n_face_contact"]))
    for q in (8, 16):
        x = e["q%d" % q]
        print("exchange q %d: outer faces %d, without ghost %d (%.1f %%), ghost weight median %.3f P10 %.3f; inner cut "
              "cells %d, pair weight median %.3f" % (q, x["outer_faces"], x["faces_without_ghost"],
                                                     100 * x["faces_without_ghost"] / x["outer_faces"],
                                                     x["ghost_weight_median"], x["ghost_weight_P10"],
                                                     x["inner_cut_cells"], x["inner_cut_pair_weight_median"]))
