#!/usr/bin/env python
"""
Hot spots of the loss outside the drive frequency f (Chris 2026-10-22; method agreed with Chris before the results):
which structure sits where the loss outside f is largest, and how it moves over a cycle. Descriptive.
Input: a run folder with config.json, geometry.vti, the accumulators acc_<i>_<stage>.npz and the snapshots
m_<stage>_c<k>_ph<deg>.vti (16 phases x 2 cycles) of the stages with snapshots (9, 50 and 150 mT).
Per stage:
  1. Hot spots. Loss outside f per cell q = p - p_f (p: Gilbert loss of the cell, p_f: the loss of its motion at f,
     as acc_reduce.py). Hot cells: q >= the 99th percentile of q over the particle cells. Regions: connected hot
     cells (26 neighbours, periodic box; the particles do not touch). Main regions: the regions with the largest q
     that together carry at least half of q of the hot cells (at most MAX_REGIONS).
  2. Structure type of a region, from the gradient tensor G_ij = sum_k (d_i m_k)(d_j m_k), summed over the region
     cells dilated by 2 cells (in the particle) and over all snapshots of the stage; eigenvalues l1 >= l2 >= l3:
       wall  l2 / l1 < T_TYPE                      (m changes along 1 direction)
       line  l2 / l1 >= T_TYPE and l3 / l1 < T_TYPE (along 2 directions: vortex core or Bloch line)
       3D    l3 / l1 >= T_TYPE                      (along 3 directions)
     BP: a Bloch point (analyze_defect.bloch_points) in the window of the region in at least one snapshot.
     T_TYPE is fixed on synthetic states (test_study.py: wall, vortex, Bloch point, each with ripple) before the run
     data. A ripple alone (no wall, no vortex) also gives 3D; thus 3D without BP means an irregular texture.
  3. Motion over the cycle, in the window (bounding box of the region + WINDOW cells): the series of the mean m of
     the region cells and of the centroid of |grad m|² in the window (the position of the structure), each projected
     on the direction of its largest variation over the snapshots.
       J    largest change between two neighbouring phases of a cycle / range in that cycle (cyclic), mean of the
            cycles. A sine sampled at 16 phases gives 0.19, one step up and one step down give 1. smooth: J <= J_SMOOTH,
            step-like: J >= J_STEP.
       rep  rms of (cycle 2 - cycle 1) / mean range of the cycles; the motion repeats if rep < REP_MAX.
       amp  range of the centroid along its direction, l_ex (mean of the cycles).
     Snapshots resolve motion down to T / 16; a faster motion shows as a step.
Reference: the same for the hot cells of the loss at f (p_f >= its 99th percentile).
Readings per run and stage, as shares of q (of p_f for the reference) of the main regions: a type "dominates" with
at least 50 %, else "mixed"; the motion is "step-like" if the regions with J >= J_STEP carry at least 50 %, "smooth"
if the regions with J <= J_SMOOTH carry at least 50 %, else "mixed"; it "repeats" if the regions with rep < REP_MAX
carry at least 50 %.

    python eval_hotspots.py [--json FILE] RUN_DIR ... > report.md
"""
import argparse
import itertools
import json
import math
import pathlib
import re

import numpy as np
from scipy import ndimage

MU0, GAMMA = 4e-7 * math.pi, 2.21276157e5
HOT_Q = 0.99                 # hot cells: the 1 % with the largest loss outside f
MAX_REGIONS = 20
DILATE = 2                   # cells around a region for its gradient tensor
WINDOW = 4                   # cells around the bounding box of a region for the centroid
T_TYPE = 0.3                 # eigenvalue ratio of the structure types
J_SMOOTH, J_STEP, REP_MAX = 0.3, 0.5, 0.2
TYPES = ("wall", "line", "3D")


# --- geometry and fields -------------------------------------------------------------------------------------------

def grad_dirs(m, ids):
    """d_i m per cell (per cell length), shape (3, N, N, N, 3): central difference where both face neighbours are in
    the same particle, else the one-sided difference that exists, else 0."""
    out = np.zeros((3,) + m.shape, np.float32)
    inside = ids >= 0
    for ax in range(3):
        mf, mb = np.roll(m, -1, ax), np.roll(m, 1, ax)
        f = (np.roll(ids, -1, ax) == ids) & inside
        b = (np.roll(ids, 1, ax) == ids) & inside
        out[ax] = np.where((f & b)[..., None], 0.5 * (mf - mb),
                           np.where(f[..., None], mf - m, np.where(b[..., None], m - mb, 0.0)))
    return out


def tensor(D, sel):
    """Gradient tensor G_ij = sum over the cells sel of sum_k (d_i m_k)(d_j m_k); D from grad_dirs."""
    v = D[:, sel, :]                                            # (3, n, 3)
    return np.einsum("ink,jnk->ij", v, v, dtype=np.float64)


def structure_type(G):
    """(type, l2 / l1, l3 / l1) from a gradient tensor."""
    lam = np.sort(np.linalg.eigvalsh(G))[::-1]
    if lam[0] <= 0:
        return "none", 0.0, 0.0
    r2, r3 = float(lam[1] / lam[0]), float(lam[2] / lam[0])
    return ("3D" if r3 >= T_TYPE else "line" if r2 >= T_TYPE else "wall"), r2, r3


def label_periodic(mask):
    """Connected components of mask with 26 neighbours in a periodic box: (labels, n)."""
    lab, n = ndimage.label(mask, np.ones((3, 3, 3), bool))
    parent = np.arange(n + 1)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for sh in itertools.product((-1, 0, 1), repeat=3):
        if sh == (0, 0, 0):
            continue
        nb = np.roll(lab, sh, axis=(0, 1, 2))
        pairs = (lab > 0) & (nb > 0) & (lab != nb)               # only pairs across the box faces remain
        for a, b in set(zip(lab[pairs].tolist(), nb[pairs].tolist())):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
    root = np.array([find(a) for a in range(n + 1)])
    _, inv = np.unique(root, return_inverse=True)
    return inv.reshape(-1)[lab].reshape(lab.shape), int(inv.max())


def dilate_in(sel, part, k):
    """sel dilated by k cells (26 neighbours, periodic), limited to part."""
    out = sel.copy()
    for _ in range(k):
        grow = out.copy()
        for sh in itertools.product((-1, 0, 1), repeat=3):
            grow |= np.roll(out, sh, axis=(0, 1, 2))
        out = grow & part
    return out


def window_index(sel, N, pad):
    """Index arrays (ix, iy, iz) of the periodic bounding box of the cells sel plus pad cells, and the offsets that
    unwrap them to a local, contiguous frame."""
    idx = np.argwhere(sel)
    ref = idx[0]
    loc = ((idx - ref + N // 2) % N) - N // 2 + ref              # minimum image around the first cell
    lo, hi = loc.min(0) - pad, loc.max(0) + pad
    return tuple(np.arange(l, h + 1) % N for l, h in zip(lo, hi)), lo


# --- series measures ---------------------------------------------------------------------------------------------

def principal(series):
    """Projection of a vector series (n, d) on the direction of its largest variation."""
    x = np.asarray(series, float)
    x = x - x.mean(0)
    if not np.any(x):
        return np.zeros(len(x))
    _, _, vt = np.linalg.svd(x, full_matrices=False)
    return x @ vt[0]


def step_measure(s, per_cycle):
    """J (mean over the cycles of the largest cyclic neighbour change / range), rep (rms of cycle 2 - cycle 1 / mean
    range) and the mean range; s: the series of all cycles, per_cycle phases each."""
    c = np.asarray(s, float).reshape(-1, per_cycle)
    rng = c.max(1) - c.min(1)
    J = [float(np.abs(np.diff(np.append(x, x[0]))).max() / r) if r > 0 else float("nan") for x, r in zip(c, rng)]
    rep = float(np.sqrt(np.mean((c[1] - c[0]) ** 2)) / rng.mean()) if len(c) > 1 and rng.mean() > 0 else float("nan")
    return float(np.nanmean(J)), rep, float(rng.mean())


# --- one stage ---------------------------------------------------------------------------------------------------

def snapshot_files(d, stage):
    """{cycle: [files ordered by phase]} of one stage."""
    out = {}
    for f in d.glob("m_%s_c*_ph*.vti" % stage):
        c, ph = map(int, re.match(r"m_.+_c(\d+)_ph(\d+)\.vti", f.name).groups())
        out.setdefault(c, []).append((ph, f))
    return {c: [f for _, f in sorted(v)] for c, v in sorted(out.items())}


def regions(q, ids, mask, hot_q=HOT_Q):
    """Main regions of the cell values q (on the particle cells): list of (selection, share of q of all particle
    cells), the hot share and the number of hot regions."""
    full = np.zeros(ids.shape)
    full[mask] = q
    hot = mask & (full >= np.quantile(q, hot_q))
    lab, n = label_periodic(hot)
    qs = ndimage.sum(np.where(hot, full, 0.0), lab, index=np.arange(1, n + 1))
    order = np.argsort(qs)[::-1]
    tot, hot_sum = q.sum(), qs.sum()
    out, cum = [], 0.0
    for i in order[:MAX_REGIONS]:
        out.append((lab == i + 1, float(qs[i] / tot)))
        cum += qs[i]
        if cum >= 0.5 * hot_sum:
            break
    return out, float(hot_sum / tot), n


def stage_record(d, stage, acc, ids, cfg, defects):
    """Hot regions of q (loss outside f) and of p_f (reference) of one stage, with type and motion."""
    import domains
    mask = ids >= 0
    N = ids.shape[0]
    alpha, Ms = cfg["alpha"], cfg["Js"] / MU0
    z = np.load(acc)
    w = 2 * math.pi * float(z["freq"])
    p = z["p"].astype(np.float64)
    pf = alpha * MU0 * Ms / GAMMA * w * w * (np.abs(z["m1"].astype(np.complex128)) ** 2).sum(-1) / 2.0
    from dissmap_batch import radii
    rR = radii(ids, cfg["dx_lex"], N)
    dist = (1.0 - rR) * 0.5 * cfg["d_lex_used"] / cfg["dx_lex"]   # cells to the surface (the spheres do not touch)
    sets = {}
    for name, val in (("nonf", p - pf), ("f", pf)):
        regs, hot_share, n_hot = regions(val, ids, mask)
        full = np.zeros(ids.shape)
        full[mask] = val
        R = []
        for sel, share in regs:
            k = int(np.bincount(ids[sel]).argmax())
            wts = full[sel]
            R.append({"sel": sel, "big": dilate_in(sel, ids == k, DILATE), "win": window_index(sel, N, WINDOW),
                      "share": share, "cells": int(sel.sum()), "particle": k,
                      "r_R": float(np.average(rR[sel], weights=np.maximum(wts, 0) + 1e-30)),
                      "surface_cells": float(np.average(dist[sel], weights=np.maximum(wts, 0) + 1e-30)),
                      "G": np.zeros((3, 3)), "m": [], "pos": [], "bp": []})
        sets[name] = {"hot_share": hot_share, "n_hot_regions": n_hot, "regions": R}
    files = snapshot_files(d, stage)
    per_cycle = {len(v) for v in files.values()}
    assert len(per_cycle) == 1, files.keys()
    per_cycle = per_cycle.pop()
    for c, fl in files.items():
        for f in fl:
            st = domains.load_state(f)
            assert np.array_equal(st["ids"], ids), f
            m = st["m"].astype(np.float32)
            D = grad_dirs(m, ids)
            e = (D ** 2).sum(axis=(0, -1))
            bps = defects.bloch_points(m, ids) if defects else []
            for s in sets.values():
                for r in s["regions"]:
                    r["G"] += tensor(D, r["big"])
                    r["m"].append(m[r["sel"]].mean(0))
                    (ix, iy, iz), lo = r["win"]
                    ew = e[np.ix_(ix, iy, iz)] * (ids[np.ix_(ix, iy, iz)] == r["particle"])
                    grid = np.stack(np.meshgrid(*[np.arange(len(a)) for a in (ix, iy, iz)], indexing="ij"), -1)
                    r["pos"].append((grid * ew[..., None]).sum((0, 1, 2)) / max(ew.sum(), 1e-30) * cfg["dx_lex"])
                    sx, sy, sz = (set(a.tolist()) for a in (ix, iy, iz))
                    r["bp"].append(sum(1 for b in bps if b[2][0] in sx and b[2][1] in sy and b[2][2] in sz))
    out = {}
    for name, s in sets.items():
        regs = []
        for r in s["regions"]:
            t, r2, r3 = structure_type(r["G"])
            Jm, repm, _ = step_measure(principal(r["m"]), per_cycle)
            Jp, repp, amp = step_measure(principal(r["pos"]), per_cycle)
            regs.append({"share": r["share"], "cells": r["cells"], "particle": r["particle"], "r_R": r["r_R"],
                         "surface_cells": r["surface_cells"], "type": t, "l2_l1": r2, "l3_l1": r3,
                         "bp_snapshots": int(sum(1 for b in r["bp"] if b)), "bp_max": int(max(r["bp"], default=0)),
                         "J_m": Jm, "rep_m": repm, "J_pos": Jp, "rep_pos": repp, "amp_pos_lex": amp})
        out[name] = {"hot_share": s["hot_share"], "n_hot_regions": s["n_hot_regions"], "regions": regs,
                     "summary": summary(regs)}
    out["snapshots"] = {"cycles": list(files), "per_cycle": per_cycle}
    return out


def summary(regs):
    """Shares of the main-region loss by type, motion and repetition, and the readings."""
    tot = sum(r["share"] for r in regs)
    sh = lambda f: float(sum(r["share"] for r in regs if f(r)) / tot) if tot > 0 else float("nan")
    s = {"type_" + t: sh(lambda r, t=t: r["type"] == t) for t in TYPES}
    s["bp"] = sh(lambda r: r["bp_snapshots"] > 0)
    for k in ("m", "pos"):
        s["step_" + k] = sh(lambda r, k=k: r["J_" + k] >= J_STEP)
        s["smooth_" + k] = sh(lambda r, k=k: r["J_" + k] <= J_SMOOTH)
        s["repeats_" + k] = sh(lambda r, k=k: r["rep_" + k] < REP_MAX)
        s["reading_motion_" + k] = ("step-like" if s["step_" + k] >= 0.5 else
                                    "smooth" if s["smooth_" + k] >= 0.5 else "mixed")
        s["reading_repeats_" + k] = s["repeats_" + k] >= 0.5
    best = max(TYPES, key=lambda t: s["type_" + t])
    s["reading_type"] = best if s["type_" + best] >= 0.5 else "mixed"
    return s


def run(d):
    import domains
    import analyze_defect
    d = pathlib.Path(d)
    cfg = json.loads((d / "config.json").read_text())
    ids = domains._ids_from_config(cfg).astype(np.int8)
    out = {"run": d.name, "stages": {}}
    for acc in sorted(d.glob("acc_*_*.npz")):
        stage = re.match(r"acc_\d+_(.+)\.npz", acc.name).group(1)
        if snapshot_files(d, stage):
            out["stages"][stage] = stage_record(d, stage, acc, ids, cfg, analyze_defect)
    return out


def report(res):
    L = ["# Hot spots of the loss outside f (eval_hotspots.py)", "",
         "Hot cells: the 1 %% of the particle cells with the largest loss outside f (q = p - p_f). Main regions: the "
         "regions that carry half of q of the hot cells. Shares are shares of q of all particle cells. Type from the "
         "gradient tensor (wall: 1 direction, line: 2, 3D: 3; threshold %.1f). J: largest step between neighbouring "
         "phases / range (sine 0.19, step 1); rep: cycle 2 - cycle 1 / range. Reference: the same for the loss at f."
         % T_TYPE, ""]
    for r in res:
        for stage, s in r["stages"].items():
            for name, lab in (("nonf", "outside f"), ("f", "at f (reference)")):
                x = s[name]
                sm = x["summary"]
                L += ["## %s %s, loss %s" % (r["run"], stage, lab), "",
                      "Hot cells carry %.2f of the loss; %d hot regions; the %d main regions carry %.2f. Readings: type "
                      "%s; motion of m %s, of the position %s; repeats (m / position) %s / %s." % (
                          x["hot_share"], x["n_hot_regions"], len(x["regions"]),
                          sum(g["share"] for g in x["regions"]), sm["reading_type"], sm["reading_motion_m"],
                          sm["reading_motion_pos"], sm["reading_repeats_m"], sm["reading_repeats_pos"]),
                      "Shares of the main-region loss: wall %.2f, line %.2f, 3D %.2f, with a Bloch point %.2f; "
                      "step-like m %.2f, position %.2f." % (sm["type_wall"], sm["type_line"], sm["type_3D"], sm["bp"],
                                                            sm["step_m"], sm["step_pos"]), "",
                      "| share | cells | particle | r/R | to surface, cells | type | l2/l1 | l3/l1 | BP snapshots | "
                      "J m | rep m | J pos | rep pos | amp pos, l_ex |", "|---" * 14 + "|"]
                for g in x["regions"]:
                    L.append("| %.3f | %d | %d | %.2f | %.1f | %s | %.2f | %.2f | %d | %.2f | %.2f | %.2f | %.2f | "
                             "%.1f |" % (g["share"], g["cells"], g["particle"], g["r_R"], g["surface_cells"], g["type"],
                                         g["l2_l1"], g["l3_l1"], g["bp_snapshots"], g["J_m"], g["rep_m"], g["J_pos"],
                                         g["rep_pos"], g["amp_pos_lex"]))
                L.append("")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=None)
    ap.add_argument("runs", nargs="+")
    a = ap.parse_args()
    res = [run(d) for d in a.runs]
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(res, indent=1))
    print(report(res))


if __name__ == "__main__":
    main()
