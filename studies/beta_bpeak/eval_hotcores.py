#!/usr/bin/env python
"""
Hot-spot follow-up (Chris 2026-10-23, option C; method agreed with Chris before the results): for the main regions of
the loss outside f (eval_hotspots.py, the same regions), (1) the structure type in each snapshot, (2) vortex core or
Bloch line, (3) the Thiele scale: the loss of the structure against the loss of its rigid motion along the measured
core path. Descriptive. The structures on the l_ex scale are not resolved at dx 3 l_ex (PROTOKOLL section 4): the
centre measure (2a) reads one cell, the ring measure (2b) reads the cells at 2 ... 3 cells from the core.
Input: run folders as for eval_hotspots.py (config.json, geometry.vti, acc_<i>_<stage>.npz, m_<stage>_c<k>_ph<deg>.vti,
16 phases per cycle); a replay folder with 1 cycle of snapshots is also accepted.
Per main region of q = p - p_f, in its window (bounding box + WINDOW cells, cells of its particle only):
  1. Type per snapshot: eval_hotspots.structure_type of the gradient tensor of the dilated region in one snapshot;
     shares of line, wall and 3D over the snapshots. Line direction t_n: eigenvector of the smallest eigenvalue; t:
     the same of the tensor summed over the snapshots (the tensor of eval_hotspots). Turn: angle between t_n and t.
  2. Core: the window is cut into slices of 1 cell along t. Per slice and snapshot, the core cell is the cell of the
     dilated region with the largest |grad m|^2; the core point is the |grad m|^2-weighted centroid of the cells of
     the slice within CORE_R cells of it.
     (a) centre: |m . t| of the core cell. (b) ring: (max - min) / 2 of m . t over the cells of the slice at RING
     cells from the core point (at least RING_MIN cells). Vortex core: m at the centre along t, m on the ring curls
     around t (one sign of m . t, small spread). Bloch line: m at the centre across t, the ring crosses the wall and
     the domains on both sides (both signs, large spread). On synthetic states at dx 3 l_ex (--selftest) the centre
     does not separate a vortex core of width 1 l_ex (0.16 ... 0.26) from a Bloch line (0.06 ... 0.31); the ring does
     (vortex core at most 0.27, Bloch line at least 0.32), and a Bloch line of width 1 ... 3 l_ex in a wall has the
     type wall in every snapshot. Class (Chris 2026-10-23, fixed before the run data), from the median ring and the
     type shares of the snapshots: vortex core: ring < T_RING and line >= 0.5; Bloch line: ring >= T_RING and
     line >= 0.5; wall or narrow Bloch line: ring >= T_RING and wall >= 0.5; else unclear. The centre is information.
  3. Thiele scale, per cycle: core path x_s(t) per slice (the part across t), velocity v_s by the spectral
     derivative (harmonics 1 ... K); a slice whose core steps more than JUMP cells between 2 phases has jumped to
     another structure and gives no velocity (jump share). Rigid motion: dm/dt = -(v_s . grad) m in every cell of the
     window, with v_s of the slice of the cell (the nearest slice with a velocity) and grad m of the snapshot. Its
     loss: P_rigid = alpha mu0 Ms / gamma |dm/dt|^2, summed over the window cells and averaged over the 16 phases
     (all frequencies: the instantaneous values); its loss at f from harmonic 1 of the Fourier series of each cell,
     q_rigid = P_rigid minus that. Measured: P = p and q = p - p_f of the accumulators, summed over the window cells.
     Readings: P / P_rigid and q / q_rigid. Near the synthetic value: the loss is the loss of the rigid motion of the
     structure along its path; well above: additional loss inside the structure (deformation, precession) or loss of
     other motion in the window. A large q at a sine path is kinematic: a localized structure that moves far gives
     each cell a non-sinusoidal signal. Also: the loss at f and at 2f ... Kf of the snapshot series ((dm/dt)_k =
     i k w m_k) over the same from the accumulators (mk2), near 1 if the 16 phases resolve the motion, and the share
     of k >= 2 in the mean square core velocity (the path itself).
     Limits: one velocity per slice (the line moves as a whole in its slice); the window holds also other textures.

    python eval_hotcores.py [--json FILE] RUN_DIR ... > report.md
    python eval_hotcores.py --selftest            (measures 2 and 3 on synthetic states, the table of test_study.py)
"""
import argparse
import json
import math
import pathlib
import re

import numpy as np

import eval_hotspots as H

MU0, GAMMA = H.MU0, H.GAMMA
K = 7                        # harmonics, as run_loops --acc_harm 7
CORE_R = 1.5                 # cells around the core cell for the core point
RING = (2.0, 3.0)            # cells from the core point: the ring
RING_MIN = 6
JUMP = 6.0                   # cells: a larger core step between 2 phases is a jump to another structure (a sine path of
#                              18 cells peak to peak, 54 l_ex as the largest in eval_hotspots, steps 3.5 cells)
T_RING = 0.30                # between vortex cores (<= 0.27) and Bloch lines (>= 0.32) of the synthetic states


def frame(t):
    """Two unit vectors that make with t a right-handed orthonormal frame."""
    a = np.array([1.0, 0.0, 0.0]) if abs(t[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = np.cross(t, a)
    e1 /= np.linalg.norm(e1)
    return e1, np.cross(t, e1)


def smallest(G):
    return np.linalg.eigh(G)[1][:, 0]


def classify(ring, line, wall):
    """Class from the median ring and the shares of the snapshots of type line and wall."""
    if not np.isfinite(ring):
        return "unclear"
    if line >= 0.5:
        return "vortex core" if ring < T_RING else "Bloch line"
    if wall >= 0.5 and ring >= T_RING:
        return "wall or narrow Bloch line"
    return "unclear"


# --- measures on the window arrays of one region -------------------------------------------------------------------

def structure(ms, Ds, part, big, core_r=None):
    """Measures 1 and 2. ms: window m per snapshot (nx, ny, nz, 3); Ds: window grad_dirs per snapshot
    (3, nx, ny, nz, 3); part: cells of the particle; big: the dilated region. Returns (record, geometry for 3)."""
    core_r = CORE_R if core_r is None else core_r
    shp = part.shape
    grid = np.stack(np.meshgrid(*[np.arange(n, dtype=float) for n in shp], indexing="ij"), -1)
    Gs = [H.tensor(D, big) for D in Ds]
    t = smallest(sum(Gs))
    types, turn = [], []
    for G in Gs:
        types.append(H.structure_type(G))
        u = smallest(G)
        turn.append(math.degrees(math.acos(min(1.0, abs(float(u @ t))))))
    e1, e2 = frame(t)
    c0 = grid[big].mean(0)
    sl = np.floor((grid - c0) @ t + 0.5).astype(int)
    pa, pb = grid @ e1, grid @ e2
    slices = [j for j in np.unique(sl[big]) if (big & (sl == j) & part).any()]
    smask = {j: part & (sl == j) for j in slices}
    centre, ring = [], []
    path = {j: np.full((len(ms), 3), np.nan) for j in slices}
    for n, (m, D) in enumerate(zip(ms, Ds)):
        e = (D ** 2).sum(axis=(0, -1))
        for j in slices:
            cand = smask[j] & big
            ic = np.unravel_index(np.argmax(np.where(cand, e, -1.0)), shp)
            d2 = (pa - pa[ic]) ** 2 + (pb - pb[ic]) ** 2
            near = smask[j] & (d2 <= core_r ** 2)
            w = e[near]
            x = (grid[near] * w[:, None]).sum(0) / w.sum() if w.sum() > 0 else grid[ic]
            path[j][n] = x - (x @ t) * t
            centre.append(abs(float(m[ic] @ t)))
            dr = np.sqrt((pa - x @ e1) ** 2 + (pb - x @ e2) ** 2)
            rc = smask[j] & (dr >= RING[0]) & (dr <= RING[1])
            if rc.sum() >= RING_MIN:
                v = m[rc] @ t
                ring.append(0.5 * float(v.max() - v.min()))
    share = lambda k: float(np.mean([ty[0] == k for ty in types]))
    rec = {"snapshots": len(ms), "line": share("line"), "wall": share("wall"), "3D": share("3D"),
           "l2_l1_med": float(np.median([ty[1] for ty in types])), "l3_l1_med": float(np.median([ty[2] for ty in types])),
           "turn_med_deg": float(np.median(turn)), "turn_max_deg": float(np.max(turn)), "slices": len(slices),
           "centre": float(np.median(centre)) if centre else float("nan"),
           "ring": float(np.median(ring)) if ring else float("nan")}
    rec["class"] = classify(rec["ring"], rec["line"], rec["wall"])
    return rec, {"t": t, "sl": sl, "slices": slices, "path": path}


def spectral_derivative(x, omega):
    """d/dt of a periodic series (n, ...) sampled at n equal phases, harmonics 1 ... K only."""
    n = x.shape[0]
    k = np.fft.fftfreq(n, 1.0 / n)
    keep = (np.abs(k) >= 1) & (np.abs(k) <= K)
    c = np.fft.fft(x, axis=0) * (1j * k * omega * keep).reshape((n,) + (1,) * (x.ndim - 1))
    return np.real(np.fft.ifft(c, axis=0))


def harmonic_power(x):
    """Mean square of harmonic k = 1 ... K of a real periodic series (n, cells, 3), summed over cells and
    components: 2 |c_k|^2 with c_k = fft / n."""
    c = np.fft.fft(x, axis=0) / x.shape[0]
    return np.array([2.0 * float((np.abs(c[k]) ** 2).sum()) for k in range(1, K + 1)])


def thiele(ms, Ds, part, geo, per_cycle, omega, coef):
    """Measure 3 for each cycle: rigid and snapshot loss per harmonic (window sums of W/m³), path measures."""
    sl, slices, path = geo["sl"], geo["slices"], geo["path"]
    out = []
    for c in range(len(ms) // per_cycle):
        r = slice(c * per_cycle, (c + 1) * per_cycle)
        vel, steps, amps, vharm, jumps = {}, [], [], [], 0
        for j in slices:
            x = path[j][r]
            if not np.all(np.isfinite(x)):
                continue
            step = float(np.linalg.norm(np.roll(x, -1, 0) - x, axis=1).max())
            if step > JUMP:                                 # the core cell jumped to another structure
                jumps += 1
                continue
            v = spectral_derivative(x, omega)
            vel[j] = v
            steps.append(step)
            amps.append(float(np.ptp(H.principal(x))))
            p = harmonic_power(v[:, None, :])
            vharm.append(p[1:].sum() / p.sum() if p.sum() > 0 else float("nan"))
        if not vel:
            continue
        have = np.array(sorted(vel))
        cells = np.argwhere(part)
        js = sl[part]
        nearest = have[np.abs(js[:, None] - have[None, :]).argmin(1)]
        rig, mea = [], []
        for n in range(r.start, r.stop):
            v = np.stack([vel[j][n - r.start] for j in have])[np.searchsorted(have, nearest)]     # (cells, 3)
            D = Ds[n][:, cells[:, 0], cells[:, 1], cells[:, 2], :]                                  # (3, cells, 3)
            rig.append(-np.einsum("ci,icj->cj", v, D))
            mea.append(ms[n][part])
        rig, mea = np.array(rig), np.array(mea)
        kk = np.arange(1, K + 1)
        out.append({"P_rigid": coef * float((rig ** 2).sum(axis=(1, 2)).mean()),
                    "f_rigid": coef * float(harmonic_power(rig)[0]),
                    "snap_k": coef * (kk * omega) ** 2 * harmonic_power(mea),
                    "step_max_cells": float(np.median(steps)), "amp_cells": float(np.median(amps)),
                    "v_harm_share": float(np.nanmedian(vharm)), "jump_share": jumps / len(slices)})
    return out


# --- run data ----------------------------------------------------------------------------------------------------

def stage_record(d, stage, acc, ids, cfg):
    import domains
    mask = ids >= 0
    N = ids.shape[0]
    alpha, Ms = cfg["alpha"], cfg["Js"] / MU0
    coef = alpha * MU0 * Ms / GAMMA
    z = np.load(acc)
    omega = 2 * math.pi * float(z["freq"])
    p = z["p"].astype(np.float64)
    pf = coef * omega ** 2 * (np.abs(z["m1"].astype(np.complex128)) ** 2).sum(-1) / 2.0
    pk = coef * omega ** 2 * np.arange(1, K + 1) ** 2 * z["mk2"][:, :K].astype(np.float64) / 2.0 if "mk2" in z.files \
        else None
    regs, hot_share, _ = H.regions(p - pf, ids, mask)

    def put(a):
        f = np.zeros(ids.shape + a.shape[1:], a.dtype)
        f[mask] = a
        return f
    P, PF = put(p), put(pf)
    PK = put(pk) if pk is not None else None
    R = []
    for sel, share in regs:
        k = int(np.bincount(ids[sel]).argmax())
        (ix, iy, iz), _ = H.window_index(sel, N, H.WINDOW)
        w = np.ix_(ix, iy, iz)
        part = ids[w] == k
        big = H.dilate_in(sel, ids == k, H.DILATE)[w]
        acc_w = {"p": float(P[w][part].sum()), "pf": float(PF[w][part].sum()),
                 "pk": PK[w][part].sum(0) if PK is not None else None}
        R.append({"share": share, "cells": int(sel.sum()), "particle": k, "w": w, "part": part, "big": big,
                  "acc": acc_w, "ms": [], "Ds": []})
    files = H.snapshot_files(d, stage)
    per_cycle = {len(v) for v in files.values()}
    assert per_cycle == {16}, files.keys()
    for c, fl in files.items():
        for f in fl:
            st = domains.load_state(f)
            assert np.array_equal(st["ids"], ids), f
            m = st["m"].astype(np.float32)
            D = H.grad_dirs(m, ids)
            for r in R:
                r["ms"].append(m[r["w"]])
                r["Ds"].append(D[(slice(None),) + r["w"]])
    out = []
    for r in R:
        rec, geo = structure(r["ms"], r["Ds"], r["part"], r["big"])
        cyc = thiele(r["ms"], r["Ds"], r["part"], geo, 16, omega, coef)
        a = r["acc"]
        sk = np.mean([x["snap_k"] for x in cyc], 0)
        ak = a["pk"] if a["pk"] is not None else np.full(K, np.nan)
        div = lambda x, y: x / y if y > 0 else float("nan")
        P_rig = float(np.mean([x["P_rigid"] for x in cyc]))
        f_rig = float(np.mean([x["f_rigid"] for x in cyc]))
        rec.update({"share": r["share"], "cells": r["cells"], "particle": r["particle"], "cycles": len(cyc),
                    "step_max_cells": float(np.mean([x["step_max_cells"] for x in cyc])),
                    "amp_lex": float(np.mean([x["amp_cells"] for x in cyc])) * cfg["dx_lex"],
                    "v_harm_share": float(np.mean([x["v_harm_share"] for x in cyc])),
                    "jump_share": float(np.mean([x["jump_share"] for x in cyc])),
                    "P_rigid": P_rig, "q_rigid": P_rig - f_rig, "p_acc": a["p"], "q_acc": a["p"] - a["pf"],
                    "snap_f": float(sk[0]), "snap_2K": float(sk[1:].sum()),
                    "acc_f": float(ak[0]), "acc_2K": float(ak[1:].sum())})
        rec["ratio_P"] = div(rec["p_acc"], rec["P_rigid"])
        rec["ratio_q"] = div(rec["q_acc"], rec["q_rigid"])
        rec["snap_over_acc_f"] = div(rec["snap_f"], rec["acc_f"])
        rec["snap_over_acc_2K"] = div(rec["snap_2K"], rec["acc_2K"])
        out.append(rec)
    return {"hot_share": hot_share, "regions": out, "cycles": list(files)}


def run(d):
    import domains
    d = pathlib.Path(d)
    cfg = json.loads((d / "config.json").read_text())
    ids = domains._ids_from_config(cfg).astype(np.int8)
    out = {"run": d.name, "stages": {}}
    for acc in sorted(d.glob("acc_*_*.npz")):
        stage = re.match(r"acc_\d+_(.+)\.npz", acc.name).group(1)
        if H.snapshot_files(d, stage):
            out["stages"][stage] = stage_record(d, stage, acc, ids, cfg)
    return out


def report(res):
    L = ["# Hot-spot follow-up: structure per snapshot, vortex core or Bloch line, Thiele scale (eval_hotcores.py)", "",
         "Main regions of the loss outside f as in eval_hotspots.py. line / wall / 3D: shares of the snapshots of each "
         "type; turn: angle of the line direction of a snapshot to the mean one (median / max). centre: |m.t| of the "
         "core cell (information); ring: (max - min) / 2 of m.t at %g ... %g cells from the core (medians); class "
         "from the ring (threshold %.2f) and the type shares. Thiele: loss of the window from the accumulators over "
         "the loss of the rigid motion along the core path, P (all) and q (P minus the loss at f); snap / acc: loss "
         "at f and at 2f ... %df from the 16 snapshot phases over the same from the accumulators; v harm: share of "
         "k >= 2 in the mean square core velocity; step: largest core step between 2 phases (cells); amp: core path, "
         "l_ex." % (RING[0], RING[1], T_RING, K), "",
         "| run | stage | share | particle | line | wall | 3D | turn, deg | centre | ring | class | amp, l_ex | step | "
         "jumps | v harm | P acc / rigid | q acc / rigid | f snap / acc | 2f ... 7f snap / acc |", "|---" * 19 + "|"]
    for r in res:
        for stage, s in r["stages"].items():
            for g in s["regions"]:
                L.append("| %s | %s | %.3f | %d | %.2f | %.2f | %.2f | %.0f / %.0f | %.2f | %.2f | %s | %.1f | %.1f | "
                         "%.2f | %.2f | %.2f | %.2f | %.2f | %.2f |" % (
                             r["run"], stage, g["share"], g["particle"], g["line"], g["wall"], g["3D"],
                             g["turn_med_deg"], g["turn_max_deg"], g["centre"], g["ring"], g["class"], g["amp_lex"],
                             g["step_max_cells"], g["jump_share"], g["v_harm_share"], g["ratio_P"], g["ratio_q"],
                             g["snap_over_acc_f"], g["snap_over_acc_2K"]))
    return "\n".join(L)


# --- synthetic states ----------------------------------------------------------------------------------------------

AXIS = np.array([0.31, 0.52, 0.80]) / np.linalg.norm([0.31, 0.52, 0.80])


def synthetic_line(kind, dx=3.0, d=300.0, core=1.0, Leff=18.0, ripple_deg=0.0, seed=1, shift=(0.0, 0.0)):
    """One sphere (diameter d, l_ex) with one line structure along an oblique axis, as domains.synthetic_state.
    Local frame (e1, e2, e3 = AXIS); shift moves the structure in the (e1, e2) plane (l_ex).
      vortex  : m curls around e3, core profile m.e3 = exp(-(rho / core)^2)
      bloch_v : Bloch line along e3 in a 180 deg wall (normal e1, domains +-e3, wall width Leff): the wall centre
                turns from -e2 to +e2 through e1 (width core); the line is along the domain magnetization
      bloch_h : the same wall, Bloch line along e2 (across the domain magnetization), wall centre turns along e3"""
    import domains
    st, _ = domains.synthetic_state("ripple", dx, Leff=Leff, d=d, ripple_deg=ripple_deg, seed=seed)
    N = st["m"].shape[0]
    x = (np.arange(N) + 0.5) * dx - N * dx / 2
    X = np.stack(np.meshgrid(x, x, x, indexing="ij"), -1)
    e3 = AXIS
    e1, e2 = frame(e3)
    a, b, c = X @ e1 - shift[0], X @ e2 - shift[1], X @ e3
    if kind == "vortex":
        rho = np.sqrt(a * a + b * b)
        mz = np.exp(-(rho / core) ** 2)
        s = np.sqrt(1 - mz ** 2) / np.maximum(rho, 1e-9)
        m = (-b * s)[..., None] * e1 + (a * s)[..., None] * e2 + mz[..., None] * e3
    elif kind in ("bloch_v", "bloch_h"):
        th = 2.0 * np.arctan(np.exp(a / Leff))
        psi = 2.0 * np.arctan(np.exp((b if kind == "bloch_v" else c) / core)) - math.pi / 2
        w = np.cos(psi)[..., None] * e1 + np.sin(psi)[..., None] * e2
        m = np.cos(th)[..., None] * e3 + np.sin(th)[..., None] * w
    else:
        raise ValueError(kind)
    if ripple_deg > 0:
        m = m + (st["m"] - np.array([1.0, 0.0, 0.0]))
    m /= np.linalg.norm(m, axis=-1, keepdims=True)
    m = np.where((st["ids"] >= 0)[..., None], m, 0.0).astype(np.float32)
    t = e2 if kind == "bloch_h" else e3
    return m, st["ids"], t


def synthetic_region(m, ids):
    """Window arrays of the main region of |grad m|^2 (as the test of eval_hotspots): (w, part, big)."""
    D = H.grad_dirs(m, ids)
    mask = ids >= 0
    regs, _, _ = H.regions((D ** 2).sum(axis=(0, -1))[mask], ids, mask)
    sel = regs[0][0]
    (ix, iy, iz), _ = H.window_index(sel, ids.shape[0], H.WINDOW)
    w = np.ix_(ix, iy, iz)
    return w, ids[w] == 0, H.dilate_in(sel, ids == 0, H.DILATE)[w]


def synthetic_structure(kind, **kw):
    m, ids, t = synthetic_line(kind, **kw)
    w, part, big = synthetic_region(m, ids)
    D = H.grad_dirs(m, ids)
    rec, geo = structure([m[w]], [D[(slice(None),) + w]], part, big)
    rec["t_err_deg"] = math.degrees(math.acos(min(1.0, abs(float(geo["t"] @ t)))))
    return rec


def synthetic_thiele(core=6.0, amp=12.0, amp3=0.0, dense=32, eps=1e-3, core_r=(CORE_R,), d=300.0):
    """Vortex moved rigidly across its axis: x(phase) = amp sin(phase) + amp3 sin(3 phase) along e1 (l_ex). Returns,
    for each core_r, the measure-3 record of the method (16 phases) and the true values as window sums with coef = 1
    and w = 2 pi: P (mean of |dm/dt|^2 over the 16 phases, dm/dt by central differences in the phase, eps) and the
    loss at f (harmonic 1 of `dense` phases)."""
    x = lambda f: amp * math.sin(f) + amp3 * math.sin(3 * f)
    st = lambda f: synthetic_line("vortex", d=d, core=core, shift=(x(f), 0.0))[0]
    ids = synthetic_line("vortex", d=d, core=core)[1]
    phases = [2 * math.pi * (n + 1) / 16 for n in range(16)]
    snaps = [st(f) for f in phases]
    hot = np.zeros(ids.shape, bool)
    for m in snaps:                                         # region: the cells the core passes (as the hot cells)
        e = (H.grad_dirs(m, ids) ** 2).sum(axis=(0, -1))
        hot |= (e >= np.quantile(e[ids >= 0], 0.995)) & (ids >= 0)
    (ix, iy, iz), _ = H.window_index(hot, ids.shape[0], H.WINDOW)
    w = np.ix_(ix, iy, iz)
    part, big = ids[w] == 0, H.dilate_in(hot, ids == 0, H.DILATE)[w]
    ms = [m[w] for m in snaps]
    Ds = [H.grad_dirs(m, ids)[(slice(None),) + w] for m in snaps]
    P_true = float(np.mean([(((st(f + eps)[w][part] - st(f - eps)[w][part]) / (2 * eps) * 2 * math.pi) ** 2).sum()
                            for f in phases]))
    c1 = np.fft.fft(np.array([st(2 * math.pi * n / dense)[w][part] for n in range(dense)]), axis=0)[1] / dense
    f_true = 2.0 * (2 * math.pi) ** 2 * float((np.abs(c1) ** 2).sum())
    out = []
    for cr in core_r:
        rec, geo = structure(ms, Ds, part, big, core_r=cr)
        out.append((cr, rec, thiele(ms, Ds, part, geo, 16, 2 * math.pi, 1.0)[0]))
    return out, {"P": P_true, "f": f_true}


def selftest():
    print("| kind | core, l_ex | ripple, deg | line | t error, deg | centre | ring |")
    print("|---|---|---|---|---|---|---|")
    for kind in ("vortex", "bloch_v", "bloch_h"):
        for core in (1.0, 3.0, 9.0):
            for rp in (0.0, 20.0):
                r = synthetic_structure(kind, core=core, ripple_deg=rp)
                print("| %s | %g | %g | %.2f | %.1f | %.2f | %.2f |" % (kind, core, rp, r["line"], r["t_err_deg"],
                                                                        r["centre"], r["ring"]), flush=True)
    print()
    print("| core, l_ex | amp, l_ex | amp3, l_ex | core_r, cells | true / rigid P | true / rigid q | snap / true f | "
          "true q / P | v harm | jump share | amp found, cells |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for core, amp, amp3 in ((3.0, 12.0, 0.0), (3.0, 12.0, 3.0), (3.0, 30.0, 0.0), (9.0, 12.0, 0.0)):
        res, tr = synthetic_thiele(core=core, amp=amp, amp3=amp3, core_r=(1.5, 2.5))
        for cr, rec, cyc in res:
            print("| %g | %g | %g | %g | %.2f | %.2f | %.2f | %.2f | %.3f | %.2f | %.1f |" % (
                core, amp, amp3, cr, tr["P"] / cyc["P_rigid"], (tr["P"] - tr["f"]) / (cyc["P_rigid"] - cyc["f_rigid"]),
                cyc["snap_k"][0] / tr["f"], (tr["P"] - tr["f"]) / tr["P"], cyc["v_harm_share"], cyc["jump_share"],
                cyc["amp_cells"]), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=None)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("runs", nargs="*")
    a = ap.parse_args()
    if a.selftest:
        selftest()
        return
    res = [run(d) for d in a.runs]
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(res, indent=1))
    print(report(res))


if __name__ == "__main__":
    main()
