#!/usr/bin/env python
"""
Wall measures of magnetization states (PLAN.md, "Domain analysis").

Question: does beta depend on the domain structure? This tool measures two
quantities of a state m(r) without a domain segmentation (a segmentation by
direction clusters was tested and archived: tag archive/domains-segmentation-2026-10-03).

  1. Texture T(w) = (1/pi) integral |grad m_w| dV  [l_ex^2]
     m_w = m smoothed by a Gaussian of width w (normalized convolution inside
     the particles, then |m_w| = 1). Across a wall m turns by the wall angle,
     independent of the wall width: a 180 deg wall gives its area, a 90 deg wall
     half its area. The ripple (scale L_eff) is suppressed by the smoothing.
     Wall measure: A_w = T(2 L_eff) / CAL (+-12 %). CAL = T(2 L_eff) / true area on
     synthetic states with 20 deg ripple (a real state, d/l_ex 96: 13 deg):
         wall spacing / curvature radius about 5 ... 8 L_eff (d = 300):  0.89 ... 0.94
         about 4 ... 6 L_eff (d = 212):                                 0.68 ... 0.88
         about 3 ... 4 L_eff (d = 150):                                 0.44 ... 0.76
     wall width 0.5 / 1 / 1.5 / 2 x pi L_eff (d = 300):     0.92 ... 0.95 / 0.89 ... 0.94 /
                                                          0.83 ... 0.91 / 0.75 ... 0.88
     dx = 3 and 1.5 agree to 0.3 %. Thus A_w is a wall area (+-10 %) only if
     the walls are >= 5 L_eff apart and not curved more strongly; else T(2 L_eff)
     is the rotation of m on scales above 2 L_eff, not a wall area. A smaller
     width does not help: at w = 1 L_eff the ripple adds 13 ... 50 % of A.
     Limit: T counts every rotation on scales above L_eff; it does not separate
     walls from a continuous rotation (vortex, curling). The localization
     (int g)^2 / (V int g^2) does not separate them either (selftest: walls
     0.61 ... 0.89, helix 0.98, vortex 0.48), because the walls (width pi L_eff)
     fill a large part of the particle.
  2. Switched volume between two states a, b of the same run:
     V_sw = integral arccos(m_a . m_b) / pi dV  with m smoothed over w_sw = 1 L_eff
     [l_ex^3]; a reversed region (180 deg) counts fully, a wall moved by s gives A s
     (linear in s; the form (1 - m_a . m_b)/2 grows only with s^2).
     V_sw counts every change of m: a uniform rotation by phi inside the domains
     gives V phi/180 deg (10 deg: 0.055 V; 3 walls moved by 0.5 L_eff: 0.077 V).
     It does not separate wall motion from rotation.
     Noise floor: a ripple pattern that changes completely between a and b
     gives V_sw = 0.06 V (selftest, 20 deg); in the model the ripple is tied to
     the fixed cubes.
  3. Correlation length l_C (1/e length of the autocorrelation of m - <m> inside each
     particle), the analogue of the spin-misalignment correlation length of magnetic SANS
     (Michels: l_C = L + l_H(H)). Synthetic: ripple 15 / 48 / 78 l_ex for a correlation of
     1 / 4 / 8 L_eff, 180 deg walls 48, bubble 63, vortex 102, helix 96 l_ex.
  4. With the 4 phases of one cycle: half-cycle switched volume and the mean wall path
     s = V_sw / A_w (Bertotti: active wall area and its path set the excess loss).
     Synthetic: 3 walls moved by 6 l_ex give s = 5.4 l_ex (-9 % from CAL).
  5. Rotation test F90: share of windows (length 1.5 pi L_eff along random lines, inside one
     particle, m smoothed over 0.5 L_eff) in which m turns by more than 90 deg. A 180 deg wall
     turns m by up to ~170 deg over this length; ripple and continuous rotation turn it less
     if the particle is much larger than the window. Synthetic F90 [%]:
                    vortex  helix  one 180 wall  three 180 walls
         d = 300     4.3     0.0      12.3          37
         d = 212     9.4     6.9      18.4          54
         d <= 150   22 ... 65 (no separation)
     Thus at d = 300: F90 well above ~5 % means 180 deg walls; at d = 212 the margin is only ~2x;
     below d ~ 200 F90 does not separate. 90 deg walls are not detected (F90 = 0.9 %); in the model
     the anisotropy is uniaxial, so 180 deg walls are the main type. dx 3 and 1.5 agree.
     Stereology: a random segment of length W hits a wall with probability W S_V / 2, thus
     A_F90 = 2 F90 V / W is a second estimate of the 180 deg wall area (synthetic, d = 300:
     three walls 1.85 vs 1.86e5, one wall -13 %, bubble +32 %, vortex 0.21e5 = false area of
     ~0.3 wall discs). A_w >> A_F90 means that a part of T is continuous rotation.
Long-wave ripple (correlation 2 ... 8 L_eff, as at small internal fields): with walls,
T(2 L_eff)/A stays 0.91 ... 0.94; without walls the ripple alone gives T(2 L_eff) =
0.04 ... 0.17 of the area of 3 walls (d = 300), a floor of up to ~0.4 wall discs.
No Porod q^-4 range exists for these walls (width pi L_eff, spacing ~5 L_eff): the
structure factor does not separate walls from a continuous rotation either.
Smoothing runs inside each particle separately: the kernel is wider than the
gaps between the particles (test: 4 uniform particles give T = 0).

Lengths in l_ex, areas in l_ex^2, volumes in l_ex^3.

    python domains.py runs/V_s901/m_amp01_c6_ph*.vti   # states of one run (in phase order)
    python domains.py runs/B300_s921/init.pt
    python domains.py --selftest                       # synthetic states with known walls
"""
import argparse
import json
import math
import pathlib
import sys

import numpy as np
from scipy import ndimage

WIDTHS = (0.0, 0.5, 1.0, 2.0, 3.0, 4.0)      # smoothing widths / L_eff
W_WALL = 2.0                                  # width of the wall measure / L_eff
W_SW = 1.0                                    # width of the switched volume / L_eff
W_ROT = 1.5 * math.pi                         # window of the rotation test / L_eff (1.5 wall widths)
CAL = 0.85                                    # T(2 L_eff) / true wall area (selftest: 0.75 ... 0.95, see below)


# --- input -----------------------------------------------------------------

def load_state(path, coarsen=1):
    """m [N,N,N,3], ids [N,N,N] (particle index, -1 = void), mask = ids >= 0, dx, Leff (l_ex)."""
    path = pathlib.Path(path)
    run = path.parent
    cfg = json.loads((run / "config.json").read_text())
    dx, Leff, N = cfg["dx_lex"], cfg["Leff_lex"], cfg["N_used"]
    if path.suffix == ".vti":
        import pyvista as pv
        m = np.asarray(pv.read(str(path)).cell_data["m"], dtype=np.float32).reshape((N, N, N, 3), order="F")
        g = run / "geometry.vti"
        if g.exists():
            ids = np.rint(np.asarray(pv.read(str(g)).cell_data["particle"])).reshape((N, N, N), order="F")
        else:
            ids = _ids_from_config(cfg)
    else:
        import torch
        m = torch.load(str(path), map_location="cpu")["m"].float().numpy()
        ids = _ids_from_config(cfg)
    ids = ids.astype(np.int8)
    if coarsen > 1:
        m, ids = coarsen_state(m, ids, coarsen)
        dx *= coarsen
    return make_state(m, ids, dx, Leff, str(path))


def make_state(m, ids, dx, Leff, name):
    return {"m": m, "ids": ids, "mask": ids >= 0, "dx": dx, "Leff": Leff, "name": name}


def _ids_from_config(cfg):
    from geometry import fcc_box
    return fcc_box(cfg["N_used"], cfg["d"], cfg["dx"])["ids"]


def coarsen_state(m, ids, c):
    """Block average by c (for example dx 1.5 -> 3: the same analysis grid on both meshes).
    A coarse cell belongs to the particle that fills at least half of it."""
    N = m.shape[0] // c
    b = lambda a: a[:N * c, :N * c, :N * c].reshape(N, c, N, c, N, c, *a.shape[3:])
    mm = b(m * (ids >= 0)[..., None]).sum(axis=(1, 3, 5))
    mm /= np.maximum(np.linalg.norm(mm, axis=-1, keepdims=True), 1e-12)
    idc = np.full((N, N, N), -1, dtype=np.int8)
    for p in np.unique(ids[ids >= 0]):
        idc[b((ids == p).astype(np.float32)).mean(axis=(1, 3, 5)) >= 0.5] = p
    return np.where((idc >= 0)[..., None], mm, 0.0).astype(np.float32), idc


# --- field operations (periodic box) ------------------------------------------

def grad_norm(f, ids, dx):
    """|grad f| per cell from differences between face neighbours in the same particle (mean of
    the forward and the backward difference that exist). f: [N,N,N,3]."""
    mask = ids >= 0
    g2 = np.zeros(mask.shape, dtype=np.float64)
    for ax in range(3):
        acc = np.zeros(mask.shape, dtype=np.float64)
        cnt = np.zeros(mask.shape, dtype=np.float64)
        for s in (-1, 1):
            ok = mask & (np.roll(ids, s, axis=ax) == ids)
            d = np.roll(f, s, axis=ax) - f
            acc += np.where(ok, (d * d).sum(-1), 0.0)
            cnt += ok
        g2 += acc / np.maximum(cnt, 1.0)
    return np.where(mask, np.sqrt(g2) / dx, 0.0)


def smooth_m(m, ids, sigma_cells):
    """Gaussian smoothing inside EACH particle separately (normalized convolution; the kernel is
    wider than the gaps, so m of a neighbour particle must not enter), then |m| = 1."""
    if sigma_cells <= 0:
        return m
    out = np.zeros_like(m, dtype=np.float32)
    for p in np.unique(ids[ids >= 0]):
        w = (ids == p).astype(np.float32)
        s = np.stack([ndimage.gaussian_filter(m[..., i] * w, sigma_cells, mode="wrap") for i in range(3)], -1)
        s /= np.maximum(np.linalg.norm(s, axis=-1, keepdims=True), 1e-12)
        out[ids == p] = s[ids == p]
    return out


def measures(st, widths=WIDTHS):
    widths = tuple(sorted(set(widths) | {W_WALL}))
    """Texture T(w), wall area A_w, localization P, ripple estimate of one state."""
    dx, L, mask, ids = st["dx"], st["Leff"], st["mask"], st["ids"]
    V = float(mask.sum() * dx ** 3)
    T = {}
    for w in widths:
        g = grad_norm(smooth_m(st["m"], ids, w * L / dx), ids, dx)
        T[w] = float(g.sum() * dx ** 3 / math.pi)
    ms = smooth_m(st["m"], ids, W_WALL * L / dx)
    ang = np.degrees(np.arccos(np.clip((st["m"] * ms).sum(-1)[mask], -1, 1)))
    rot = rotation_angles(st, window_leff=W_ROT)
    return {"name": st["name"], "V": V, "T": T, "A_w": T[W_WALL] / CAL,
            "ripple_rms_deg": float(np.sqrt(np.mean(ang ** 2))), "l_C": correlation_length(st),
            "F90": float((rot > 90.0).mean()), "F60": float((rot > 60.0).mean()),
            "A_F90": float(2.0 * (rot > 90.0).mean() * V / (W_ROT * L))}


def correlation_length(st):
    """l_C [l_ex]: distance at which the autocorrelation of dm = m - <m>_particle (inside each
    particle, normalized by the overlap of the particle with itself) falls to 1/e. The analogue of
    the correlation length of the spin misalignment in magnetic SANS (Michels)."""
    ids, dx = st["ids"], st["dx"]
    N = ids.shape[0]
    num = np.zeros((N, N, N))
    den = np.zeros((N, N, N))
    for p in np.unique(ids[ids >= 0]):
        w = (ids == p).astype(np.float64)
        dm = (st["m"] - st["m"][ids == p].mean(0)) * w[..., None]
        for i in range(3):
            F = np.fft.rfftn(dm[..., i])
            num += np.fft.irfftn(F * np.conj(F), s=(N, N, N), axes=(0, 1, 2))
        Fw = np.fft.rfftn(w)
        den += np.fft.irfftn(Fw * np.conj(Fw), s=(N, N, N), axes=(0, 1, 2))
    k = np.minimum(np.arange(N), N - np.arange(N)) * dx
    r = np.sqrt(k[:, None, None] ** 2 + k[None, :, None] ** 2 + k[None, None, :] ** 2).ravel()
    ok = den.ravel() > 0.2 * den.flat[0]                   # enough overlap for a stable estimate
    rb = np.arange(0.0, r[ok].max(), dx)
    idx = np.digitize(r[ok], rb)
    C = np.bincount(idx, weights=num.ravel()[ok]) / np.maximum(np.bincount(idx, weights=den.ravel()[ok]), 1e-30)
    C = C[1:len(rb)] / C[1]
    below = np.nonzero(C < math.exp(-1.0))[0]
    return float(rb[below[0]]) if len(below) else float("nan")


def rotation_angles(st, window_leff=math.pi, smooth_leff=0.5, n_lines=3000, seed=0):
    """Net rotation angle of m over a window of length window_leff * L_eff along random lines
    (the window lies inside one particle). m is smoothed over smooth_leff * L_eff first.
    A wall turns m by its wall angle over about one wall width (pi L_eff); ripple, a helix or a
    vortex away from its core turn m much less over the same length. Returns angles in deg."""
    ids, dx, L = st["ids"], st["dx"], st["Leff"]
    ms = smooth_m(st["m"], ids, smooth_leff * L / dx)
    N = ids.shape[0]
    rng = np.random.default_rng(seed)
    step = 0.5
    nw = int(round(window_leff * L / dx / step))
    n_s = 4 * N
    out = []
    for _ in range(n_lines):
        v = rng.normal(size=3)
        v /= np.linalg.norm(v)
        pts = rng.uniform(0, N, size=3)[None, :] + np.arange(n_s)[:, None] * step * v[None, :]
        idx = np.floor(pts).astype(int) % N
        pid = ids[idx[:, 0], idx[:, 1], idx[:, 2]]
        mm = ms[idx[:, 0], idx[:, 1], idx[:, 2]]
        same = pid[:-nw] >= 0
        # the window must stay inside one particle: no void sample between both ends
        bad = np.convolve((pid < 0).astype(int), np.ones(nw + 1, dtype=int), "valid") > 0
        ok = same & ~bad & (pid[:-nw] == pid[nw:])
        c = np.clip((mm[:-nw] * mm[nw:]).sum(-1), -1.0, 1.0)
        out.append(np.degrees(np.arccos(c[ok])))
    return np.concatenate(out)


def switched_volume(sa, sb, w=W_SW):
    ma = smooth_m(sa["m"], sa["ids"], w * sa["Leff"] / sa["dx"])
    mb = smooth_m(sb["m"], sb["ids"], w * sb["Leff"] / sb["dx"])
    ang = np.arccos(np.clip((ma * mb).sum(-1), -1.0, 1.0))
    return float((ang / math.pi)[sa["mask"]].sum() * sa["dx"] ** 3)


def analyze(states):
    """With the 4 phases of one cycle (H max, H = 0 falling, H min, H = 0 rising): half-cycle
    switched volume (H max -> H min and back) and the mean wall path s = V_sw / A_w (Bertotti:
    the active wall area and its path set the excess loss)."""
    res = {"states": [measures(st) for st in states]}
    n = len(states)
    res["switched"] = [switched_volume(states[i], states[(i + 1) % n]) for i in range(n)] if n > 1 else []
    if n == 4:
        r = res["states"]
        half = [switched_volume(states[0], states[2]), switched_volume(states[1], states[3])]
        A = 0.5 * (r[0]["A_w"] + r[2]["A_w"])
        res["half_cycle"] = {"V_sw": half, "A_w_mean": A,
                             "wall_path": [v / A if A > 0 else float("nan") for v in half]}
    return res


def report(res):
    for r in res["states"]:
        print("%s\n   V %.4g l_ex^3   A_w %.4g l_ex^2 (A_w/V %.3g 1/l_ex)   ripple %.1f deg\n"
              "   F90 %.1f %%   A_F90 %.4g l_ex^2   A_w/A_F90 %.2f"
              % (r["name"], r["V"], r["A_w"], r["A_w"] / r["V"], r["ripple_rms_deg"], 100 * r["F90"],
                 r["A_F90"], r["A_w"] / r["A_F90"] if r["A_F90"] > 0 else float("inf")))
        print("   T(w) [l_ex^2]: " + "  ".join("w=%g: %.4g" % (w, t) for w, t in r["T"].items())
              + "   l_C %.1f l_ex" % r["l_C"])
    if res["switched"]:
        V = res["states"][0]["V"]
        print("switched volume to the next state (cyclic): " +
              "  ".join("%.4g (%.1f %%)" % (v, 100 * v / V) for v in res["switched"]))
    if "half_cycle" in res:
        h = res["half_cycle"]
        print("half cycle (H max <-> H min, H = 0 falling <-> rising): V_sw %s l_ex^3, mean wall path "
              "s = V_sw/A_w = %s l_ex" % (" / ".join("%.4g" % v for v in h["V_sw"]),
                                         " / ".join("%.2f" % x for x in h["wall_path"])))


# --- synthetic states with known walls ------------------------------------------

def synthetic_state(kind, dx, Leff=12.0, d=300.0, ripple_deg=0.0, seed=0, shift=0.0,
                    normal=(0.31, 0.52, 0.80), wall_factor=1.0, ripple_corr=1.0):
    """One sphere (diameter d) in a box. Walls have the profile theta = 2 atan(exp(s / L_eff))
    (width pi L_eff; wall_factor scales it). The ripple is a random field on cubes of edge L_eff,
    smoothed by a Gaussian of width ripple_corr L_eff / 2 (ripple_corr > 1: long-wave ripple, as at
    small internal fields, l_C = L + l_H). shift moves all walls (switched-volume test). Returns (state, truth) with
    truth = {"A": wall area weighted by angle/180 deg (= expected T), "A_geo": wall area,
    "V_rev": reversed volume relative to shift = 0}.
      ripple  : uniform m, no wall           slab180: 3 parallel 180 deg walls (oblique)
      slab90  : 3 parallel 90 deg walls      bubble : spherical inner domain, radius 0.3 d
      helix   : m turns continuously (one turn over d) - no wall, rotation everywhere
      vortex  : m curls around an axis (core radius L_eff) - no wall"""
    a = d + 24.0
    N = int(round(a / dx))
    x = (np.arange(N) + 0.5) * dx - a / 2
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R = d / 2
    rr = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    mask = rr <= R
    n = np.array(normal) / np.linalg.norm(normal)
    u = np.array([1.0, 0.0, 0.0])
    v = np.cross(n, u)
    v /= np.linalg.norm(v)
    prof = lambda s: 2.0 * np.arctan(np.exp(s / (wall_factor * Leff)))
    cap = lambda p: math.pi * (R * R * p - p ** 3 / 3.0)          # sphere volume between 0 and p (|p| <= R)
    truth = {"A": 0.0, "A_geo": 0.0, "V_rev": 0.0}
    m = None
    if kind in ("ripple", "helix", "vortex"):
        if kind == "ripple":
            m = np.broadcast_to(u, X.shape + (3,)).copy()
        elif kind == "helix":
            ph = 2 * math.pi * Z / d
            m = np.stack([np.cos(ph), np.sin(ph), np.zeros_like(ph)], -1)
        else:
            rho = np.sqrt(X ** 2 + Y ** 2)
            mz = np.exp(-(rho / Leff) ** 2)
            t = np.sqrt(1 - mz ** 2) / np.maximum(rho, 1e-9)
            m = np.stack([-Y * t, X * t, mz], -1)
    elif kind in ("slab180", "slab90", "single180"):
        s = X * n[0] + Y * n[1] + Z * n[2]
        pos = (0.0,) if kind == "single180" else (-0.45 * R, 0.0, 0.4 * R)
        f = 0.5 if kind == "slab90" else 1.0
        th = sum(f * prof(s - p - shift) for p in pos)
        truth["A_geo"] = sum(math.pi * (R ** 2 - (p + shift) ** 2) for p in pos)
        truth["A"] = f * truth["A_geo"]
        if kind != "slab90":
            truth["V_rev"] = sum(abs(cap(p + shift) - cap(p)) for p in pos)
    elif kind == "bubble":
        r0 = 0.3 * d + shift
        th = prof(r0 - rr)
        truth["A_geo"] = truth["A"] = 4 * math.pi * r0 ** 2
        truth["V_rev"] = 4.0 / 3.0 * math.pi * abs(r0 ** 3 - (0.3 * d) ** 3)
    else:
        raise ValueError(kind)
    if m is None:
        m = np.cos(th)[..., None] * u + np.sin(th)[..., None] * v
    if ripple_deg > 0:                                                # random cubes of edge L_eff
        rng = np.random.default_rng(seed)
        nc = int(math.ceil(a / Leff))
        ic = np.minimum((np.arange(N) * dx / Leff).astype(int), nc - 1)
        g = rng.normal(size=(nc, nc, nc, 3))[np.ix_(ic, ic, ic)]
        g = np.stack([ndimage.gaussian_filter(g[..., i], 0.5 * ripple_corr * Leff / dx, mode="wrap")
                      for i in range(3)], -1)
        g -= (g * m).sum(-1, keepdims=True) * m
        g /= np.sqrt(np.mean((g ** 2).sum(-1)[mask]))
        m = m + math.tan(math.radians(ripple_deg)) * g
    m /= np.linalg.norm(m, axis=-1, keepdims=True)
    m = np.where(mask[..., None], m, 0.0).astype(np.float32)
    return make_state(m, np.where(mask, 0, -1).astype(np.int8), dx, Leff,
                      "%s dx=%g ripple=%g shift=%g" % (kind, dx, ripple_deg, shift)), truth


def selftest(d=300.0, dxs=(3.0, 1.5), ripples=(0.0, 20.0),
             kinds=("ripple", "slab180", "slab90", "bubble", "helix", "vortex"), shift=6.0):
    """Wall measure, localization and switched volume against the known truth. The fine mesh is
    analysed on the dx = 3 grid (coarsen 2), as for dx = 1.5 runs."""
    rows = []
    for kind in kinds:
        for dx in dxs:
            for rp in ripples:
                states = []
                for sh in (0.0, shift):
                    st, tr = synthetic_state(kind, dx, d=d, ripple_deg=rp, shift=sh)
                    if dx < 3.0:
                        c = int(round(3.0 / dx))
                        m, ic = coarsen_state(st["m"], st["ids"], c)
                        st = make_state(m, ic, dx * c, st["Leff"], st["name"])
                    states.append((st, tr))
                r = measures(states[0][0])
                tr0, tr1 = states[0][1], states[1][1]
                vsw = switched_volume(states[0][0], states[1][0]) if kind in ("slab180", "bubble") else float("nan")
                row = {"kind": kind, "dx": dx, "ripple": rp, "A_true": tr0["A"], "T2": r["T"][W_WALL],
                       "T2_over_true": r["T"][W_WALL] / tr0["A"] if tr0["A"] else float("nan"),
                       "ripple_est": r["ripple_rms_deg"],
                       "Vsw_true": tr1["V_rev"], "Vsw": vsw,
                       "Vsw_over_true": vsw / tr1["V_rev"] if tr1["V_rev"] else float("nan"),
                       "V": r["V"], **{"T%g" % w: t for w, t in r["T"].items()}}
                rows.append(row)
                print("%-8s dx %-4g ripple %4g | T(2)/A %6.3f  T(2)/V %.3g  ripple_est %5.1f | V_sw/true %6.3f"
                      % (kind, dx, rp, row["T2_over_true"], row["T2"] / row["V"], r["ripple_rms_deg"],
                         row["Vsw_over_true"]), flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("states", nargs="*", help="m_*.vti or init.pt / checkpoint.pt of ONE run")
    ap.add_argument("--coarsen", type=int, default=1, help="block average (2 for dx 1.5 -> analysis grid dx 3)")
    ap.add_argument("--json", default=None, help="write the result as json")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--selftest_csv", default=None)
    a = ap.parse_args()
    if a.selftest:
        rows = selftest()
        if a.selftest_csv:
            import csv
            with open(a.selftest_csv, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
                w.writeheader()
                w.writerows(rows)
        return
    if not a.states:
        sys.exit("no states given")
    res = analyze([load_state(p, a.coarsen) for p in a.states])
    report(res)
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
