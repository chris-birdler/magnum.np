#!/usr/bin/env python
"""
Domain analysis of magnetization states (PLAN.md, "Domain analysis").

Question: does beta depend on the domain structure (wall area, switched
volume)? This tool finds domains and walls in a state m(r) of the powder model
and measures them. It does not use a fixed angle threshold.

Method
  1. Direction distribution: m (smoothed over sigma_fit = 1 L_eff to suppress
     the ripple) is fitted with a mixture of von Mises-Fisher distributions on
     the unit sphere, K = 1 ... K_max; the BIC selects K. If K = 1, there are no
     separate domain directions.
  2. Segmentation: hidden Markov random field (mean field). Each cell gets a
     probability q_k for each direction cluster k. Data term: vMF likelihood
     of the raw m. Prior: Potts coupling lam between face neighbours in the
     same particle (a boundary costs lam per neighbour pair).
  3. Domains: connected regions (periodic box) of one label. A region is a
     domain only if it is thicker than a wall (largest inscribed radius >=
     r_min = 2 L_eff; the wall width is pi L_eff). Thinner regions (for
     example the centre of a 180 deg wall, where m points sideways) are not
     allowed for that label, and the field is segmented again.
  4. Wall area, two independent estimates:
     (i)  co-area: A = 1/2 sum_k integral |grad q_k| dV (each wall counts once,
          from both sides; independent of the wall width);
     (ii) stereology: random lines, S_V = 2 P_L (P_L = label changes per
          line length inside the particles), A = S_V V.
  5. Texture (needs no domains): T(w) = (1/pi) integral |grad m_w| dV with m
     smoothed by a Gaussian of width w (normalized convolution inside the
     particles). A 180 deg wall gives its area, a 90 deg wall half its area.
  6. Switched volume between two states with the same clusters:
     V_sw = integral 1/2 sum_k |q_k(a) - q_k(b)| dV.

Lengths are in l_ex, areas in l_ex^2, volumes in l_ex^3.

    python domains.py runs/V_s901/m_amp00_c6_ph090.vti ...   # states of one run
    python domains.py runs/B300_s921/init.pt
    python domains.py --selftest                              # synthetic states with known walls
"""
import argparse
import json
import math
import pathlib
import sys

import numpy as np
from scipy import ndimage

HERE = pathlib.Path(__file__).resolve().parent
LAMS = (0.25, 0.5, 1.0, 2.0)                 # Potts couplings of the stability scan
WIDTHS = (0.0, 0.5, 1.0, 2.0, 3.0, 4.0)      # texture smoothing widths / L_eff


# --- input -----------------------------------------------------------------

def load_state(path, coarsen=1):
    """m [N,N,N,3], mask [N,N,N] (bool), ids [N,N,N] (-1 void), dx_lex, Leff_lex."""
    path = pathlib.Path(path)
    run = path.parent
    cfg = json.loads((run / "config.json").read_text())
    dx_lex, Leff = cfg["dx_lex"], cfg["Leff_lex"]
    N = cfg["N_used"]
    if path.suffix == ".vti":
        import pyvista as pv
        m = np.asarray(pv.read(str(path)).cell_data["m"], dtype=np.float32).reshape((N, N, N, 3), order="F")
        g = run / "geometry.vti"
        if g.exists():
            ids = np.asarray(pv.read(str(g)).cell_data["particle"]).reshape((N, N, N), order="F")
            ids = np.rint(ids).astype(np.int8)
        else:
            ids = _ids_from_config(cfg)
    else:
        import torch
        m = torch.load(str(path), map_location="cpu")["m"].float().numpy()
        ids = _ids_from_config(cfg)
    mask = ids >= 0
    if coarsen > 1:
        m, mask, ids = coarsen_state(m, mask, ids, coarsen)
        dx_lex *= coarsen
    return {"m": m, "mask": mask, "ids": ids, "dx": dx_lex, "Leff": Leff, "name": str(path)}


def _ids_from_config(cfg):
    from geometry import fcc_box
    return fcc_box(cfg["N_used"], cfg["d"], cfg["dx"])["ids"]


def coarsen_state(m, mask, ids, c):
    """Block average by c (for example dx 1.5 -> 3: the same analysis grid on both meshes)."""
    N = m.shape[0] // c
    b = lambda a: a[:N * c, :N * c, :N * c].reshape(N, c, N, c, N, c, *a.shape[3:])
    mm = (b(m * mask[..., None]).sum(axis=(1, 3, 5)))
    frac = b(mask.astype(np.float32)).mean(axis=(1, 3, 5))
    mk = frac >= 0.5
    mm /= np.maximum(np.linalg.norm(mm, axis=-1, keepdims=True), 1e-12)
    idc = np.full((N, N, N), -1, dtype=np.int8)
    idb = b(ids)
    for p in range(4):
        idc[(b((ids == p).astype(np.float32)).mean(axis=(1, 3, 5)) >= 0.5)] = p
    idc[~mk] = -1
    return mm.astype(np.float32), mk & (idc >= 0), idc


# --- field operations (periodic box) ------------------------------------------

def grad_norm(f, mask, dx):
    """|grad f| per cell from differences between face neighbours in the same particle
    (mean of the forward and backward difference that exist). f: [N,N,N] or [N,N,N,C]."""
    vec = f.ndim == 4
    g2 = np.zeros(mask.shape, dtype=np.float64)
    for ax in range(3):
        acc = np.zeros(mask.shape, dtype=np.float64)
        cnt = np.zeros(mask.shape, dtype=np.float64)
        for s in (-1, 1):
            nb = np.roll(f, s, axis=ax)
            ok = mask & np.roll(mask, s, axis=ax)
            d = (nb - f)
            d2 = (d * d).sum(-1) if vec else d * d
            acc += np.where(ok, d2, 0.0)
            cnt += ok
        g2 += acc / np.maximum(cnt, 1.0)
    return np.where(mask, np.sqrt(g2) / dx, 0.0)


def smooth_m(m, mask, sigma_cells):
    """Gaussian smoothing inside the particles (normalized convolution), then |m| = 1."""
    if sigma_cells <= 0:
        return m.copy()
    w = mask.astype(np.float32)
    num = np.stack([ndimage.gaussian_filter(m[..., i] * w, sigma_cells, mode="wrap") for i in range(3)], -1)
    num /= np.maximum(np.linalg.norm(num, axis=-1, keepdims=True), 1e-12)
    return np.where(mask[..., None], num, 0.0).astype(np.float32)


def texture(st, widths=WIDTHS):
    """T(w) = (1/pi) integral |grad m_w| dV [l_ex^2] for each smoothing width w (in L_eff)."""
    out = {}
    for w in widths:
        ms = smooth_m(st["m"], st["mask"], w * st["Leff"] / st["dx"])
        out[w] = float(grad_norm(ms, st["mask"], st["dx"]).sum() * st["dx"] ** 3 / math.pi)
    return out


def neighbour_sum(q, mask):
    """Sum of q over the 6 face neighbours in the same particle. q: [N,N,N,K]."""
    s = np.zeros_like(q)
    for ax in range(3):
        for sh in (-1, 1):
            ok = (mask & np.roll(mask, sh, axis=ax))[..., None]
            s += np.where(ok, np.roll(q, sh, axis=ax), 0.0)
    return s


def periodic_label(binary):
    """Connected components (face connectivity) in a periodic box."""
    lab, n = ndimage.label(binary)
    if n == 0:
        return lab, 0
    parent = np.arange(n + 1)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for ax in range(3):
        a = np.take(lab, 0, axis=ax).ravel()
        b = np.take(lab, -1, axis=ax).ravel()
        for x, y in set(zip(a[(a > 0) & (b > 0)].tolist(), b[(a > 0) & (b > 0)].tolist())):
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[max(rx, ry)] = min(rx, ry)
    root = np.array([find(i) for i in range(n + 1)])
    uniq, new = np.unique(root, return_inverse=True)
    return new[lab].reshape(lab.shape), len(uniq) - 1


def inscribed_radius(region, dx, margin):
    """Distance of each cell to the nearest cell outside the region (periodic, up to margin cells)."""
    p = np.pad(region, margin, mode="wrap")
    d = ndimage.distance_transform_edt(p)
    sl = tuple(slice(margin, -margin) for _ in range(3))
    return d[sl] * dx


# --- direction mixture (von Mises-Fisher on S^2) -------------------------------

def _logC(k):
    k = np.maximum(k, 1e-8)
    logsinh = np.where(k > 1e-3, k + np.log1p(-np.exp(-2 * k)) - math.log(2.0), np.log(k))
    return np.log(k) - math.log(4 * math.pi) - logsinh


def _loglik(X, mu, kap):
    return _logC(kap)[None, :] + kap[None, :] * (X @ mu.T)


def fit_vmf(X, K, rng, n_iter=200, restarts=3):
    """EM for a K-component vMF mixture. Returns (weights, mu [K,3], kappa [K], loglik)."""
    best = None
    n = len(X)
    for _ in range(restarts):
        mu = [X[rng.integers(n)]]
        for _k in range(1, K):                       # k-means++ on the sphere
            d = np.clip(1.0 - np.max(X @ np.array(mu).T, axis=1), 0.0, None)
            mu.append(X[rng.choice(n, p=d / d.sum())] if d.sum() > 0 else X[rng.integers(n)])
        mu = np.array(mu)
        kap = np.full(K, 10.0)
        pi = np.full(K, 1.0 / K)
        ll_old = -np.inf
        for _it in range(n_iter):
            lp = np.log(pi)[None, :] + _loglik(X, mu, kap)
            mx = lp.max(1, keepdims=True)
            lse = mx[:, 0] + np.log(np.exp(lp - mx).sum(1))
            r = np.exp(lp - lse[:, None])
            ll = float(lse.sum())
            nk = r.sum(0) + 1e-12
            pi = nk / n
            s = r.T @ X
            rn = np.linalg.norm(s, axis=1)
            mu = s / np.maximum(rn[:, None], 1e-12)
            rb = np.clip(rn / nk, 1e-6, 1 - 1e-6)
            kap = np.clip(rb * (3 - rb ** 2) / (1 - rb ** 2), 1e-3, 1e4)
            if ll - ll_old < 1e-7 * abs(ll):
                break
            ll_old = ll
        if best is None or ll > best[3]:
            best = (pi, mu, kap, ll)
    return best


def select_mixture(X, k_max, rng):
    """vMF mixtures K = 1 ... k_max, BIC per K; the smallest BIC wins."""
    n = len(X)
    fits, bic = {}, {}
    for K in range(1, k_max + 1):
        f = fit_vmf(X, K, rng)
        fits[K] = f
        bic[K] = -2.0 * f[3] + (4 * K - 1) * math.log(n)
    K = min(bic, key=bic.get)
    return K, fits[K], bic


# --- segmentation ------------------------------------------------------------

def hmrf(m, mask, mu, kap, logpi, lam, allowed=None, q=None, n_iter=40):
    """Mean-field hidden Markov random field. Returns q [N,N,N,K] (0 outside the particles)."""
    K = len(kap)
    X = m.reshape(-1, 3)
    data = (_logC(kap)[None, :] + kap[None, :] * (X @ mu.T)).reshape(*mask.shape, K) + logpi
    if allowed is not None:
        data = np.where(allowed, data, -np.inf)
    if q is None:
        q = _softmax(data)
    for _ in range(n_iter):
        qn = _softmax(data + lam * neighbour_sum(q, mask))
        qn = np.where(mask[..., None], qn, 0.0)
        if np.abs(qn - q).max() < 1e-4:
            q = qn
            break
        q = 0.5 * q + 0.5 * qn                        # damped update
    return q.astype(np.float32)


def _softmax(a):
    mx = np.max(a, axis=-1, keepdims=True)
    mx = np.where(np.isfinite(mx), mx, 0.0)
    e = np.exp(a - mx)
    return e / np.maximum(e.sum(-1, keepdims=True), 1e-30)


def segment(st, mix, lam, r_min_leff=1.25, max_rounds=12):
    """HMRF with the thickness rule: a domain is thicker than a wall (largest inscribed radius
    >= r_min = 1.25 L_eff, thickness >= 2.5 L_eff; the centre layer of a 180 deg wall, where m
    points sideways, is 1.8 L_eff thick).
    Stage 1: drop whole direction clusters whose thickest region is thinner than r_min, one per
             round, the thinnest first, and segment again (the remaining domains grow back).
    Stage 2: forbid the thin regions (specks) of the remaining clusters.
    Returns q [N,N,N,K] for all K clusters of mix (0 for dropped clusters)."""
    pi, mu, kap = (np.asarray(x, dtype=np.float64) for x in mix)
    mask, dx = st["mask"], st["dx"]
    r_min = r_min_leff * st["Leff"]
    margin = int(math.ceil(r_min / dx)) + 2
    keep = np.ones(len(kap), dtype=bool)
    for _ in range(max_rounds):
        sel = np.nonzero(keep)[0]
        q = hmrf(st["m"], mask, mu[sel], kap[sel], np.log(np.maximum(pi[sel], 1e-12)), lam)
        lab = np.argmax(q, axis=-1)
        rmax_k, thin = [], []
        for j in range(len(sel)):
            reg = mask & (lab == j)
            if not reg.any():
                rmax_k.append(0.0)
                thin.append(None)
                continue
            comp, n = periodic_label(reg)
            rad = inscribed_radius(reg, dx, margin)
            rmax = np.atleast_1d(np.asarray(ndimage.maximum(rad, comp, index=np.arange(1, n + 1))))
            rmax_k.append(float(rmax.max()))
            thin.append(np.isin(comp, np.nonzero(rmax < r_min)[0] + 1))
        rmax_k = np.array(rmax_k)
        if len(sel) > 1 and (rmax_k < r_min).any():
            keep[sel[np.argmin(rmax_k)]] = False
            continue
        break
    allowed = np.ones(mask.shape + (len(sel),), dtype=bool)
    for j, t in enumerate(thin):
        if t is not None:
            allowed[..., j] &= ~t
    none = ~allowed.any(-1) & mask                     # keep at least the best label everywhere
    if none.any():
        ii = np.nonzero(none)
        allowed[ii + (np.argmax(q[ii], -1),)] = True
    q = hmrf(st["m"], mask, mu[sel], kap[sel], np.log(np.maximum(pi[sel], 1e-12)), lam, allowed, q)
    full = np.zeros(mask.shape + (len(kap),), dtype=np.float32)
    full[..., sel] = q
    return full


def domains_of(q, mask, dx):
    """Domains (connected label regions): list of (label, volume l_ex^3)."""
    lab = np.argmax(q, axis=-1)
    out = []
    for k in range(q.shape[-1]):
        comp, n = periodic_label(mask & (lab == k))
        if n:
            vols = np.bincount(comp.ravel())[1:] * dx ** 3
            out += [(k, float(v)) for v in vols]
    return out


def wall_area_coarea(q, mask, dx):
    return float(0.5 * sum(grad_norm(q[..., k], mask, dx).sum() for k in range(q.shape[-1])) * dx ** 3)


def wall_area_lines(q, mask, dx, rng, n_lines=4000, persist_cells=2.0):
    """Stereology: S_V = 2 P_L on random lines; a label change counts if the new label holds
    for persist_cells cells (suppresses flicker at a boundary)."""
    N = mask.shape[0]
    lab = np.argmax(q, axis=-1)
    step = 0.5
    L = 2 * N
    n_s = int(L / step)
    hits, length = 0, 0.0
    hold = max(1, int(round(persist_cells / step)))
    for _ in range(n_lines):
        v = rng.normal(size=3)
        v /= np.linalg.norm(v)
        p0 = rng.uniform(0, N, size=3)
        pts = p0[None, :] + np.arange(n_s)[:, None] * step * v[None, :]
        idx = np.floor(pts).astype(int) % N
        inside = mask[idx[:, 0], idx[:, 1], idx[:, 2]]
        ll = lab[idx[:, 0], idx[:, 1], idx[:, 2]]
        length += inside.sum() * step
        cur, run_lab, run_n = None, None, 0
        for i in range(n_s):
            if not inside[i]:
                cur, run_lab, run_n = None, None, 0
                continue
            if cur is None:
                cur = ll[i]
                continue
            if ll[i] == cur:
                run_lab, run_n = None, 0
                continue
            if ll[i] == run_lab:
                run_n += 1
            else:
                run_lab, run_n = ll[i], 1
            if run_n >= hold:
                hits += 1
                cur, run_lab, run_n = ll[i], None, 0
    P_L = hits / max(length * dx, 1e-12)
    return float(2.0 * P_L * mask.sum() * dx ** 3)


def switched_volume(qa, qb, mask, dx):
    return float(0.5 * np.abs(qa - qb).sum(-1)[mask].sum() * dx ** 3)


def ripple_rms_deg(st, q):
    """rms angle between m and m smoothed over 2 L_eff, in domain cores (q_max > 0.95)."""
    ms = smooth_m(st["m"], st["mask"], 2.0 * st["Leff"] / st["dx"])
    core = st["mask"] & (q.max(-1) > 0.95)
    if not core.any():
        return float("nan")
    c = np.clip((st["m"][core] * ms[core]).sum(-1), -1, 1)
    return float(np.degrees(np.sqrt(np.mean(np.arccos(c) ** 2))))


# --- analysis of one or several states of one run -------------------------------

def fit_directions(states, rng, k_max=6, n_fit=40000, sigma_fit_leff=1.0):
    X = []
    for st in states:
        ms = smooth_m(st["m"], st["mask"], sigma_fit_leff * st["Leff"] / st["dx"])
        x = ms[st["mask"]]
        X.append(x[rng.choice(len(x), size=min(len(x), n_fit // len(states)), replace=False)])
    X = np.concatenate(X).astype(np.float64)
    return select_mixture(X, k_max, rng)


def analyze(states, lams=LAMS, k_max=6, seed=0, lines=True):
    """States of one run share the direction clusters (labels comparable between phases)."""
    rng = np.random.default_rng(seed)
    K, (pi, mu, kap, ll), bic = fit_directions(states, rng, k_max)
    res = {"K": K, "bic": bic, "mu": mu.tolist(), "kappa": kap.tolist(), "weights": pi.tolist(),
           "states": []}
    for st in states:
        V = float(st["mask"].sum() * st["dx"] ** 3)
        r = {"name": st["name"], "V": V, "texture": texture(st), "lam": {}}
        for lam in lams:
            q = segment(st, (pi, mu, kap), lam)
            d = domains_of(q, st["mask"], st["dx"])
            r["lam"][lam] = {"A_coarea": wall_area_coarea(q, st["mask"], st["dx"]),
                             "A_lines": wall_area_lines(q, st["mask"], st["dx"], rng) if lines else float("nan"),
                             "n_domains": len(d), "domain_volumes": sorted((v for _, v in d), reverse=True),
                             "certainty": float(q.max(-1)[st["mask"]].mean())}
            if lam == lams[len(lams) // 2]:
                r["ripple_rms_deg"] = ripple_rms_deg(st, q)
                r["_q"] = q
        res["states"].append(r)
    qs = [r.pop("_q") for r in res["states"]]
    res["switched"] = [switched_volume(qs[i], qs[(i + 1) % len(qs)], states[i]["mask"], states[i]["dx"])
                       for i in range(len(qs))] if len(qs) > 1 else []
    return res


def report(res):
    print("direction clusters: K = %d (BIC %s)" % (res["K"], ", ".join(
        "K%d %.4g" % (k, v) for k, v in res["bic"].items())))
    for w, mu, ka in zip(res["weights"], res["mu"], res["kappa"]):
        print("   weight %.2f  mu (%+.2f %+.2f %+.2f)  kappa %.1f" % (w, *mu, ka))
    for r in res["states"]:
        print("\n%s   V = %.4g l_ex^3   ripple rms %.1f deg" % (r["name"], r["V"], r.get("ripple_rms_deg", float("nan"))))
        print("   texture T(w) [l_ex^2]: " + "  ".join("w=%g: %.4g" % (w, t) for w, t in r["texture"].items()))
        print("   %6s %12s %12s %9s %10s" % ("lam", "A_coarea", "A_lines", "domains", "certainty"))
        for lam, x in r["lam"].items():
            print("   %6g %12.4g %12.4g %9d %10.3f" % (lam, x["A_coarea"], x["A_lines"], x["n_domains"], x["certainty"]))
    if res["switched"]:
        print("\nswitched volume between consecutive states [l_ex^3]: " +
              "  ".join("%.4g" % v for v in res["switched"]))


# --- synthetic states with known walls ------------------------------------------

def synthetic_state(kind, dx, Leff=12.0, d=300.0, ripple_deg=0.0, seed=0, normal=(0.31, 0.52, 0.80)):
    """One sphere (diameter d) in a box; walls with the profile theta = 2 atan(exp(s/Leff))
    (width pi L_eff). Returns (state, true wall area, true wall angles).
      'ripple' : uniform m + ripple, no wall
      'slab180': 3 parallel 180 deg walls (oblique normal)
      'slab90' : 3 parallel 90 deg walls
      'bubble' : a spherical inner domain (radius 0.3 d), 180 deg wall"""
    a = d + 2 * 12.0
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
    prof = lambda s: 2.0 * np.arctan(np.exp(s / Leff))          # 0 -> pi across a wall at s = 0
    area = 0.0
    if kind == "ripple":
        th = np.zeros_like(X)
    elif kind in ("slab180", "slab90"):
        s = X * n[0] + Y * n[1] + Z * n[2]
        pos = (-0.45 * R, 0.0, 0.4 * R)
        f = 1.0 if kind == "slab180" else 0.5
        th = sum(f * prof(s - p) for p in pos)
        area = sum(math.pi * (R ** 2 - p ** 2) for p in pos)
    elif kind == "bubble":
        r0 = 0.3 * d
        th = prof(r0 - rr)
        area = 4 * math.pi * r0 ** 2
    else:
        raise ValueError(kind)
    m = np.cos(th)[..., None] * u + np.sin(th)[..., None] * v
    if ripple_deg > 0:
        rng = np.random.default_rng(seed)
        nc = int(math.ceil(a / Leff))
        cube = rng.normal(size=(nc, nc, nc, 3))
        ic = np.minimum((np.arange(N) * dx / Leff).astype(int), nc - 1)
        g = cube[np.ix_(ic, ic, ic)]
        g = np.stack([ndimage.gaussian_filter(g[..., i], 0.5 * Leff / dx, mode="wrap") for i in range(3)], -1)
        g -= (g * m).sum(-1, keepdims=True) * m                    # perpendicular to m
        g /= np.sqrt(np.mean((g ** 2).sum(-1)[mask]))
        m = m + math.tan(math.radians(ripple_deg)) * g
        m /= np.linalg.norm(m, axis=-1, keepdims=True)
    m = np.where(mask[..., None], m, 0.0).astype(np.float32)
    ids = np.where(mask, 0, -1).astype(np.int8)
    st = {"m": m, "mask": mask, "ids": ids, "dx": dx, "Leff": Leff,
          "name": "%s dx=%g ripple=%g" % (kind, dx, ripple_deg)}
    return st, area


def selftest(dxs=(3.0, 1.5), ripples=(0.0, 10.0, 20.0, 30.0), kinds=("ripple", "slab180", "slab90", "bubble"),
             d=300.0, lams=LAMS, coarsen_fine=True):
    rows = []
    for kind in kinds:
        for dx in dxs:
            for rp in ripples:
                st, A = synthetic_state(kind, dx, d=d, ripple_deg=rp)
                if coarsen_fine and dx < 3.0:              # same analysis grid as dx = 3
                    c = int(round(3.0 / dx))
                    m, mk, ids = coarsen_state(st["m"], st["mask"], st["ids"], c)
                    st = dict(st, m=m, mask=mk, ids=ids, dx=dx * c)
                res = analyze([st], lams=lams, lines=True)
                r = res["states"][0]
                for lam, x in r["lam"].items():
                    rows.append((kind, dx, rp, lam, res["K"], A, x["A_coarea"], x["A_lines"], x["n_domains"],
                                 r["texture"][2.0], r.get("ripple_rms_deg", float("nan"))))
                print("%-8s dx %-4g ripple %4g  K %d  A_true %9.4g  A_coarea(lam) %s  A_lines %s  T(2) %.4g  ripple_est %.1f"
                      % (kind, dx, rp, res["K"], A,
                         " ".join("%.4g" % x["A_coarea"] for x in r["lam"].values()),
                         " ".join("%.4g" % x["A_lines"] for x in r["lam"].values()),
                         r["texture"][2.0], r.get("ripple_rms_deg", float("nan"))), flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("states", nargs="*", help="m_*.vti or init.pt / checkpoint.pt of ONE run")
    ap.add_argument("--coarsen", type=int, default=1, help="block average (2 for dx 1.5 -> analysis grid dx 3)")
    ap.add_argument("--k_max", type=int, default=6)
    ap.add_argument("--json", default=None, help="write the result as json")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--selftest_csv", default=None)
    a = ap.parse_args()
    if a.selftest:
        rows = selftest()
        if a.selftest_csv:
            import csv
            with open(a.selftest_csv, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["kind", "dx", "ripple_deg", "lam", "K", "A_true", "A_coarea", "A_lines",
                            "n_domains", "T_2Leff", "ripple_est_deg"])
                w.writerows(rows)
        return
    if not a.states:
        sys.exit("no states given")
    states = [load_state(p, a.coarsen) for p in a.states]
    res = analyze(states, k_max=a.k_max)
    report(res)
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
