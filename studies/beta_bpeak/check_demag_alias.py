"""Measuring-method check of PROTOKOLL 7.105 (python check_demag_alias.py > results/demag_alias.md; needs
runs/B0b_I3_s941/init.pt and runs/B0b_T3_s941/init.pt). Independent exact demag energy of piecewise-constant M on cubic cells in a periodic box (zero mean field):
E / (K_d V) = 1 / (N^6 Ms^2) sum_q M~(q)^H K(q) M~(q),  K(q) = sum_G sinc^2-form factor * t t^T / |t|^2,
t_a = pi (m_a / N + G_a), G truncated at |G_a| <= Gm (k = 0 excluded). Compared with DemagFieldPBC (FD) and Newell with
image shells + compensation, on: (1) aligned cube f = 1/8, (2) B0b uniform (111), (3) relaxed B0b states I3, T3 (s941),
(4) a checkerboard pattern m = (1, -1, 0)/sqrt2 * (-1)^(i+j) in a full box."""
import sys, time, math, numpy as np, torch
import pathlib
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))
sys.path.insert(0, str(HERE))
from magnumnp import Mesh, State, DemagFieldPBC, set_precision
set_precision("double")
import check_mesh_energy as C
from anisotropy_model import box_edge
from geometry import fcc_box
from units import units
U = units(1.5, 1e-11); Ms, Kd, lex = U["Ms"], U["Kd"], U["l_ex"]

def alias_energy(M, Gm):
    """M: [N,N,N,3] in A/m. Returns E/(K_d V)."""
    N = M.shape[0]
    Mt = np.fft.fftn(M, axes=(0, 1, 2)).reshape(-1, 3)
    m1 = np.fft.fftfreq(N, 1.0 / N) / N                      # m_a / N in [-1/2, 1/2)
    QX, QY, QZ = [g.reshape(-1) for g in np.meshgrid(m1, m1, m1, indexing="ij")]
    K = np.zeros((QX.size, 6))
    for gx in range(-Gm, Gm + 1):
        tx = QX + gx; sx = np.sinc(tx) ** 2
        for gy in range(-Gm, Gm + 1):
            ty = QY + gy; sxy = sx * np.sinc(ty) ** 2
            for gz in range(-Gm, Gm + 1):
                tz = QZ + gz; w = sxy * np.sinc(tz) ** 2
                t2 = tx * tx + ty * ty + tz * tz
                ok = t2 > 0
                w = np.where(ok, w / np.where(ok, t2, 1.0), 0.0)
                K[:, 0] += w * tx * tx; K[:, 1] += w * ty * ty; K[:, 2] += w * tz * tz
                K[:, 3] += w * tx * ty; K[:, 4] += w * tx * tz; K[:, 5] += w * ty * tz
    Mx, My, Mz = Mt[:, 0], Mt[:, 1], Mt[:, 2]
    e = (K[:, 0] * abs(Mx) ** 2 + K[:, 1] * abs(My) ** 2 + K[:, 2] * abs(Mz) ** 2
         + 2 * K[:, 3] * (np.conj(Mx) * My).real + 2 * K[:, 4] * (np.conj(Mx) * Mz).real
         + 2 * K[:, 5] * (np.conj(My) * Mz).real).sum()
    return float(e) / (N ** 6 * Ms ** 2)

def fd_energy(ms, m, dx):
    N = ms.shape[0]
    st = State(Mesh((N, N, N), (dx,) * 3)); st.material = {"Ms": torch.tensor(ms[..., None])}
    st.m = torch.tensor(m)
    return float(DemagFieldPBC().E(st)) / (Kd * dx ** 3 * N ** 3)

def newell_energy(ms, m, dx, shells):
    st, term = C.newell_pbc(ms.shape[0], dx, torch.tensor(ms[..., None]), shells)
    return C.energy_newell_pbc(st, term, torch.tensor(m), Ms, Kd)

def report(name, ms, m, dx, exact=None, shells=(1, 2), gms=(4, 8, 16)):
    t = time.time()
    M = ms[..., None] * m
    al = [alias_energy(M, g) for g in gms]
    ext = al[-1] + (al[-1] - al[-2]) * gms[-2] / (gms[-1] - gms[-2])     # 1/Gm extrapolation
    fd = fd_energy(ms, m, dx)
    nw = [newell_energy(ms, m, dx, s) for s in shells]
    print("%-28s alias Gm %s: %s -> extrap. %.7f | FD %.7f | Newell %s | exact %s  [%.0f s]" % (
        name, gms, " ".join("%.7f" % a for a in al), ext, fd, " ".join("%.7f" % x for x in nw),
        "%.7f" % exact if exact is not None else "-", time.time() - t), flush=True)

def legacy():
    """The report of PROTOKOLL 7.105 (results/demag_alias.md)."""
    # (1) aligned cube
    n = 16; ms = np.zeros((n, n, n)); ms[:8, :8, :8] = Ms
    report("cube f=1/8, m (111)", ms, np.broadcast_to(np.ones(3) / math.sqrt(3), (n, n, n, 3)).copy(), lex, exact=7 / 192)
    # (4) checkerboard in a full box
    n = 16; ms = np.full((n, n, n), Ms); i, j, k = np.meshgrid(*[np.arange(n)] * 3, indexing="ij")
    s = (-1.0) ** (i + j); m = np.stack([s, -s, 0 * s], -1) / math.sqrt(2)
    report("checkerboard (1,-1,0)", ms, m, lex)
    # (2), (3) B0b dx 3
    a_lex, d_used = box_edge(96.0, 0.65, 12.0); N = 48; dx = 3 * lex
    chi = fcc_box(N, d_used * lex, dx)["ids"] >= 0; ms = np.where(chi, Ms, 1e-8 * Ms)
    f = chi.mean()
    report("B0b uniform (111)", ms, np.broadcast_to(np.ones(3) / math.sqrt(3), (N, N, N, 3)).copy(), dx, exact=f * (1 - f) / 3)
    for arm in ("I3", "T3"):
        m = torch.load(HERE / "runs" / ("B0b_%s_s941" % arm) / "init.pt", map_location="cpu", weights_only=False)["m"].double().numpy()
        report("B0b %s s941 (box: FD/N4)" % arm, ms, m, dx, shells=(1, 2, 4))


# ---------------------------------------------------------------------------------------------------------------------
# PLAN 4.12, V1 (decision Chris 2026-10-17): validation of the final stray-field setting against the alias sum.
# Readings fixed before the run (PLAN 4.12):
#   calibration  "reference calibrated" if the cube and the uniform B0b state agree with the exact energies within
#                0.01 %, and the checkerboard energy and field within 0.05 % and 0.1 % RMS; else stop;
#   candidate    "valid" if |dE| + 2 u <= 0.025 % and dH + 2 u_H <= 0.25 % for all 12 N0b states on dx 3;
#                "not valid" if |dE| - 2 u > 0.025 % or dH - 2 u_H > 0.25 % for one state or more;
#                else "not determined"; the final setting is the most accurate valid candidate.
# dE = (E - E_ref) / E_ref; dH = RMS over the magnetic cells of |H - H_ref| / RMS of |H_ref|; E_ref, H_ref: fit
# c0 + c1 / Gm + c2 / Gm^2 over Gm = 8, 16, 32; u, u_H: difference to the linear extrapolation 2 v(32) - v(16).
#     python check_demag_alias.py v1_ref N0b              # reference (CPU), runs/v1_ref_N0b.npz
#     CUDA_DEVICE=0 python check_demag_alias.py v1_newell N0b 1 2 4 2x 4x 8     # candidates, runs/v1_newell_N0b_<c>.npz
#     python check_demag_alias.py v1_report N0b > runs/v1_N0b.md
GMS_V1 = (8, 16, 32)
LIM_E, LIM_H = 0.025e-2, 0.25e-2
CANDIDATES = ("1", "2", "4", "8", "2x", "4x")      # "kx": k shells, kernel extrapolated with the inner k // 2 shells


def alias_kernel(N, Gm):
    """Cell-averaged periodic demag tensor K(q) [N^3, 6] (xx, yy, zz, xy, xz, yz) of the alias sum, |G_a| <= Gm, k = 0
    excluded. It does not depend on M; the gz sum is a vector axis."""
    m1 = np.fft.fftfreq(N, 1.0 / N) / N
    QX, QY, QZ = [g.reshape(-1) for g in np.meshgrid(m1, m1, m1, indexing="ij")]
    gz = np.arange(-Gm, Gm + 1, dtype=float)[:, None]
    tz = QZ[None, :] + gz
    sz = np.sinc(tz) ** 2
    K = np.zeros((QX.size, 6))
    for gx in range(-Gm, Gm + 1):
        tx = QX + gx
        sx = np.sinc(tx) ** 2
        for gy in range(-Gm, Gm + 1):
            ty = QY + gy
            sxy = sx * np.sinc(ty) ** 2
            t2 = tx ** 2 + ty ** 2 + tz ** 2
            ok = t2 > 0
            w = np.where(ok, sxy * sz / np.where(ok, t2, 1.0), 0.0)
            K[:, 0] += (w * tx ** 2).sum(0)
            K[:, 1] += (w * ty ** 2).sum(0)
            K[:, 2] += (w * tz ** 2).sum(0)
            K[:, 3] += (w * tx * ty).sum(0)
            K[:, 4] += (w * tx * tz).sum(0)
            K[:, 5] += (w * ty * tz).sum(0)
    return K


def kernel_trace_error(N, Gm, K):
    """max |tr K(q) - prod_a sum_G sinc^2(t_a)| / prod over q != 0 (the box sum factorises)."""
    m1 = np.fft.fftfreq(N, 1.0 / N) / N
    g = np.arange(-Gm, Gm + 1, dtype=float)
    s1 = (np.sinc(m1[:, None] + g[None, :]) ** 2).sum(1)
    P = (s1[:, None, None] * s1[None, :, None] * s1[None, None, :]).reshape(-1)
    tr = K[:, 0] + K[:, 1] + K[:, 2]
    return float(np.max(np.abs(tr[1:] - P[1:]) / P[1:]))


def alias_energy_field(M, K):
    """(E / (K_d V), H [N, N, N, 3] in A/m) of M [N, N, N, 3] (A/m) with the kernel K; H = −K M~ in Fourier space."""
    N = M.shape[0]
    Mt = np.fft.fftn(M, axes=(0, 1, 2)).reshape(-1, 3)
    Ht = -np.stack([K[:, 0] * Mt[:, 0] + K[:, 3] * Mt[:, 1] + K[:, 4] * Mt[:, 2],
                    K[:, 3] * Mt[:, 0] + K[:, 1] * Mt[:, 1] + K[:, 5] * Mt[:, 2],
                    K[:, 4] * Mt[:, 0] + K[:, 5] * Mt[:, 1] + K[:, 2] * Mt[:, 2]], -1)
    e = -(np.conj(Mt) * Ht).sum().real / (N ** 6 * Ms ** 2)
    H = np.fft.ifftn(Ht.reshape(N, N, N, 3), axes=(0, 1, 2)).real
    return float(e), H


def gm_fit(v8, v16, v32):
    """(fit c0 of c0 + c1 / Gm + c2 / Gm^2 over Gm = 8, 16, 32; linear extrapolation 2 v32 − v16), elementwise."""
    c0 = (v8 - 6.0 * v16 + 8.0 * v32) / 3.0
    return c0, 2.0 * v32 - v16


def rms_rel(H, Href, mask):
    d = np.sqrt(((H - Href) ** 2).sum(-1))[mask]
    r = np.sqrt((Href ** 2).sum(-1))[mask]
    return float(np.sqrt((d ** 2).mean()) / np.sqrt((r ** 2).mean()))


def v1_states(prefix):
    """[(name, ms [N,N,N] A/m, m [N,N,N,3], mask, dx, exact energy or None)]: the calibration cases and the 12 states."""
    import run_loops
    out = []
    n = 16
    ms = np.zeros((n, n, n))
    ms[:8, :8, :8] = Ms
    out.append(("cube f=1/8 (111)", ms, np.broadcast_to(np.ones(3) / math.sqrt(3), (n, n, n, 3)).copy(),
                ms > 0, lex, 7 / 192))
    ms = np.full((n, n, n), Ms)
    i, j, k = np.meshgrid(*[np.arange(n)] * 3, indexing="ij")
    s = (-1.0) ** (i + j)
    out.append(("checkerboard (1,-1,0)", ms, np.stack([s, -s, 0 * s], -1) / math.sqrt(2), ms > 0, lex, 0.5))
    S3, S075 = C.setup(96.0, 3.0, C.SEED), C.setup(96.0, 0.75, C.SEED)
    ms3 = S3[0].material["Ms"][..., 0].double().cpu().numpy()
    dx3 = float(S3[0].mesh.dx[0])
    mag3 = S3[3] >= 0
    out.append(("B0b geometry uniform (111)", ms3, np.broadcast_to(np.ones(3) / math.sqrt(3), ms3.shape + (3,)).copy(),
                mag3, dx3, float((ms3 / Ms).var()) / 3.0))
    for seed in (941, 942, 943, 944):
        mF = torch.load(HERE / "runs" / ("%s_F_s%d" % (prefix, seed)) / "init.pt", map_location="cpu",
                        weights_only=False)["m"].double().numpy()
        mc, _ = run_loops.transfer_m(mF, S075[3], S3[3])
        out.append(("%s F→3 s%d" % (prefix, seed), ms3, np.asarray(mc, dtype=np.float64), mag3, dx3, None))
        for arm in ("I3", "T3"):
            m = torch.load(HERE / "runs" / ("%s_%s_s%d" % (prefix, arm, seed)) / "init.pt", map_location="cpu",
                           weights_only=False)["m"].double().numpy()
            out.append(("%s %s s%d" % (prefix, arm, seed), ms3, m, mag3, dx3, None))
    return out


def v1_ref(prefix, gms=GMS_V1, device=None):
    """Alias reference: per state E and H at the Gm set (runs/v1_ref_<prefix>[_<gms>].npz), the code checks. device:
    a torch device for the kernel (alias_kernel_torch), else numpy."""
    t0 = time.time()
    states = v1_states(prefix)
    kern, res = {}, {}
    for name, ms, m, mask, dx, exact in states:
        N = ms.shape[0]
        M = ms[..., None] * m
        vals = []
        for g in gms:
            if (N, g) not in kern:
                kern[(N, g)] = alias_kernel_torch(N, g, device) if device else alias_kernel(N, g)
                print("kernel N %d Gm %d: trace error %.1e [%.0f s]" % (N, g, kernel_trace_error(N, g, kern[(N, g)]),
                                                                        time.time() - t0), file=sys.stderr, flush=True)
            e, H = alias_energy_field(M, kern[(N, g)])
            e_mh = float(-(M * H).sum() / (Ms ** 2 * N ** 3))          # -1/2 mu0 <M.H> / K_d
            vals.append((e, H, e_mh))
        res[name] = vals
        print("%s: E %s, -1/2<M.H>/Kd - E = %.1e [%.0f s]" % (name, " ".join("%.8f" % v[0] for v in vals),
              max(abs(v[2] - v[0]) for v in vals), time.time() - t0), file=sys.stderr, flush=True)
    tag = "" if tuple(gms) == GMS_V1 else "_" + "_".join(map(str, gms))
    np.savez(HERE / "runs" / ("v1_ref_%s%s.npz" % (prefix, tag)),
             **{"E_%d_%d" % (i, g): res[st[0]][j][0] for i, st in enumerate(states) for j, g in enumerate(gms)},
             **{"H_%d_%d" % (i, g): res[st[0]][j][1] for i, st in enumerate(states) for j, g in enumerate(gms)},
             names=np.array([st[0] for st in states]),
             trace_err=np.array([kernel_trace_error(N, g, K) for (N, g), K in kern.items()]))


def v1_newell(prefix, cands):
    """Newell candidates (DemagField(true_pbc=True[, extrapolate])): E and H per state (runs/v1_newell_<prefix>_<c>.npz)."""
    from magnumnp import DemagField
    states = v1_states(prefix)
    for c in cands:
        t0 = time.time()
        k, ext = int(c.rstrip("x")), c.endswith("x")
        terms, out = {}, {}
        for i, (name, ms, m, mask, dx, exact) in enumerate(states):
            N = ms.shape[0]
            if (N, dx) not in terms:
                st = State(Mesh((N, N, N), (dx,) * 3, pbc=(k, k, k)))
                terms[(N, dx)] = (st, DemagField(true_pbc=True, extrapolate=ext))
            st, term = terms[(N, dx)]
            st.material = {"Ms": torch.tensor(ms[..., None], dtype=st.dtype, device=st.device)}
            st.m = torch.tensor(m, dtype=st.dtype, device=st.device)
            out["E_%d" % i] = float(term.E(st)) / (Kd * dx ** 3 * N ** 3)
            out["H_%d" % i] = term.h(st).double().cpu().numpy()
        np.savez(HERE / "runs" / ("v1_newell_%s_%s.npz" % (prefix, c)), **out)
        print("candidate %s done [%.0f s]" % (c, time.time() - t0), file=sys.stderr, flush=True)


def v1_classify(dE, uE, dH, uH):
    """Reading of one candidate from the per-state arrays (PLAN 4.12)."""
    if np.all(np.abs(dE) + 2 * uE <= LIM_E) and np.all(dH + 2 * uH <= LIM_H):
        return "valid"
    if np.any(np.abs(dE) - 2 * uE > LIM_E) or np.any(dH - 2 * uH > LIM_H):
        return "not valid"
    return "not determined"


def v1_report(prefix, gms=GMS_V1):
    states = v1_states(prefix)
    tag = "" if tuple(gms) == GMS_V1 else "_" + "_".join(map(str, gms))
    R = np.load(HERE / "runs" / ("v1_ref_%s%s.npz" % (prefix, tag)))
    L = ["# PLAN 4.12 V1: Newell candidates against the alias sum (check_demag_alias.py v1_report %s)" % prefix, "",
         "Energies per box volume in K_d. Reference: fit c0 + c1 / Gm + c2 / Gm² over Gm = %s; u: difference to "
         "the linear extrapolation from the two largest Gm. dE = (E − E_ref) / E_ref; dH = RMS over the magnetic cells of "
         "|H − H_ref| / RMS |H_ref|. Largest trace error of the alias kernels: %.1e." % (
             ", ".join(map(str, gms)), float(R["trace_err"].max())), "",
         "## Calibration", "", "| case | E(8) | E(16) | E(32) | E_ref | u | exact | E_ref / exact − 1 | dH vs exact |",
         "|---|---|---|---|---|---|---|---|---|"]
    ref = []
    for i, (name, ms, m, mask, dx, exact) in enumerate(states):
        E = [float(R["E_%d_%d" % (i, g)]) for g in gms]
        H = [R["H_%d_%d" % (i, g)] for g in gms]
        e0, el = gm_fit(*E)
        h0, hl = gm_fit(*H)
        ref.append((e0, abs(e0 - el) / abs(e0), h0, rms_rel(hl, h0, mask)))
        if exact is not None:
            dh = "–"
            if name.startswith("checkerboard"):
                dh = "%.4f %%" % (100 * rms_rel(h0, -0.5 * ms[..., None] * m, mask))
            L.append("| %s | %.7f | %.7f | %.7f | %.7f | %.1e | %.7f | %+.4f %% | %s |" % (
                name, E[0], E[1], E[2], e0, ref[-1][1], exact, 100 * (e0 / exact - 1), dh))
    cal = [(s, r) for s, r in zip(states, ref) if s[5] is not None]
    ok = all(abs(r[0] / s[5] - 1) <= (5e-4 if s[0].startswith("checkerboard") else 1e-4) for s, r in cal)
    cb = [(s, r) for s, r in cal if s[0].startswith("checkerboard")][0]
    ok = ok and rms_rel(cb[1][2], -0.5 * cb[0][1][..., None] * cb[0][2], cb[0][3]) <= 1e-3
    L += ["", "**Calibration: %s.**" % ("reference calibrated" if ok else "not calibrated: stop"), ""]
    idx = [i for i, s in enumerate(states) if s[5] is None]
    L += ["## Candidates (%d states)" % len(idx), "",
          "| candidate | max dE | max (dE + 2u) | max dH | max (dH + 2u_H) | reading |", "|---|---|---|---|---|---|"]
    out = {}
    for c in CANDIDATES:
        f = HERE / "runs" / ("v1_newell_%s_%s.npz" % (prefix, c))
        if not f.exists():
            L.append("| %s | not computed | | | | |" % c)
            continue
        Z = np.load(f)
        dE = np.array([(float(Z["E_%d" % i]) - ref[i][0]) / ref[i][0] for i in idx])
        uE = np.array([ref[i][1] for i in idx])
        dH = np.array([rms_rel(Z["H_%d" % i], ref[i][2], states[i][3]) for i in idx])
        uH = np.array([ref[i][3] for i in idx])
        cls = v1_classify(dE, uE, dH, uH)
        out[c] = {"dE": dE.tolist(), "dH": dH.tolist(), "reading": cls}
        L.append("| %s | %+.4f %% | %.4f %% | %.4f %% | %.4f %% | %s |" % (
            c, 100 * dE[np.argmax(np.abs(dE))], 100 * np.max(np.abs(dE) + 2 * uE), 100 * dH.max(),
            100 * np.max(dH + 2 * uH), cls))
    valid = [c for c in out if out[c]["reading"] == "valid"]
    best = min(valid, key=lambda c: max(np.max(np.abs(out[c]["dE"])) / LIM_E, np.max(out[c]["dH"]) / LIM_H)) \
        if (valid and ok) else None
    L += ["", "**Final setting: %s.**" % (best if best else "none (stop and report)")]
    import json
    (HERE / "runs" / ("v1_%s%s.json" % (prefix, tag))).write_text(json.dumps(
        {"calibrated": ok, "candidates": out, "final": best}, indent=1))
    print("\n".join(L))


# PLAN 4.12 V1, finer reference (decision Chris 2026-10-17, PROTOKOLL 7.109 option 1). Reading of the checkerboard
# check, fixed before the run: "finer reference works" if u of the checkerboard with Gm = 32, 64, 128 is at least 8 x
# smaller than with Gm = 8, 16, 32 (prediction about 16 x, u ~ c2 / (2 Gm^2)) and |E_ref / 0.5 − 1| ≤ 0.005 %; else stop.
# Then the 12 states with Gm = 32, 64, 128 on a GPU (v1_ref <prefix> 32 64 128), the readings of PLAN 4.12 unchanged.
GMS_FINE = (32, 64, 128)


def alias_kernel_torch(N, Gm, device="cuda"):
    """alias_kernel on a torch device (float64); the gy and gz sums are vector axes."""
    m1 = torch.fft.fftfreq(N, 1.0 / N, dtype=torch.float64, device=device) / N
    QX, QY, QZ = [g.reshape(-1) for g in torch.meshgrid(m1, m1, m1, indexing="ij")]
    g = torch.arange(-Gm, Gm + 1, dtype=torch.float64, device=device)[:, None]
    tz = QZ[None, :] + g
    sz = torch.sinc(tz) ** 2
    K = torch.zeros((QX.numel(), 6), dtype=torch.float64, device=device)
    for gx in range(-Gm, Gm + 1):
        tx = QX + gx
        sx = torch.sinc(tx) ** 2
        for gy in range(-Gm, Gm + 1):
            ty = QY + gy
            sxy = sx * torch.sinc(ty) ** 2
            t2 = tx ** 2 + ty ** 2 + tz ** 2
            w = torch.where(t2 > 0, sxy * sz / torch.where(t2 > 0, t2, torch.ones_like(t2)), torch.zeros_like(t2))
            K[:, 0] += (w * tx ** 2).sum(0)
            K[:, 1] += (w * ty ** 2).sum(0)
            K[:, 2] += (w * tz ** 2).sum(0)
            K[:, 3] += (w * tx * ty).sum(0)
            K[:, 4] += (w * tx * tz).sum(0)
            K[:, 5] += (w * ty * tz).sum(0)
    return K.cpu().numpy()


def v1_checkerboard(gms_list=(GMS_V1, GMS_FINE)):
    """E_ref, u and the field error of the 16^3 checkerboard for the Gm sets; the reading of the finer reference."""
    n = 16
    i, j, k = np.meshgrid(*[np.arange(n)] * 3, indexing="ij")
    s = (-1.0) ** (i + j)
    M = Ms * np.stack([s, -s, 0 * s], -1) / math.sqrt(2)
    mask = np.ones((n, n, n), bool)
    res = {}
    for gms in gms_list:
        vals = [alias_energy_field(M, alias_kernel(n, g)) for g in gms]
        e0, el = gm_fit(*[v[0] for v in vals])
        h0, hl = gm_fit(*[v[1] for v in vals])
        res[gms] = (e0, abs(e0 - el) / abs(e0), rms_rel(h0, -0.5 * M, mask))
        print("Gm %s: E %s E_ref %.8f (E_ref / 0.5 − 1 %+.5f %%), u %.5f %%, field vs exact %.5f %%" % (
            gms, " ".join("%.8f" % v[0] for v in vals), e0, 100 * (e0 / 0.5 - 1), 100 * res[gms][1],
            100 * res[gms][2]), flush=True)
    r = res[gms_list[0]][1] / res[gms_list[-1]][1]
    ok = r >= 8 and abs(res[gms_list[-1]][0] / 0.5 - 1) <= 5e-5
    print("u ratio %.1f. **Reading: %s.**" % (r, "finer reference works" if ok else "stop"))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "v1_cb":
        v1_checkerboard()
    elif len(sys.argv) > 2 and sys.argv[1] == "v1_ref":
        a = sys.argv[3:]
        dev = a.pop() if a and not a[-1].isdigit() else None
        v1_ref(sys.argv[2], tuple(int(x) for x in a) or GMS_V1, dev)
    elif len(sys.argv) > 2 and sys.argv[1] == "v1_newell":
        v1_newell(sys.argv[2], sys.argv[3:] or CANDIDATES)
    elif len(sys.argv) > 2 and sys.argv[1] == "v1_report":
        v1_report(sys.argv[2], tuple(int(x) for x in sys.argv[3:]) or GMS_V1)
    else:
        legacy()
