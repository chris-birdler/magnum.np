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
