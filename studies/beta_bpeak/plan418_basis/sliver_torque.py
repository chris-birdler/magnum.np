"""Torque balance on the cells of the large-angle pairs (relaxed state, zero applied field).
In equilibrium m x (h_ex + h_demag + h_K) = 0, so |h_demag, perp| = |h_ex, perp| within |h_K| <= 2 Q_eff.
h_ex in units of Ms from the fraction-scaled exchange of magnumnp (A_ij harmonic mean of f_i A and f_j A, divided by
f_i Ms): h_ex,i = (l_ex/dx)^2 sum_j 2 f_j/(f_i + f_j) (m_j - m_i). Usage: sliver_torque.py <state> <run dir> ..."""
import json, os, sys
import numpy as np, torch
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from geometry import fcc_box
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pairs_loc import sub_of

state = sys.argv[1]
for run in sys.argv[2:]:
    cfg = json.load(open(run + "/config.json"))
    N = cfg["N_used"]; g = fcc_box(N, cfg["d"], cfg["dx"], sub=sub_of(cfg))
    ids, f = g["ids"].astype(np.int64), g["frac"].astype(np.float64)
    m = torch.load(run + "/" + state, map_location="cpu", weights_only=False)["m"].double().numpy()
    c2 = (cfg["units"]["l_ex"] / cfg["dx"]) ** 2
    h = np.zeros_like(m)
    for ax in range(3):
        for sh in (1, -1):
            fj = np.roll(f, sh, axis=ax); same = (ids >= 0) & (ids == np.roll(ids, sh, axis=ax))
            w = np.where(same & (f > 0), 2 * fj / np.maximum(f + fj, 1e-30), 0.0)
            h += w[..., None] * (np.roll(m, sh, axis=ax) - m)
    h *= c2
    hperp = np.linalg.norm(h - (h * m).sum(-1, keepdims=True) * m, axis=-1)
    flag = np.zeros(m.shape[:3], bool)
    for ax in range(3):
        same = (ids >= 0) & (ids == np.roll(ids, 1, axis=ax))
        c = (m * np.roll(m, 1, axis=ax)).sum(-1)
        sel = same & (c < 0.0)
        fo = np.roll(f, 1, axis=ax)
        small_here = sel & (f <= fo)
        small_other = sel & (f > fo)
        flag |= small_here
        flag |= np.roll(small_other, -1, axis=ax)
    mag = ids >= 0
    bulk = mag & (f >= 1)
    sl = mag & (f > 0) & (f < 0.3)
    q = lambda a: "median %.3f  p90 %.3f  max %.3f" % (np.median(a), np.percentile(a, 90), a.max()) if a.size else "-"
    print("%s  Q_eff %.4f (|h_K| <= %.3f Ms)" % (run.split("/")[-1], cfg["anisotropy"]["Q_eff"], 2 * cfg["anisotropy"]["Q_eff"]))
    print("  |h_ex,perp| / Ms   flagged slivers (n %d, f %s): %s" % (flag.sum(), np.round(np.sort(f[flag]), 3).tolist(), q(hperp[flag])))
    print("  |h_ex,perp| / Ms   all cells f < 0.3 (n %d): %s" % (sl.sum(), q(hperp[sl])))
    print("  |h_ex,perp| / Ms   bulk cells f = 1 (n %d): %s" % (bulk.sum(), q(hperp[bulk])))
