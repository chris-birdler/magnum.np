"""Tilt of the cut cells against the bulk: angle between m_i and the sum of m over the full cells (f = 1) of the same
particle within 2 cells (5 x 5 x 5 block). Usage: tilt_ref.py <state> <run dir> ..."""
import json, os, sys
import numpy as np, torch
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from geometry import fcc_box
from pairs_loc import sub_of
state = sys.argv[1]
for run in sys.argv[2:]:
    cfg = json.load(open(run + "/config.json"))
    N = cfg["N_used"]; g = fcc_box(N, cfg["d"], cfg["dx"], sub=sub_of(cfg))
    ids, f = g["ids"].astype(int), g["frac"].astype(float)
    mag = ids >= 0; full = mag & (f >= 1); lim = cfg["units"]["l_ex"] / cfg["dx"]
    m = torch.load(run + "/" + state, map_location="cpu", weights_only=False)["m"].double().numpy()
    s = np.zeros_like(m)
    for a in range(-2, 3):
        for b in range(-2, 3):
            for c in range(-2, 3):
                sh = (a, b, c); nb = np.roll(full, sh, (0, 1, 2)) & (np.roll(ids, sh, (0, 1, 2)) == ids)
                s += nb[..., None] * np.roll(m, sh, (0, 1, 2))
    n = np.linalg.norm(s, axis=-1)
    ang = np.degrees(np.arccos(np.clip((m * s).sum(-1) / np.maximum(n, 1e-30), -1, 1)))
    out = [run.rstrip("/").split("/")[-1]]
    for lab, sel in (("thin", mag & (f < lim) & (n > 0)), ("cut", mag & (f < 1) & (n > 0)), ("full", full & (n > 0))):
        a_ = ang[sel]
        out.append("%s n %d med %.1f p99 %.1f max %.1f >45 %d >90 %d" % (lab, a_.size, np.median(a_), np.percentile(a_, 99),
                                                                       a_.max(), (a_ > 45).sum(), (a_ > 90).sum()) if a_.size else lab + " -")
    print(" | ".join(out))
