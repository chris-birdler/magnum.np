"""Base rate for reading 3a: share of the counted cubes (all 8 corners in one particle, run_loops.bloch_mask_t) that
contain a cut cell (0 < f < 1) or a thin cell (f < l_ex/dx). Usage: cube_base.py <run dir with config.json>"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from geometry import fcc_box
cfg = json.load(open(sys.argv[1] + "/config.json"))
N = cfg["N_used"]; q = cfg.get("staircase_q") if cfg.get("staircase", "none") != "none" else 8
g = fcc_box(N, cfg["d"], cfg["dx"], sub=q)
ids, f = g["ids"].astype(np.int64), g["frac"].astype(np.float64)
lim = cfg["units"]["l_ex"] / cfg["dx"]
ok = ids >= 0; cut = np.zeros_like(ok); thin = np.zeros_like(ok)
for s in ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)):
    r = lambda a: np.roll(a, shift=(-s[0], -s[1], -s[2]), axis=(0, 1, 2))
    if s != (0, 0, 0):
        ok &= r(ids) == ids
    cut |= r((ids >= 0) & (f < 1)); thin |= r((ids >= 0) & (f < lim))
print("N %d q %d counted cubes %d | with a cut cell %.2f %% | with a thin cell %.2f %% | cut-cell share of magnetic cells %.2f %%"
      % (N, q, ok.sum(), 100 * (ok & cut).sum() / ok.sum(), 100 * (ok & thin).sum() / ok.sum(),
         100 * ((ids >= 0) & (f < 1)).sum() / (ids >= 0).sum()))
