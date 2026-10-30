"""Pair coupling c_ij = (l_ex/dx)^2 * 2 f_i f_j / (f_i + f_j) (the field per unit |m_j - m_i| on a cell of f = 1 scale;
here the symmetric coupling coefficient A_ij/A in units of (l_ex/dx)^2) of the pairs above 90 deg against all pairs
with a cut cell. Usage: pair_coupling.py <state> <run dir> ..."""
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
    m = torch.load(run + "/" + state, map_location="cpu", weights_only=False)["m"].double().numpy()
    big, allc = [], []
    for ax in range(3):
        fo = np.roll(f, 1, axis=ax); same = (ids >= 0) & (ids == np.roll(ids, 1, axis=ax))
        cut = same & ((f < 1) | (fo < 1))
        a = 2 * f * fo / np.maximum(f + fo, 1e-30)
        c = (m * np.roll(m, 1, axis=ax)).sum(-1)
        allc.append(a[cut]); big.append(a[same & (c < 0)])
    allc, big = np.concatenate(allc), np.concatenate(big)
    pct = [(allc <= b).mean() * 100 for b in big]
    print("%s pairs>90: %d  A_ij/A of these: %s  | percentile among %d cut pairs: %s  | all cut pairs: median %.3f p5 %.3f"
          % (run.rstrip("/").split("/")[-1], len(big), np.round(np.sort(big), 3).tolist(), len(allc),
             np.round(sorted(pct), 1).tolist(), np.median(allc), np.percentile(allc, 5)))
