"""K1 investigation (PROTOKOLL 7.135): per-particle descriptors of the prepared states (init.pt) and of the state
after the 9 mT stage (ckpt_00_amp00.pt) of the K1 runs. Runs on the instance (CPU), reads only.
    python3 eval_k1_states.py <runs dir> <code dir> K1_0_sat K1_1_sat ...
"""
import json
import math
import pathlib
import sys

import numpy as np
import torch

runs, code = pathlib.Path(sys.argv[1]), sys.argv[2]
sys.path.insert(0, code)
from geometry import fcc_box, FCC_SITES          # noqa: E402
from anisotropy_model import particle_axes       # noqa: E402

out = {}
for name in sys.argv[3:]:
    cfg = json.loads((runs / name / "config.json").read_text())
    N, d, dx = cfg["N_used"], cfg["d"], cfg["dx"]
    ids = fcc_box(N, d, dx)["ids"]
    u = particle_axes(cfg.get("direction", [1.0, 0.0, 0.0]))
    a = N * dx
    x = (np.arange(N) + 0.5) * dx
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    rec = {}
    for f in ("init.pt", "ckpt_00_amp00.pt"):
        m = torch.load(str(runs / name / f), map_location="cpu", weights_only=True)["m"].double().numpy()
        parts = []
        for p in range(4):
            s = FCC_SITES[p] * a
            dd = [g - c for g, c in zip((X, Y, Z), s)]
            dd = [v - a * np.round(v / a) for v in dd]
            r = np.sqrt(sum(v**2 for v in dd)) / (0.5 * d)
            sel = ids == p
            core = sel & (r < 0.5)
            mu = m[sel] @ u[p]
            muc = m[core] @ u[p]
            mp = m[sel].mean(0)
            parts.append({
                "cos_axis_drive": float(abs(u[p][0])),
                "m_mean": [round(float(v), 4) for v in mp],
                "m_u_mean": round(float(mu.mean()), 4),
                "share_u_pos": round(float((mu > 0).mean()), 4),
                "share_wall": round(float((np.abs(mu) < 0.5).mean()), 4),
                "core_m_u_mean": round(float(muc.mean()), 4),
                "core_share_u_pos": round(float((muc > 0).mean()), 4),
                "core_share_wall": round(float((np.abs(muc) < 0.5).mean()), 4),
                "core_mx_mean": round(float(m[core][:, 0].mean()), 4),
            })
        rec[f] = parts
    out[name] = rec
print(json.dumps(out))
