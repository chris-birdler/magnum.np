"""
R3 of PLAN.md: smallest and largest eigenvalue of the corrected stray-field operator x -> −v ⊙ h(x) on the magnetic
V1 cells (h = FFT field of the moments v m plus the near-field correction). The operator is positive semi-definite
in the exact case. Plain Lanczos with a fixed number of steps (D7): "not usable" if the smallest Ritz value (an
upper bound of λ_min) is below −10⁻⁴ λ_max; "passed" if it is above and its residual is at most 10⁻⁴ λ_max; else
"not determined".

    python check_sign.py --run K4_1_nx_pd300 --q 2 --p 1 2 3 --forms P L [--shells 1] --out results
"""
import argparse
import json
import os
import pathlib
import sys
import time

import numpy as np
import torch

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "beta_bpeak"))
os.environ.setdefault("CUDA_DEVICE", "-1")

from magnumnp import set_precision  # noqa: E402
set_precision("double")
torch.set_num_threads(int(os.environ.get("STAIRCASE_THREADS", os.cpu_count())))

from anisotropy_model import box_edge  # noqa: E402
from geometry import fcc_box  # noqa: E402
from nearfield import NearCorrection, FFTField, v1_cells, centroid_offsets, centroid_field  # noqa: E402
from check_static import STATES, CACHE, log  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="a run of the geometry (only its config is used)")
    ap.add_argument("--q", type=int, default=2)
    ap.add_argument("--p", type=int, nargs="+", default=[1, 2, 3])
    ap.add_argument("--forms", nargs="+", default=["P", "L"])
    ap.add_argument("--shells", type=int, default=None)
    ap.add_argument("--steps", type=int, default=300, help="Lanczos steps (field evaluations)")
    ap.add_argument("--centroid", action="store_true", help="far field with the first moment (PROTOKOLL 3.7)")
    ap.add_argument("--out", default=str(HERE / "results"))
    a = ap.parse_args()
    cfg = json.loads((STATES / a.run / "config.json").read_text())
    a_lex, d_lex = box_edge(cfg["d_lex"], cfg["phi"], cfg["a_quant"])
    N = int(round(a_lex / cfg["dx_lex"]))
    q = a.q
    n, pid, occ = v1_cells(fcc_box(q * N, d_lex, cfg["dx_lex"] / q)["ids"], q)
    v = n / q ** 3
    mag = torch.as_tensor(n.reshape(-1) > 0)
    vt = torch.as_tensor(v.reshape(-1))[mag][:, None]
    nm = int(mag.sum())
    setting = ("newell_ext", cfg["strayfield_shells"]) if a.shells is None else ("newell", a.shells)
    coarse = FFTField(N, setting, cache_dir=CACHE)
    geo_key = "d%.0f_phi%.2f" % (cfg["d_lex"], cfg["phi"])
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    fname = out / ("sign_%s_q%d%s%s.json" % (geo_key, q, "" if a.shells is None else "_k%d" % a.shells,
                                             "_c" if a.centroid else ""))
    vfull = torch.as_tensor(v)
    dcen = centroid_offsets(occ, q) if a.centroid else None
    res = json.loads(fname.read_text()) if fname.exists() else {}
    for p in a.p:
        nc = NearCorrection(n, occ, q, p, periodic=True, centroid=a.centroid)
        for form in a.forms:
            corr = {"P": nc.h_pair, "L": nc.h_local, "N": torch.zeros_like}[form]   # N: no window term (reference)
            count = [0]

            def mv(x):
                count[0] += 1
                m = torch.zeros(N ** 3, 3, dtype=torch.float64)
                m[mag] = torch.as_tensor(np.asarray(x, dtype=np.float64)).reshape(-1, 3)
                m = m.reshape(N, N, N, 3)
                h = coarse.h(v, m) + corr(m)
                if a.centroid:
                    h = h + centroid_field(coarse.term._N, vfull, dcen, m, shape=coarse.fft_shape)
                return (-(vt * h.reshape(-1, 3)[mag])).reshape(-1).numpy()

            t = time.time()
            key = "%s_q%d_%s%d" % (geo_key, q, form, p)
            # plain Lanczos (no reorthogonalisation) with a fixed number of steps: the extreme Ritz values of the
            # tridiagonal matrix converge first; the smallest Ritz value is an upper bound of λ_min, its residual
            # |β_m s_m| measures the convergence (D7, PROTOKOLL 3.10)
            rng = np.random.default_rng(1)
            x = rng.normal(size=3 * nm)
            x /= np.linalg.norm(x)
            x_old = np.zeros_like(x)
            al, be = [], [0.0]
            for k in range(a.steps):
                w = mv(x) - be[-1] * x_old
                al.append(float(w @ x))
                w -= al[-1] * x
                b_new = float(np.linalg.norm(w))
                be.append(b_new)
                x_old, x = x, w / b_new
                if (k + 1) % 25 == 0:
                    ev = np.linalg.eigvalsh(np.diag(al) + np.diag(be[1:-1], 1) + np.diag(be[1:-1], -1))
                    log(key, "step", k + 1, "ritz min %.4e max %.4e" % (ev[0], ev[-1]))
            T = np.diag(al) + np.diag(be[1:-1], 1) + np.diag(be[1:-1], -1)
            ev, S = np.linalg.eigh(T)
            lmin, lmax = float(ev[0]), float(ev[-1])
            r_min, r_max = abs(be[-1] * S[-1, 0]), abs(be[-1] * S[-1, -1])
            ratio = lmin / lmax
            conv = r_min <= 1e-4 * lmax
            if ratio < -1e-4:
                reading = "not usable"                # an upper bound of λ_min below the limit is a proof
            elif conv:
                reading = "passed"
            else:
                reading = "not determined"
            res[key] = {"lmin_ritz": lmin, "lmax_ritz": lmax, "ratio": ratio, "resid_min": r_min, "resid_max": r_max,
                        "converged": bool(conv), "reading": reading, "matvecs": count[0], "time_s": time.time() - t,
                        "setting": list(setting), "method": "Lanczos %d steps" % a.steps}
            log(key, res[key])
            fname.write_text(json.dumps(res, indent=1))
        del nc


if __name__ == "__main__":
    main()
