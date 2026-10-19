"""
R3 of PLAN.md: smallest and largest eigenvalue of the corrected stray-field operator x -> −v ⊙ h(x) on the magnetic
V1 cells (h = FFT field of the moments v m plus the near-field correction). The operator is positive semi-definite
in the exact case. Reading: "passed" if λ_min ≥ −10⁻⁴ λ_max; "not usable" if below; "not determined" if Lanczos
(scipy eigsh, ARPACK) does not converge.

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
from scipy.sparse.linalg import LinearOperator, eigsh, ArpackNoConvergence

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "beta_bpeak"))
os.environ.setdefault("CUDA_DEVICE", "-1")

from magnumnp import set_precision  # noqa: E402
set_precision("double")
torch.set_num_threads(int(os.environ.get("STAIRCASE_THREADS", os.cpu_count())))

from anisotropy_model import box_edge  # noqa: E402
from geometry import fcc_box  # noqa: E402
from nearfield import NearCorrection, FFTField, v1_cells  # noqa: E402
from check_static import STATES, CACHE, log  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="a run of the geometry (only its config is used)")
    ap.add_argument("--q", type=int, default=2)
    ap.add_argument("--p", type=int, nargs="+", default=[1, 2, 3])
    ap.add_argument("--forms", nargs="+", default=["P", "L"])
    ap.add_argument("--shells", type=int, default=None)
    ap.add_argument("--maxiter", type=int, default=400)
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
    fname = out / ("sign_%s_q%d%s.json" % (geo_key, q, "" if a.shells is None else "_k%d" % a.shells))
    res = json.loads(fname.read_text()) if fname.exists() else {}
    for p in a.p:
        nc = NearCorrection(n, occ, q, p, periodic=True)
        for form in a.forms:
            corr = nc.h_pair if form == "P" else nc.h_local
            count = [0]

            def mv(x):
                count[0] += 1
                m = torch.zeros(N ** 3, 3, dtype=torch.float64)
                m[mag] = torch.as_tensor(np.asarray(x, dtype=np.float64)).reshape(-1, 3)
                m = m.reshape(N, N, N, 3)
                h = coarse.h(v, m) + corr(m)
                return (-(vt * h.reshape(-1, 3)[mag])).reshape(-1).numpy()

            op = LinearOperator((3 * nm, 3 * nm), matvec=mv, dtype=np.float64)
            t = time.time()
            lmax = float(eigsh(op, k=1, which="LA", tol=1e-3, maxiter=a.maxiter, return_eigenvectors=False)[0])
            try:
                lmin = float(eigsh(op, k=1, which="SA", tol=1e-3, ncv=40, maxiter=a.maxiter,
                                   return_eigenvectors=False)[0])
                conv = True
            except ArpackNoConvergence as e:
                lmin = float(np.min(e.eigenvalues)) if len(e.eigenvalues) else float("nan")
                conv = False
            ratio = lmin / lmax
            reading = ("not determined" if not conv else ("passed" if ratio >= -1e-4 else "not usable"))
            key = "%s_q%d_%s%d" % (geo_key, q, form, p)
            res[key] = {"lmin": lmin, "lmax": lmax, "ratio": ratio, "converged": conv, "reading": reading,
                        "matvecs": count[0], "time_s": time.time() - t, "setting": list(setting)}
            log(key, res[key])
            fname.write_text(json.dumps(res, indent=1))
        del nc


if __name__ == "__main__":
    main()
