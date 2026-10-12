"""Design data for the audit: the 116 alpha 0.1 design runs (S000 ... S095, C0 ... C7, K1_0 ... K4_2) under several
cycle rules. Independent of analyze.py / analyze_design.py.

rules:
  steady2_mc2   main result (7.65): steady tail, >= 2 calm cycles, finished reruns X_mc2 replace X
  steady2       the same without reruns (7.54 check; 13 stages missing)
  steady1       steady tail, >= 1 calm cycle, no reruns (7.54 first variant)
  all           all measured cycles (2 ... 6), no reruns (pre-registered PLAN 4.0)
  last2_mc2     steady tail with reruns, only the last 2 calm cycles
  wdis_mc2      as steady2_mc2 but w = w_dis (LLG dissipation)
"""
import math
import re

import numpy as np

import audit_lib as al

AMPS = np.array(al.AMPS)
LNA = np.log(AMPS)
PAT = re.compile(r"^(S\d{3}|C\d|K\d_\d)$")
RANGES = {"lnQ": (math.log(1 / 30.0 ** 2), math.log(1 / 12.0 ** 2)), "r_p": (0.0, 3.0), "phi": (0.55, 0.69),
          "lnd": (math.log(212.0), math.log(300.0))}
FACT = ("lnQ", "r_p", "phi", "lnd")

RULES = {
    "steady2_mc2": dict(rule="steady", min_calm=2, reruns=True, key="w_loop", last=None),
    "steady2": dict(rule="steady", min_calm=2, reruns=False, key="w_loop", last=None),
    "steady1": dict(rule="steady", min_calm=1, reruns=False, key="w_loop", last=None),
    "all": dict(rule="all", min_calm=1, reruns=False, key="w_loop", last=None),
    "last2_mc2": dict(rule="steady", min_calm=2, reruns=True, key="w_loop", last=2),
    "wdis_mc2": dict(rule="steady", min_calm=2, reruns=True, key="w_dis", last=None),
}


def design_names():
    return sorted(d.name for d in al.RUNS.iterdir() if d.is_dir() and PAT.match(d.name) and (d / "DONE").exists()
                  and abs(al.read(d)[1]["alpha"] - 0.1) < 1e-12)


def to_nominal(B, w):
    """ln w at the 7 nominal amplitudes: linear in ln B between present stages, line through the two outer present
    stages at the ends; nan where the stage itself is missing."""
    B, w = np.asarray(B, float), np.asarray(w, float)
    ok = np.isfinite(w) & (w > 0)
    lb, lw = np.log(B[ok]), np.log(w[ok])
    out = np.interp(LNA, lb, lw)
    lo = LNA < lb[0]
    hi = LNA > lb[-1]
    out[lo] = lw[0] + (lw[1] - lw[0]) / (lb[1] - lb[0]) * (LNA[lo] - lb[0])
    out[hi] = lw[-1] + (lw[-1] - lw[-2]) / (lb[-1] - lb[-2]) * (LNA[hi] - lb[-1])
    # stage present?  (one stage per nominal amplitude)
    nom = [al.nominal(b) for b in B]
    for k, a in enumerate(AMPS):
        idx = [i for i, x in enumerate(nom) if x == a]
        if not idx or not np.any(ok[idx]):
            out[k] = np.nan
    return out


def load(rule="steady2_mc2"):
    p = RULES[rule]
    recs = []
    for n in design_names():
        use = al.use_name(n, p["reruns"])
        cfg, summ, st = al.stages(use, rule=p["rule"], min_calm=p["min_calm"], key=p["key"], last=p["last"])
        if p["key"] == "w_dis":
            for s in st:
                s["w"] = s["w_dis"]
        B = [s["B_mT"] for s in st]
        w = [s["w"] if not s["missing"] else np.nan for s in st]
        if len(st) != 7:
            raise RuntimeError("%s: %d stages" % (use, len(st)))
        lnw = to_nominal(B, w)
        Leff = cfg["Leff_lex"]
        rec = {"name": n, "folder": use, "seed": cfg["seed"], "lnQ": math.log(cfg.get("Q_eff") or 1.0 / Leff ** 2),
               "r_p": cfg["r_p"], "phi": cfg["phi_vox"], "lnd": math.log(cfg["d_lex_used"]), "Leff": Leff,
               "d_lex": cfg["d_lex_used"], "lnw": lnw, "B": np.array(B), "w_stage": np.array(w),
               "n_kept": np.array([s["n_kept"] for s in st]), "missing": np.array([s["missing"] for s in st]),
               "pf": al.p_factor(cfg), "kind": "corner" if n.startswith("K") else ("centre" if n.startswith("C") else "sobol")}
        # reading 1 (d = 1 µm): l_ex = 1 µm / (d/l_ex), A = K_d l_ex^2
        rec["A_pJm"] = 1e12 * cfg["units"]["Kd"] * (1e-6 / cfg["d_lex_used"]) ** 2
        rec["d_nm_A10"] = cfg["d_lex_used"] * 3.3421710327935706       # particle size at A = 10 pJ/m (reading 2)
        rec["K_eff"] = cfg["units"]["Kd"] / Leff ** 2 if not cfg.get("Q_eff") else cfg["units"]["Kd"] * cfg["Q_eff"]
        recs.append(rec)
    return recs


def coded(recs):
    return np.array([[2.0 * (r[f] - RANGES[f][0]) / (RANGES[f][1] - RANGES[f][0]) - 1.0 for f in FACT] for r in recs])


def terms(X, kind="full"):
    n, k = X.shape
    cols = [np.ones(n)] + [X[:, i] for i in range(k)]
    if kind == "full":
        cols += [X[:, i] ** 2 for i in range(k)] + [X[:, i] * X[:, j] for i in range(k) for j in range(i + 1, k)]
    return np.column_stack(cols)


def slope_w(idx):
    x = LNA[list(idx)]
    return (x - x.mean()) / ((x - x.mean()) ** 2).sum()


def interp_w(mT):
    w = np.zeros(7)
    j = int(np.clip(np.searchsorted(LNA, math.log(mT)), 1, 6))
    t = (math.log(mT) - LNA[j - 1]) / (LNA[j] - LNA[j - 1])
    w[j - 1], w[j] = 1 - t, t
    return w


def ix(*mts):
    return [list(AMPS).index(m) for m in mts]


# outputs: linear functionals of ln w at the 7 nominal amplitudes
OUT = {"beta_all": slope_w(range(7))}
for win in ((9, 20, 35), (35, 50, 70), (70, 100, 150), (20, 35, 50), (50, 70, 100), (9, 20, 35, 50), (50, 70, 100, 150),
            (9, 20, 35, 50, 70, 100), (20, 35, 50, 70, 100, 150)):
    w = np.zeros(7)
    w[ix(*win)] = slope_w(ix(*win))
    OUT["beta@%d-%dmT" % (win[0], win[-1]) + ("" if len(win) == 3 else "(%d)" % len(win))] = w
for mT in (9, 10, 20, 35, 50, 70, 100, 150):
    OUT["lnw@%gmT" % mT] = interp_w(mT)


def per_run_outputs(lnw):
    """per-run outputs; with a missing amplitude the slopes use the remaining amplitudes of the range (>= 2)."""
    ok = np.isfinite(lnw)
    out = {}
    for k, wv in OUT.items():
        if k.startswith("beta"):
            idx = [i for i in np.nonzero(wv)[0] if ok[i]]
            out[k] = float(np.polyfit(LNA[idx], lnw[idx], 1)[0]) if len(idx) >= 2 else np.nan
        else:
            nz = wv != 0
            out[k] = float(lnw[nz] @ wv[nz]) if np.all(ok[nz]) else np.nan
    return out
