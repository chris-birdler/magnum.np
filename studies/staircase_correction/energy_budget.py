#!/usr/bin/env python
"""
Design B of the dynamic test (PROTOKOLL 3.17): energy of each term and each particle in the start states of the four
arms, each state in the model of its own arm (built as run_loops builds it, double precision).

Energy of a linear term X in particle k: E_X,k = sum over the cells of k of −½ μ0 M_s,cell m · h_X V_cell, in units of
K_d V_box (as the term energies of run_loops). Terms: demag (with the near-field correction "demag_corr" of the arm
with the staircase correction), exchange, aniso_cube.

Geometry effects of each term, as Δ_c,k and Δ_g,k of 3.16 but as energy differences:
    δ_c,X,k = E_X,k(T3Pc) − E_X,k(T3)      the correction at dx 3
    δ_g,X,k = E_X,k(D15) − E_X,k(T15g)     the fine staircase at dx 1.5
The start states of the four arms are nearly the same state (3.15: overlap 0.98 … 0.9995), so these differences show
how each term of each model sees the geometry.

    python energy_budget.py --runs RUNS [--seeds 941 944] --out results   (writes energy_budget_B.json and .md)
RUNS/<arm>_s<seed>/init.pt for the arms T3, T3Pc, D15 and T15g.
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
sys.path.insert(0, str(HERE.parent / "beta_bpeak"))
os.environ.setdefault("CUDA_DEVICE", "-1")

from magnumnp import set_precision  # noqa: E402
set_precision("double")
torch.set_num_threads(int(os.environ.get("STAIRCASE_THREADS", os.cpu_count())))

import run_loops as RL  # noqa: E402
from anisotropy_model import box_edge  # noqa: E402
from geometry import fcc_box  # noqa: E402
from units import units  # noqa: E402

ARMS = ("T3", "T3Pc", "D15", "T15g")
CACHE = os.environ.get("STAIRCASE_KERNEL_CACHE")
MU0 = 4e-7 * np.pi


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def build(cfg, m):
    """State and terms of the arm (material and terms as run_loops.main, no drive). Returns (state, terms, ids, Kd)."""
    from magnumnp import Mesh, State, ExchangeField, UniaxialAnisotropyField, DemagField
    args = argparse.Namespace(**{**RL.NEW_KEYS, **cfg})
    U = units(args.Js, args.A)
    Ms, Kd, lex = U["Ms"], U["Kd"], U["l_ex"]
    a_lex, d_lex = box_edge(args.d_lex, args.phi, args.a_quant)
    N = int(round(a_lex / args.dx_lex))
    d, dx = d_lex * lex, args.dx_lex * lex
    rep = RL.geometry_rep(args)
    sub = (args.staircase_q if args.staircase != "none" else 8) if args.surface == "fraction" else 0
    geo = RL.locked_geometry(N, d, dx, rep, sub=sub)
    ids = geo["ids"]
    assert args.strayfield == "newell_ext" and args.r_p == 0.0
    stray = DemagField(true_pbc=True, extrapolate=True, cache_dir=CACHE)     # = RL.strayfield_setup, with a cache
    state = State(Mesh((N, N, N), (dx, dx, dx), pbc=(args.strayfield_shells,) * 3))
    T = lambda a: torch.tensor(a, dtype=state.dtype)  # noqa: E731
    mag4 = (ids >= 0)[..., None]
    frac4 = geo["frac"].astype(np.float64)[..., None]
    state.material = {"alpha": args.alpha}
    state.material["Ms"] = T(np.where(mag4, np.maximum(frac4 * Ms, args.ms_void * Ms), args.ms_void * Ms))
    state.material["A"] = T(np.where(mag4, frac4 * args.A, 0.0))
    terms = {"demag": stray}
    if args.staircase != "none":
        from magnumnp import DemagStaircaseCorrection
        from magnumnp.field_terms.demag_staircase import staircase_cells
        q_s = args.staircase_q
        n_s, pid_s, occ_s = staircase_cells(fcc_box(q_s * N, d, dx / q_s)["ids"], q_s)
        assert np.array_equal(pid_s, ids.astype(np.int64)) and np.allclose(n_s / q_s ** 3, geo["frac"], atol=1e-6)
        terms["demag_corr"] = DemagStaircaseCorrection(n_s, occ_s, q_s, args.staircase_p, Ms, form=args.staircase,
                                                       centroid=args.staircase_centroid, demag=stray)
    terms["exchange"] = ExchangeField()
    L = args.Leff_lex
    Q_eff = args.Q_eff if args.Q_eff is not None else 1.0 / L ** 2
    ax_c, _ = RL.locked_cube_axes(ids, rep, args.dx_lex, a_lex, d_lex, L, args.seed)
    state.material["Ku"] = T(np.where(mag4, frac4 * Q_eff * Kd, 0.0))
    state.material["Ku_axis"] = T(ax_c)
    terms["aniso_cube"] = UniaxialAnisotropyField()
    state.m = m.to(state.dtype)
    return state, terms, ids, Kd


@torch.no_grad()
def budget(path):
    st = torch.load(str(path), map_location="cpu", weights_only=False)
    cfg = json.loads((path.parent / "config.json").read_text())            # the full config of the run
    diff = [k for k, v in st["config"].items() if cfg.get(k, RL.NEW_KEYS.get(k)) != v]
    assert not diff, "config.json and init.pt of %s differ: %s" % (path.parent, diff)
    m = st["m"]
    state, terms, ids, Kd = build(cfg, m)
    N = ids.shape[0]
    vcell = float(np.prod(state.mesh.dx))
    norm = Kd * N ** 3 * vcell
    idt = torch.as_tensor(ids.astype(np.int64)).reshape(-1)
    out, tot = {}, {}
    for name, t in terms.items():
        h = t.h(state)
        e = (-0.5 * MU0 * state.material["Ms"] * (state.m * h).sum(-1, keepdim=True) * vcell).reshape(-1)
        out[name] = [float(e[idt == k].sum() / norm) for k in range(4)]
        tot[name] = float(t.E(state)) / norm
    ini = st.get("init", {})
    hist = ini.get("relax_history") or []
    # record: the energy of the prepared state (after the AC preparation if there was one, else the relaxation)
    rec = (ini.get("predemag") or {}).get("E_over_Kd") or (hist[-1]["E_over_Kd"] if hist else None)
    return {"per_particle": out, "total_E": tot, "sum_total": sum(tot.values()), "record_E": rec,
            "void_share": {k: tot[k] - sum(out[k]) for k in out}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True)
    ap.add_argument("--seeds", type=int, nargs="+", default=[941, 942, 943, 944])
    ap.add_argument("--out", default=str(HERE / "results"))
    a = ap.parse_args()
    root, out = pathlib.Path(a.runs), pathlib.Path(a.out)
    fj = out / "energy_budget_B.json"
    res = json.loads(fj.read_text()) if fj.exists() else {}
    for s in a.seeds:
        for arm in ARMS:
            key = "%s_s%d" % (arm, s)
            if key in res:
                continue
            t = time.time()
            res[key] = budget(root / key / "init.pt")
            log(key, "sum %.9f record %s (%.0f s)" % (res[key]["sum_total"], res[key]["record_E"], time.time() - t))
            fj.write_text(json.dumps(res, indent=1))
    lines = ["# Design B: energy of each term and particle in the start states (PROTOKOLL 3.17)\n",
             "Units: 10⁻⁶ K_d V_box. δ_c = E(T3Pc) − E(T3), δ_g = E(D15) − E(T15g); demag of T3Pc includes demag_corr.\n",
             "| seed | term | δ_c per particle 0 / 1 / 2 / 3 | δ_g per particle 0 / 1 / 2 / 3 | δ_c − δ_g |",
             "|---|---|---|---|---|"]
    seeds = sorted({int(k.split("_s")[1]) for k in res})
    for s in seeds:
        if not all("%s_s%d" % (arm, s) in res for arm in ARMS):
            continue
        E = {arm: res["%s_s%d" % (arm, s)]["per_particle"] for arm in ARMS}
        for X in ("demag", "exchange", "aniso_cube"):
            g = lambda arm: np.array(E[arm][X]) + (np.array(E[arm]["demag_corr"]) if X == "demag" and "demag_corr" in E[arm] else 0)  # noqa: E731
            dc, dg = g("T3Pc") - g("T3"), g("D15") - g("T15g")
            f = lambda v: " / ".join("%+.1f" % (1e6 * x) for x in v)  # noqa: E731
            lines.append("| %d | %s | %s | %s | %s |" % (s, X, f(dc), f(dg), f(dc - dg)))
    lines += ["", "Check: the sum of the term energies against the last relaxation record of the run (single precision "
              "in the run):", ""]
    for k, v in sorted(res.items()):
        lines.append("- %s: sum %.9f, record %s" % (k, v["sum_total"], v["record_E"]))
    (out / "energy_budget_B.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
