"""
Stage A of PLAN.md: static test of the near-field correction (V1 cells after Donahue & McMichael 2007) on the
production geometry, with the clean prepared states of the Newell solver (studies/beta_bpeak/runs/states_pd300).

Usage: python check_static.py --runs K4_1_nx_pd300 [K2_0_nx_pd300 K2_1_nx_pd300 ...] --q 2 --p 1 2 3 --out results
All runs of one call must have the same geometry. Writes <out>/static_A_<geometry>.json.
"""
import argparse
import json
import os
import pathlib
import re
import resource
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

import run_loops as RL  # noqa: E402
from anisotropy_model import box_edge, particle_axes  # noqa: E402
from units import units  # noqa: E402
from geometry import fcc_box  # noqa: E402
from nearfield import NearCorrection, FFTField, v1_cells, fine_owner_m, coarse_mean  # noqa: E402

STATES = pathlib.Path(os.environ.get("STAIRCASE_STATES",
                                     "/home/christoph/Documents/code/magnum.np-beta_bpeak/studies/beta_bpeak/runs"))
CACHE = os.environ.get("STAIRCASE_KERNEL_CACHE")   # optional folder for the demag kernels (local raw data)
LAYER = 3                       # PLAN: voxel cells within 3 cells of a modified cell, the same for all p


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def load_run(run):
    cfg = json.loads((STATES / run / "config.json").read_text())
    rec = None
    for line in (STATES / run / "stdout.log").read_text().splitlines():
        mm = re.search(r"\[relax2\].*E/K_d\s+(-?[0-9.eE+-]+)", line)
        if mm:
            rec = float(mm.group(1))
    st = torch.load(STATES / "states_pd300" / run / "init.pt", map_location="cpu", weights_only=False)
    return cfg, rec, st["m"].to(torch.float64)


def production(cfg, m, stray):
    """Material and terms as run_loops.main (no drive); returns (energies / (K_d V_box), ids, N, d_lex, dx_lex, a_lex)."""
    from magnumnp import Mesh, State, ExchangeField, UniaxialAnisotropyField
    args = argparse.Namespace(**cfg)
    U = units(args.Js, args.A)
    Ms, Kd, lex = U["Ms"], U["Kd"], U["l_ex"]
    a_lex, d_lex = box_edge(args.d_lex, args.phi, args.a_quant)
    N = int(round(a_lex / args.dx_lex))
    d, dx = d_lex * lex, args.dx_lex * lex
    rep = RL.geometry_rep(args)
    geo = RL.locked_geometry(N, d, dx, rep, sub=0)
    ids = geo["ids"]
    pbc = (args.strayfield_shells,) * 3
    state = State(Mesh((N, N, N), (dx, dx, dx), pbc=pbc))
    T = lambda a: torch.tensor(a, dtype=state.dtype)  # noqa: E731
    mag4 = (ids >= 0)[..., None]
    frac4 = geo["frac"].astype(np.float64)[..., None]
    state.material = {"alpha": args.alpha}
    state.material["Ms"] = T(np.where(mag4, np.maximum(frac4 * Ms, args.ms_void * Ms), args.ms_void * Ms))
    state.material["A"] = T(np.where(mag4, frac4 * args.A, 0.0))
    e_d = np.asarray(args.direction, float)
    e_d /= np.linalg.norm(e_d)
    L = args.Leff_lex
    Q_eff = args.Q_eff if args.Q_eff is not None else 1.0 / L ** 2
    K_eff = Q_eff * Kd
    K_p = args.r_p * K_eff
    terms = {"demag": stray, "exchange": ExchangeField()}
    if K_eff > 0.0:
        ax_c, _ = RL.locked_cube_axes(ids, rep, args.dx_lex, a_lex, d_lex, L, args.seed)
        state.material["Ku"] = T(np.where(mag4, frac4 * K_eff, 0.0))
        state.material["Ku_axis"] = T(ax_c)
        terms["aniso_cube"] = UniaxialAnisotropyField()
    ax_p = particle_axes(e_d)
    if K_p > 0.0:
        kp_ax = np.zeros((N, N, N, 3))
        for p in range(4):
            kp_ax[ids == p] = ax_p[p]
        state.material["Kp"] = T(np.where(mag4, frac4 * K_p, 0.0))
        state.material["Kp_axis"] = T(kp_ax)
        terms["aniso_particle"] = UniaxialAnisotropyField(Ku="Kp", Ku_axis="Kp_axis")
    state.m = m.to(state.dtype)
    with torch.no_grad():
        E = {k: float(t.E(state)) / (Kd * N ** 3 * dx ** 3) for k, t in terms.items()}
    return E, ids, N, d_lex, args.dx_lex, a_lex


def bloch_count(m, ids):
    q = RL.cube_charge_t(m)
    sel = RL.bloch_mask_t(torch.as_tensor(ids.astype(np.int64))) & (q.abs() > 0.5)
    return int((sel & (q > 0)).sum()), int((sel & (q < 0)).sum())


def nearest_voxel_map(ids_vox, pid, n):
    """For each magnetic V1 cell that is not a voxel cell: flat index of the nearest voxel cell of the same particle
    (periodic centre distance; ties by the lower cell index)."""
    from scipy.spatial import cKDTree
    N = ids_vox.shape[0]
    new = (n > 0) & (ids_vox < 0)
    src = np.flatnonzero(new)
    dst = np.empty_like(src)
    xyz = np.stack(np.unravel_index(np.arange(N ** 3), (N,) * 3), -1).astype(np.float64)
    for p in range(4):
        vp = np.flatnonzero(ids_vox.reshape(-1) == p)
        sp = src[pid.reshape(-1)[src] == p]
        if len(sp) == 0:
            continue
        tree = cKDTree(xyz[vp], boxsize=N)
        dist, k = tree.query(xyz[sp], k=1)
        dst[np.searchsorted(src, sp)] = vp[k]
    return src, dst


def energy_coarse(w, m, h):
    N = m.shape[0]
    return float(-(torch.as_tensor(w, dtype=torch.float64).reshape(N, N, N, 1) * m * h).sum() / N ** 3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--q", type=int, default=2)
    ap.add_argument("--p", type=int, nargs="+", default=[1, 2, 3])
    ap.add_argument("--out", default=str(HERE / "results"))
    a = ap.parse_args()
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    q = a.q
    res = {"plan": "PLAN.md stage A", "q": q, "p": a.p, "layer": LAYER, "runs": {}, "times_s": {}}
    t0 = time.time()

    runs = [(r,) + load_run(r) for r in a.runs]
    cfg0 = runs[0][1]
    shells = cfg0["strayfield_shells"]
    log("runs", a.runs, "strayfield", cfg0["strayfield"], shells)
    assert all(r[1]["strayfield"] == "newell_ext" for r in runs), "only states of the Newell solver (rule of Chris)"

    # coarse field (production setting); its kernel is shared with the energy check of step 2
    coarse = None
    for run, cfg, rec, m in runs:
        from magnumnp import DemagField
        if coarse is None:
            stray = DemagField(true_pbc=True, extrapolate=True, cache_dir=CACHE)
        t = time.time()
        E, ids_vox, N, d_lex, dx_lex, a_lex = production(cfg, m, stray)
        if coarse is None:
            coarse = FFTField(N, ("newell_ext", shells))
            coarse.term = stray                                   # same kernel (scale-free)
            res["times_s"]["coarse_kernel_and_check"] = time.time() - t
        norm_dev = float((m.norm(dim=-1)[torch.as_tensor(ids_vox >= 0)] - 1).abs().max())
        bp = bloch_count(m.to(torch.float32), ids_vox)
        Etot = sum(E.values())
        ok = (rec is not None and abs(Etot - rec) <= 1e-6 and norm_dev < 1e-5 and bp == (0, 0)
              and tuple(m.shape[:3]) == (N, N, N))
        res["runs"][run] = {"record_E": rec, "E_terms": E, "E_total": Etot, "E_K": E.get("aniso_cube", 0.0)
                            + E.get("aniso_particle", 0.0), "norm_dev": norm_dev, "bloch": bp, "check_ok": ok}
        log(run, "step 2: E_total %.10f record %s |m|-1 %.1e BP %s -> %s" % (Etot, rec, norm_dev, bp, ok))
    geo_key = "d%.0f_phi%.2f" % (cfg0["d_lex"], cfg0["phi"])
    res["geometry"] = {"N": N, "d_lex": d_lex, "dx_lex": dx_lex, "a_lex": a_lex, "key": geo_key}

    # step 3: fine body and V1 cells (R0)
    ids_f = fcc_box(q * N, d_lex, dx_lex / q)["ids"]
    body = ids_f >= 0
    n, pid, occ = v1_cells(ids_f, q)
    vox = ids_vox >= 0
    mag = n > 0
    nq = q ** 3
    r0 = {"same_particle": bool(np.all(pid[vox] == ids_vox[vox])), "voxel_v0": int((vox & ~mag).sum()),
          "n_voxel": int(vox.sum()), "n_magnetic": int(mag.sum()), "n_new": int((mag & ~vox).sum()),
          "n_modified": int((mag & (n < nq)).sum()), "v_hist": np.bincount(n[mag], minlength=nq + 1).tolist(),
          "phi_vox": float(vox.mean()), "phi_f": float(body.mean())}
    r0["ok"] = r0["same_particle"] and r0["voxel_v0"] == 0
    res["R0"] = r0
    log("R0", {k: r0[k] for k in r0 if k != "v_hist"})
    if not r0["ok"]:
        (out / ("static_A_%s.json" % geo_key)).write_text(json.dumps(res, indent=1))
        sys.exit("R0 failed: stop")
    src, dst = nearest_voxel_map(ids_vox, pid, n)
    v = n / nq
    mod = mag & (n < nq)
    from scipy.ndimage import maximum_filter
    layer = maximum_filter(mod, size=2 * LAYER + 1, mode="wrap") & vox
    newc = mag & ~vox
    lay_t, new_t = torch.as_tensor(layer), torch.as_tensor(newc)
    fine = FFTField(q * N, ("newell_ext", shells), cache_dir=CACHE)
    body_t = torch.as_tensor(body)

    def v1_m(m):
        mv = m.clone().reshape(-1, 3)
        mv[torch.as_tensor(src)] = mv[torch.as_tensor(dst)]
        mv = mv.reshape(N, N, N, 3)
        return torch.where(torch.as_tensor(mag)[..., None], mv, torch.zeros_like(mv))

    def reference(mv):
        t = time.time()
        mf = fine_owner_m(mv, q)
        hf = fine.h(body.astype(np.float64), mf)
        E = float(-(body_t[..., None] * mf * hf).sum() / (q * N) ** 3)
        h = coarse_mean(hf, body, q)
        del mf, hf
        res["times_s"].setdefault("reference", []).append(time.time() - t)
        return E, h

    # fields of the uncorrected models and of the reference
    fields = {}
    tests = [("uniform_100", (1, 0, 0)), ("uniform_110", (1, 1, 0)), ("uniform_111", (1, 1, 1))]
    for name, e in tests:
        e = torch.tensor(e, dtype=torch.float64) / np.linalg.norm(e)
        m = torch.zeros(N, N, N, 3, dtype=torch.float64)
        m[...] = e
        fields[name] = m
    for run, cfg, rec, m in runs:
        if res["runs"][run]["check_ok"]:
            fields[run] = m
        else:
            log(run, "not used: a state check failed")
    stat = {}
    for name, m in fields.items():
        mv = v1_m(m)
        mvox = torch.where(torch.as_tensor(vox)[..., None], m, torch.zeros_like(m))
        h_vox = coarse.h(vox.astype(np.float64), mvox)
        h_frac = coarse.h(v, mv)
        E_ref, h_ref = reference(mv)
        stat[name] = {"mv": mv, "h_vox": h_vox, "h_frac": h_frac, "h_ref": h_ref}
        d = {"E_ref": E_ref, "E_vox": energy_coarse(vox, mvox, h_vox), "E_frac": energy_coarse(v, mv, h_frac)}
        d["S"] = d["E_vox"] - d["E_ref"]
        dv = (h_vox - h_ref)[lay_t]
        d["rms_layer_vox"] = float(dv.pow(2).sum(-1).mean().sqrt())
        d["rms_layer_frac"] = float((h_frac - h_ref)[lay_t].pow(2).sum(-1).mean().sqrt())
        d["rms_new_frac"] = float((h_frac - h_ref)[new_t].pow(2).sum(-1).mean().sqrt())
        res.setdefault("fields", {})[name] = d
        log(name, {k: d[k] for k in ("E_ref", "E_vox", "E_frac", "S", "rms_layer_vox")}, "rss %.1f GB" % rss_gb())
    phi_f, phi_vox = r0["phi_f"], r0["phi_vox"]
    r1 = {}
    for name, _ in tests:
        d = res["fields"][name]
        r1[name] = {"ref_rel": d["E_ref"] / (phi_f * (1 - phi_f) / 3) - 1, "vox_rel": d["E_vox"] / (phi_vox * (1 - phi_vox) / 3) - 1}
    res["R1"] = {"cases": r1, "ok": all(abs(c["ref_rel"]) <= 5e-4 and abs(c["vox_rel"]) <= 5e-4 for c in r1.values())}
    log("R1", res["R1"])
    if not res["R1"]["ok"]:
        (out / ("static_A_%s.json" % geo_key)).write_text(json.dumps(res, indent=1))
        sys.exit("R1 failed: computation not valid, stop")

    # corrected models, window by window
    for p in a.p:
        t = time.time()
        nc = NearCorrection(n, occ, q, p, periodic=True)
        tb = time.time() - t
        res["times_s"]["build_p%d" % p] = tb
        log("p", p, "modified cells", len(nc.J), "build %.0f s" % tb, "rss %.1f GB" % rss_gb())
        for name, s in stat.items():
            t = time.time()
            hP = s["h_frac"] + nc.h_pair(s["mv"])
            tP = time.time() - t
            hL = s["h_frac"] + nc.h_local(s["mv"])
            d = res["fields"][name]
            for form, h in (("P", hP), ("L", hL)):
                E = nc.energy(s["mv"], h)
                d["E_%s%d" % (form, p)] = E
                d["rms_layer_%s%d" % (form, p)] = float((h - s["h_ref"])[lay_t].pow(2).sum(-1).mean().sqrt())
                d["rms_new_%s%d" % (form, p)] = float((h - s["h_ref"])[new_t].pow(2).sum(-1).mean().sqrt())
            res["times_s"]["field_P_p%d" % p] = tP
            log(name, "p", p, "E_P %.8f E_L %.8f E_ref %.8f" % (d["E_P%d" % p], d["E_L%d" % p], d["E_ref"]))
        del nc
        (out / ("static_A_%s.json" % geo_key)).write_text(json.dumps(res, indent=1))
    res["times_s"]["total"] = time.time() - t0
    res["max_rss_gb"] = rss_gb()
    (out / ("static_A_%s.json" % geo_key)).write_text(json.dumps(res, indent=1))
    log("done", out / ("static_A_%s.json" % geo_key))


if __name__ == "__main__":
    main()
