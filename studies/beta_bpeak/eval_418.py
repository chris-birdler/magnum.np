#!/usr/bin/env python
"""PLAN 4.18 evaluation: tilt of the cut cells, large-angle pairs, overlap, field ratio rho, and the readings and the
verdict of step 1.

Terms (PLAN 4.18): a cut cell is a magnetic cell with 0 < f < 1. The tilt of a cell is the angle between its m and the
sum of m over the full cells (f = 1) of the same particle within 2 cells (run_loops.full_cell_sum). A large-angle pair
is a pair of neighbour cells of one particle with an angle above 90 deg. rho_i = |h_d,perp,i| / sum_j (l_ex/dx)^2
2 f_j / (f_i + f_j) in the aligned state (diag_fields.npz of run_loops --diag_fields).

    eval_418.py table <run dir>[:<state file>] ...
        tilt table of relaxed states (default init.pt)
    eval_418.py step1 <runs dir> --ref <dir with T3Pc8_s941, T3Pc8_s944> --d15 <dir with D15_s<seed>/init.pt>
                [--out results/plan418_step1]
        readings 1.0, 1a ... 1d, the verdict and k* (arms B8, C8t, C8h, T4, T16 in <runs dir>/<arm>_s<seed>)
"""
import argparse
import json
import math
import pathlib
import sys
import types

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

SEEDS = (941, 942, 943, 944)
ARMS = ("B8", "C8t", "C8h", "T4", "T16")
REF_SEEDS = (941, 944)                    # GPU states T3Pc8 (staircase box 6) for the method check 1.0
SLOPE_FIT, SLOPE_CHECK = (941, 942), (943, 944)


def geometry(cfg):
    """(ids, f) of the mesh of a run (run_loops.locked_geometry with the sub-cells of run_loops.surface_sub)."""
    from run_loops import locked_geometry, surface_sub, NEW_KEYS
    ns = types.SimpleNamespace(surface=cfg.get("surface", "voxel"), staircase=cfg.get("staircase", "none"),
                               staircase_q=cfg.get("staircase_q", NEW_KEYS["staircase_q"]))
    g = locked_geometry(cfg["N_used"], cfg["d"], cfg["dx"], 1, sub=surface_sub(ns))
    return g["ids"].astype(np.int64), g["frac"].astype(np.float64)


def load_run(run, state="init.pt"):
    """(config, m [N, N, N, 3] float64, the init record of the state file or {})."""
    import torch
    run = pathlib.Path(run)
    cfg = json.loads((run / "config.json").read_text())
    c = torch.load(str(run / state), map_location="cpu", weights_only=False)
    return cfg, c["m"].double().numpy(), c.get("init", {})


def tilt_deg(m, ids, frac):
    """(tilt of each cell in deg, mask of the cells with at least one full cell in the block)."""
    from run_loops import full_cell_sum
    s = full_cell_sum(m, ids, frac)
    n = np.linalg.norm(s, axis=-1)
    c = (m * s).sum(-1) / np.maximum(n, 1e-30)
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0))), n > 0


def pair_angles(m, ids, keep=None):
    """Angles (deg) of the neighbour pairs of one particle; keep: optional cell mask, a pair counts if both cells are
    in it."""
    out = []
    for ax in range(3):
        same = (ids >= 0) & (ids == np.roll(ids, 1, axis=ax))
        if keep is not None:
            same &= keep & np.roll(keep, 1, axis=ax)
        c = (m * np.roll(m, 1, axis=ax)).sum(-1)[same]
        out.append(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))
    return np.concatenate(out)


def state_readings(m, ids, frac, keep=None):
    """Tilt readings of the cut cells (largest, 99th percentile, number above 45 and 90 deg), the largest tilt of the
    full cells, and the large-angle pairs (number, largest angle). keep: optional mask of the cells to evaluate."""
    th, ref = tilt_deg(m, ids, frac)
    mag = ids >= 0
    if keep is None:
        keep = np.ones(ids.shape, bool)
    cut = mag & (frac < 1.0) & ref & keep
    full = mag & (frac >= 1.0) & ref & keep
    a = th[cut]
    pa = pair_angles(m, ids, keep)
    return {"n_cut": int(cut.sum()), "cut_max": float(a.max()), "cut_p99": float(np.percentile(a, 99)),
            "cut_gt45": int((a > 45.0).sum()), "cut_gt90": int((a > 90.0).sum()),
            "full_max": float(th[full].max()), "pairs_gt90": int((pa > 90.0).sum()),
            "pair_max": float(pa.max()) if pa.size else None}


def overlap(m, m_ref, ids):
    """Mean of m . m_ref over the cells of each particle."""
    d = (m * m_ref).sum(-1)
    return [float(d[ids == p].mean()) for p in range(int(ids.max()) + 1)]


def rho_theta(npz):
    """(theta in rad, rho) of the cut cells with a full cell in the block, from diag_fields.npz: theta from the relaxed m,
    rho from the stray-field terms (names that start with "demag") in the aligned state."""
    z = np.load(npz)
    ids, f = z["ids"].astype(np.int64), z["frac"].astype(np.float64)
    names = [str(x) for x in z["term_names"]]
    hd = sum(z["h_aligned_%s" % nm].astype(np.float64) for nm in names if nm.startswith("demag"))
    u = z["m_aligned"].astype(np.float64)
    hp = np.linalg.norm(hd - (hd * u).sum(-1, keepdims=True) * u, axis=-1) / float(z["Ms"])
    c2 = (float(z["l_ex"]) / float(z["dx"])) ** 2
    den = np.zeros(ids.shape)
    for ax in range(3):
        for sh in (1, -1):
            fj = np.roll(f, sh, axis=ax)
            same = (ids >= 0) & (ids == np.roll(ids, sh, axis=ax))
            den += np.where(same, 2.0 * fj / np.maximum(f + fj, 1e-30), 0.0)
    th, ref = tilt_deg(z["m_relaxed"].astype(np.float64), ids, f)
    cut = (ids >= 0) & (f < 1.0) & ref & (den > 0)
    return np.radians(th[cut]), hp[cut] / (c2 * den[cut])


def spearman(x, y):
    from scipy.stats import spearmanr
    return float(spearmanr(x, y).statistic)


def start_state(cfg, d15_init):
    """The transfer of the D15 state onto the mesh of the run (the start of B8 and C8t) and the mask of the cells with
    a voxel source (run_loops.transfer_m: a coarse cell without a fine sub-cell of its own particle is left out)."""
    import torch
    from run_loops import locked_geometry, transfer_m
    ids, _ = geometry(cfg)
    cf = torch.load(str(d15_init), map_location="cpu", weights_only=False)
    dxf = float(cf["config"]["dx_lex"])
    Nf = int(round(cfg["a_lex"] / dxf))
    ids_f = locked_geometry(Nf, cfg["d"], cfg["dx"] * dxf / cfg["dx_lex"], 1)["ids"]
    m_c, st = transfer_m(cf["m"].numpy(), ids_f, ids)
    N, k = ids.shape[0], Nf // ids.shape[0]
    same = ids_f.reshape(N, k, N, k, N, k) == ids[:, None, :, None, :, None]
    return m_c, same.sum(axis=(1, 3, 5)) > 0, st


def arm_values(runs, arm, seeds=SEEDS):
    """{seed: (config, m, init record)} of the arm runs with a converged relaxation; the others are listed in skipped."""
    vals, skipped = {}, {}
    for s in seeds:
        run = pathlib.Path(runs) / ("%s_s%d" % (arm, s))
        if not (run / "init.pt").exists():
            skipped[s] = "no init.pt"
            continue
        cfg, m, rec = load_run(run)
        if not rec.get("relax_converged"):
            skipped[s] = "not converged after %s iterations" % rec.get("relax_iterations")
            continue
        vals[s] = (cfg, m, rec)
    return vals, skipped


def method_check(b8, ref_runs):
    """1.0: B8 against the GPU states T3Pc8 (same seed): overlap at least 0.99 in each particle, largest cut-cell tilt
    and its 99th percentile within 5 deg, number of cut cells above 45 deg within 2, large-angle pairs within 1."""
    out = {}
    for s in REF_SEEDS:
        ref = pathlib.Path(ref_runs) / ("T3Pc8_s%d" % s)
        if s not in b8 or not (ref / "init.pt").exists():
            out[s] = {"pass": False, "reason": "no B8 value or no GPU state"}
            continue
        cfg, m, _ = b8[s]
        _, m_r, _ = load_run(ref)
        ids, f = geometry(cfg)
        a, r = state_readings(m, ids, f), state_readings(m_r, ids, f)
        ov = overlap(m, m_r, ids)
        ok = {"overlap": min(ov) >= 0.99, "cut_max": abs(a["cut_max"] - r["cut_max"]) <= 5.0,
              "cut_p99": abs(a["cut_p99"] - r["cut_p99"]) <= 5.0, "cut_gt45": abs(a["cut_gt45"] - r["cut_gt45"]) <= 2,
              "pairs_gt90": abs(a["pairs_gt90"] - r["pairs_gt90"]) <= 1}
        out[s] = {"pass": all(ok.values()), "checks": ok, "overlap": ov, "cpu": a, "gpu": r}
    return {"pass": all(v["pass"] for v in out.values()), "seeds": out}


def predictions(R):
    """1a and 1b per seed from the readings R[arm][seed]. A seed takes part in 1a if B8 has at least 3 cut cells above
    45 deg. 1a: C8t has at most 1 cut cell above 45 deg and no large-angle pair. T4: at most half the number of B8
    (same seed) and no pair. T16: at most 1 cut cell above 45 deg and no pair. None: no value of the arm."""
    part = [s for s in SEEDS if s in R["B8"] and R["B8"][s]["cut_gt45"] >= 3]
    h1a = {s: (R["C8t"][s]["cut_gt45"] <= 1 and R["C8t"][s]["pairs_gt90"] == 0) if s in R["C8t"] else None
           for s in part}
    t4 = {s: (R["T4"][s]["cut_gt45"] <= 0.5 * R["B8"][s]["cut_gt45"] and R["T4"][s]["pairs_gt90"] == 0)
          for s in SEEDS if s in R["T4"] and s in R["B8"]}
    t16 = {s: (R["T16"][s]["cut_gt45"] <= 1 and R["T16"][s]["pairs_gt90"] == 0) for s in SEEDS if s in R["T16"]}
    return {"part": part, "1a": h1a, "T4": t4, "T16": t16}


def verdict(R, method_pass):
    """Verdict of step 1 (PLAN 4.18): "confirmed" if 1.0 passes, at least 3 seeds take part, 1a holds in each of them
    and T16 meets its prediction in at least 3 seeds; "refuted" if C8t or T16 keeps at least half the number of cut
    cells above 45 deg of B8 (same seed) in at least 2 of the seeds that take part; else "not determined". An arm with
    fewer than 3 values gives "not determined" (never "confirmed"). k* = 4 if T4 meets its prediction in at least 3
    seeds, else 16."""
    p1 = predictions(R)
    part = p1["part"]
    keeps = {arm: sum(1 for s in part if s in R[arm] and R[arm][s]["cut_gt45"] >= 0.5 * R["B8"][s]["cut_gt45"])
             for arm in ("C8t", "T16")}
    few = [a for a in ("B8", "C8t", "T16") if len(R[a]) < 3]
    n_t16 = sum(1 for v in p1["T16"].values() if v)
    if keeps["C8t"] >= 2 or keeps["T16"] >= 2:
        v = "refuted"
    elif method_pass and not few and len(part) >= 3 and all(x is True for x in p1["1a"].values()) and n_t16 >= 3:
        v = "confirmed"
    else:
        v = "not determined"
    return {"step1": v, "participating": part, "C8t_keeps_half": keeps["C8t"], "T16_keeps_half": keeps["T16"],
            "T16_holds": n_t16, "arms_with_fewer_than_3_values": few,
            "k_star": 4 if sum(1 for x in p1["T4"].values() if x) >= 3 else 16}


def step1(runs, ref_runs, d15_dir):
    arms, skipped = {}, {}
    for a in ARMS:
        arms[a], skipped[a] = arm_values(runs, a)
    R = {a: {} for a in ARMS}
    geo = {}
    for a in ARMS:
        for s, (cfg, m, rec) in arms[a].items():
            if s not in geo:
                geo[s] = geometry(cfg)
            ids, f = geo[s]
            R[a][s] = dict(state_readings(m, ids, f), relax_iterations=rec.get("relax_iterations"),
                           E_over_Kd=rec.get("E_over_Kd"))
            if a != "B8" and s in arms["B8"]:
                R[a][s]["overlap_vs_B8"] = overlap(m, arms["B8"][s][1], ids)
    res = {"readings": R, "skipped": skipped}
    res["1.0"] = method_check(arms["B8"], ref_runs)

    p1 = predictions(R)
    res["1a"] = {"participating": p1["part"], "holds": p1["1a"], "C8h_info": {s: R["C8h"].get(s) for s in p1["part"]}}
    res["1b"] = {"T4": p1["T4"], "T16": p1["T16"]}

    # 1c: formation (information): the start state (transfer, cells without a voxel source left out), B8 and C8t
    res["1c"] = {}
    for s in SEEDS:
        src = pathlib.Path(d15_dir) / ("D15_s%d" % s) / "init.pt"
        cfg = (arms["B8"].get(s) or arms["C8t"].get(s) or (None,))[0]
        if cfg is None or not src.exists():
            continue
        m0, has_src, st = start_state(cfg, src)
        ids, f = geo[s]
        res["1c"][s] = {"start": state_readings(m0, ids, f, keep=has_src), "transfer": st,
                        "B8": R["B8"].get(s), "C8t": R["C8t"].get(s)}

    # 1d: field ratio (support)
    rt = {}
    for s in SEEDS:
        npz = pathlib.Path(runs) / ("B8_s%d" % s) / "diag_fields.npz"
        if s in arms["B8"] and npz.exists():
            rt[s] = rho_theta(npz)
    d1 = {"spearman": {s: spearman(*rt[s]) for s in rt}}
    fit = [s for s in SLOPE_FIT if s in rt]
    if fit:
        th = np.concatenate([rt[s][0] for s in fit])
        rh = np.concatenate([rt[s][1] for s in fit])
        sel = rh >= 0.05
        a_sl = float(np.median(th[sel] / rh[sel])) if sel.any() else None
        d1["slope_a"] = a_sl
        d1["check"] = {}
        for s in SLOPE_CHECK:
            if s in rt and a_sl:
                t, r = rt[s]
                q = float(np.median(t[r >= 0.05] / (a_sl * r[r >= 0.05])))
                d1["check"][s] = {"median_ratio": q, "in_0.67_1.5": 0.67 <= q <= 1.5}
    d1["prediction_rank_ge_0.5"] = {s: v >= 0.5 for s, v in d1["spearman"].items()}
    res["1d"] = d1

    res["verdict"] = verdict(R, res["1.0"]["pass"])
    return res


def fmt_table(rows):
    lines = ["| run | cut cells | cut-cell tilt: largest | 99th percentile | cut cells > 45° | large-angle pairs | "
             "full cells: largest |", "|---|---|---|---|---|---|---|"]
    for name, r in rows:
        lines.append("| %s | %d | %.1f° | %.1f° | %d | %d | %.1f° |" % (
            name, r["n_cut"], r["cut_max"], r["cut_p99"], r["cut_gt45"], r["pairs_gt90"], r["full_max"]))
    return "\n".join(lines)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("table")
    t.add_argument("runs", nargs="+")
    s1 = sub.add_parser("step1")
    s1.add_argument("runs")
    s1.add_argument("--ref", required=True)
    s1.add_argument("--d15", required=True)
    s1.add_argument("--out", default=None)
    a = p.parse_args(argv)
    if a.cmd == "table":
        rows = []
        for r in a.runs:
            run, _, state = r.partition(":")
            cfg, m, _ = load_run(run, state or "init.pt")
            ids, f = geometry(cfg)
            rows.append((pathlib.Path(run).name, state_readings(m, ids, f)))
        print(fmt_table(rows))
        return
    res = step1(a.runs, a.ref, a.d15)
    txt = json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    if a.out:
        out = pathlib.Path(a.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.with_suffix(".json").write_text(txt)
        rows = [("%s_s%d" % (arm, s), r) for arm in ARMS for s, r in sorted(res["readings"][arm].items())]
        out.with_suffix(".md").write_text(fmt_table(rows) + "\n\nVerdict of step 1: %s; k* = %d\n" % (
            res["verdict"]["step1"], res["verdict"]["k_star"]))
    print(txt)


if __name__ == "__main__":
    main()
