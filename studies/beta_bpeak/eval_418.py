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
    eval_418.py step3 <runs dir> [--slope_a A] [--ref416 results/hf416.json] [--no_fields] [--out results/plan418_step3]
        readings 3a ... 3g and the verdict of step 3 (arms C, D, Dk in <runs dir>/S3<arm>_s<seed>)
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


S3_SEEDS = (921, 923, 925)               # 924 replaced by 925 (its preparation has a Bloch point)
S3_ARMS = ("C", "D", "Dk")
THIN = 1.0 / 3.0                          # thin cell: f below l_ex/dx (dx = 3 l_ex; amendment 1)


def read_csv_rows(path):
    """The rows of a run_loops csv file as lists of floats (comment lines left out); [] if the file is missing."""
    p = pathlib.Path(path)
    if not p.exists():
        return []
    return [[float(x) for x in ln.split(",")] for ln in p.read_text().splitlines()
            if ln.strip() and not ln.startswith("#")]


def load_step3_run(run):
    """Step-3 data of one run (stage amp00), or None without a DONE file: the summary stage, the Bloch counts per sample
    (cycle, sample, n_pos, n_neg), the Bloch rows (cycle, sample, i, j, k, Q, particle, largest cube angle, f_edge_lo,
    f_edge_hi, f_cube_min) and the large-angle pairs per sample (cycle, sample, n_pairs_gt90, max_angle)."""
    run = pathlib.Path(run)
    if not (run / "DONE").exists():
        return None
    s = json.loads((run / "summary.json").read_text())
    return {"init": s.get("init", {}), "stage": s["stages"][0], "bloch": read_csv_rows(run / "bloch_00_amp00.csv"),
            "rows": read_csv_rows(run / "bloch_pos_00_amp00.csv"), "n90": read_csv_rows(run / "diag_00_amp00.csv")}


def reading_3a(D):
    """3a on D {seed: run data}: the samples with a Bloch point in the cycles >= 1. A sample counts as "cut" if the cube
    of each of its Bloch points has a cut cell (smallest f of the 8 corners below 1). Holds if the share is at least
    90 %, pooled and in each seed with at least 5 such samples; None without such samples. Also the share of the
    samples with a large-angle pair (an edge above 90 deg) in a cube of the sample."""
    per = {}
    for s, r in D.items():
        smp = {}
        for x in r["rows"]:
            if x[0] >= 1:
                smp.setdefault((int(x[0]), int(x[1])), []).append(x)
        cut = sum(1 for v in smp.values() if all(x[10] < 1.0 for x in v))
        large = sum(1 for v in smp.values() if any(x[7] > 90.0 for x in v))
        per[s] = {"samples": len(smp), "samples_bloch_csv": sum(1 for b in r["bloch"] if b[0] >= 1 and b[2] + b[3] > 0),
                  "cut": cut, "large_angle": large, "share_cut": cut / len(smp) if smp else None,
                  "share_large": large / len(smp) if smp else None}
    n = sum(v["samples"] for v in per.values())
    big = [s for s, v in per.items() if v["samples"] >= 5]
    pooled = sum(v["cut"] for v in per.values()) / n if n else None
    holds = None if not n else (pooled >= 0.9 and all(per[s]["share_cut"] >= 0.9 for s in big))
    return {"seeds": per, "samples": n, "seeds_ge5": big, "pooled_share_cut": pooled,
            "pooled_share_large": sum(v["large_angle"] for v in per.values()) / n if n else None, "holds": holds}


def has_bloch(r, c0):
    return any(b[0] >= c0 and b[2] + b[3] > 0 for b in r["bloch"])


def reading_3b(D, Dk):
    """3b: in each seed in which D has a Bloch point in the cycles >= 2, Dk has none in its cycles >= 2. Holds (True or
    False) only with at least 2 such seeds with a Dk value, else None."""
    per = {s: (not has_bloch(Dk[s], 2) if s in Dk else None) for s in D if has_bloch(D[s], 2)}
    vals = [v for v in per.values() if v is not None]
    return {"seeds": per, "fails": sum(1 for v in vals if v is False), "holds": all(vals) if len(vals) >= 2 else None}


def per_cycle_max(n90, c0=1):
    """{cycle: largest number of large-angle pairs of a sample of the cycle}, cycles >= c0."""
    out = {}
    for x in n90:
        if x[0] >= c0:
            out[int(x[0])] = max(out.get(int(x[0]), 0), int(x[2]))
    return out


def reading_3c(C, D, Dk):
    """3c: a seed takes part if the largest per-cycle number of large-angle pairs of D (cycles >= 1) is above that of C.
    In each such seed the per-cycle number of Dk is at most that of C plus 1 in each cycle >= 1 (same cycle index; a
    cycle of Dk after the last cycle of C is compared with the last cycle of C). Holds only with at least 2 such seeds
    with a Dk value, else None."""
    per = {}
    for s in D:
        if s not in C:
            continue
        mC, mD = per_cycle_max(C[s]["n90"]), per_cycle_max(D[s]["n90"])
        if not mC or not mD or max(mD.values()) <= max(mC.values()):
            continue
        if s not in Dk:
            per[s] = {"pass": None}
            continue
        mK = per_cycle_max(Dk[s]["n90"])
        last = max(mC)
        ok = all(v <= mC.get(c, mC[last]) + 1 for c, v in mK.items()) if mK else None
        per[s] = {"pass": ok, "C": mC, "D": mD, "Dk": mK}
    vals = [v["pass"] for v in per.values() if v["pass"] is not None]
    return {"seeds": per, "holds": all(vals) if len(vals) >= 2 else None}


def row_type(x):
    """3g type of a Bloch row: "thin" (a cell of the largest-angle edge has f below 1/3), "full" (both cells f = 1),
    else "cut"."""
    if min(x[8], x[9]) < THIN:
        return "thin"
    if x[8] == 1.0 and x[9] == 1.0:
        return "full"
    return "cut"


def n_events(rows, n_per_cycle):
    """Rows of the same cube at successive samples count as one event."""
    n, last = 0, {}
    for i, j, k, g in sorted((int(x[2]), int(x[3]), int(x[4]), int(x[0]) * n_per_cycle + int(x[1])) for x in rows):
        if last.get((i, j, k)) != g - 1:
            n += 1
        last[(i, j, k)] = g
    return n


def reading_3g(arms):
    """3g: the type counts of the Bloch rows in the cycles >= 1 per arm and seed, and pooled (with the number of events);
    the place of D: needs at least 2 seeds with at least 5 rows, else "not determined"; "place 1" if at least 80 % of
    the rows are "thin", pooled and in each seed with at least 5 rows; "place 2" the same with "full"; else "mixed"."""
    out = {}
    for arm, runs in arms.items():
        per = {}
        for s, r in runs.items():
            rows = [x for x in r["rows"] if x[0] >= 1]
            npc = int(max((x[1] for x in r["n90"]), default=0)) + 1
            cnt = {t: sum(1 for x in rows if row_type(x) == t) for t in ("thin", "full", "cut")}
            per[s] = dict(cnt, rows=len(rows), events=n_events(rows, npc))
        pooled = {t: sum(v[t] for v in per.values()) for t in ("thin", "full", "cut", "rows", "events")}
        out[arm] = {"seeds": per, "pooled": pooled}
    D = out.get("D", {"seeds": {}, "pooled": {"rows": 0}})
    big = [s for s, v in D["seeds"].items() if v["rows"] >= 5]
    place = "not determined"
    if len(big) >= 2:
        for t, name in (("thin", "place 1"), ("full", "place 2")):
            if (D["pooled"][t] >= 0.8 * D["pooled"]["rows"]
                    and all(D["seeds"][s][t] >= 0.8 * D["seeds"][s]["rows"] for s in big)):
                place = name
                break
        else:
            place = "mixed"
    out["place_D"] = place
    out["seeds_ge5_D"] = big
    return out


def w_tail(stage):
    """w of a run as in 4.16: mean w_loop over the calm tail (the last calm_tail cycles) of a complete stage; None
    without a calm tail. The steady flag is no gate (4.16 gave a value for all 36 jobs, also with steady False)."""
    c = stage["cycles"]
    t = int(stage.get("calm_tail") or 0)
    return float(np.mean([x["w_loop"] for x in c[-t:]])) if t and stage.get("complete") else None


def log_factor(pairs):
    """Mean of ln(factor) +- SE with the 95 % t-interval over the seeds; pairs: {seed: (w_num, w_den)}."""
    from scipy.stats import t as tdist
    x = np.array([math.log(a / b) for a, b in pairs.values()])
    n = len(x)
    if n < 2:
        return {"n": n, "ln": [float(v) for v in x]}
    m, se = float(x.mean()), float(x.std(ddof=1) / math.sqrt(n))
    h = float(tdist.ppf(0.975, n - 1)) * se
    return {"n": n, "per_seed": {s: a / b for s, (a, b) in pairs.items()}, "mean_ln": m, "se": se,
            "ci95_ln": [m - h, m + h], "factor": math.exp(m), "ci95_factor": [math.exp(m - h), math.exp(m + h)]}


def verdict3(a, b, c, few):
    """Verdict of step 3: refuted if 3a is below 50 % or 3b fails in at least 2 seeds; confirmed if 3a, 3b and 3c hold
    and no arm has fewer than 3 values; else not determined."""
    if (a["pooled_share_cut"] is not None and a["pooled_share_cut"] < 0.5) or b["fails"] >= 2:
        return "refuted"
    if a["holds"] and b["holds"] and c["holds"] and not few:
        return "confirmed"
    return "not determined"


def step3(runs, seeds=S3_SEEDS, slope_a=None, ref416=None, fields=True):
    """Readings 3a ... 3g and the verdict of step 3 (PLAN 4.18 with amendment 1); runs/S3<arm>_s<seed>."""
    data, skipped = {a: {} for a in S3_ARMS}, {a: {} for a in S3_ARMS}
    for a in S3_ARMS:
        for s in seeds:
            r = load_step3_run(pathlib.Path(runs) / ("S3%s_s%d" % (a, s)))
            if r is None:
                skipped[a][s] = "no DONE"
            elif not r["init"].get("relax_converged"):
                skipped[a][s] = "relaxation not converged"
            else:
                data[a][s] = r
    few = [a for a in S3_ARMS if len(data[a]) < 3]
    res = {"skipped": skipped, "arms_with_fewer_than_3_values": few,
           "relax": {a: {s: {"iterations": r["init"].get("relax_iterations"), "E_over_Kd": r["init"].get("E_over_Kd")}
                         for s, r in data[a].items()} for a in S3_ARMS},
           "bloch_samples": {a: {s: {"cycles_ge1": sum(1 for b in r["bloch"] if b[0] >= 1 and b[2] + b[3] > 0),
                                     "cycles_ge2": sum(1 for b in r["bloch"] if b[0] >= 2 and b[2] + b[3] > 0),
                                     "samples_ge1": sum(1 for b in r["bloch"] if b[0] >= 1)}
                                 for s, r in data[a].items()} for a in S3_ARMS}}
    res["3a"] = reading_3a(data["D"])
    res["3b"] = reading_3b(data["D"], data["Dk"])
    res["3c"] = reading_3c(data["C"], data["D"], data["Dk"])
    res["3d"] = {a: {s: float(np.mean([x["mean_dt"] for x in r["stage"]["cycles"] if x["cycle"] >= 1]))
                     for s, r in data[a].items()} for a in S3_ARMS}
    if fields:
        e = {}
        for a in ("D", "Dk"):
            for s in data[a]:
                run = pathlib.Path(runs) / ("S3%s_s%d" % (a, s))
                cfg, m, _ = load_run(run)
                ids, f = geometry(cfg)
                v = {"state": state_readings(m, ids, f)}
                npz = run / "diag_fields.npz"
                if npz.exists():
                    th, rh = rho_theta(npz)
                    v["spearman"] = spearman(th, rh)
                    sel = rh >= 0.05
                    if slope_a and sel.any():
                        v["median_ratio_to_a_rho"] = float(np.median(th[sel] / (slope_a * rh[sel])))
                e.setdefault(a, {})[s] = v
        res["3e"] = {"slope_a": slope_a, "arms": e}
    w = {a: {s: w_tail(r["stage"]) for s, r in data[a].items()} for a in S3_ARMS}
    res["3f"] = {"w": w, "flags": {a: {s: {k: r["stage"].get(k) for k in ("steady", "calm_tail", "complete")}
                                        for s, r in data[a].items()} for a in S3_ARMS}}
    for num, den in (("D", "C"), ("Dk", "C"), ("Dk", "D")):
        pairs = {s: (w[num][s], w[den][s]) for s in seeds if w[num].get(s) and w[den].get(s)}
        res["3f"]["%s/%s" % (num, den)] = log_factor(pairs)
    if ref416:
        j = json.loads(pathlib.Path(ref416).read_text())["jobs"]
        res["3f"]["vs_4.16"] = {"%s_s%d" % (a, s): w[a][s] / j["HF%s_s%d_b100" % (a, s)]["w"]
                                for a in ("C", "D") for s in seeds
                                if w[a].get(s) and "HF%s_s%d_b100" % (a, s) in j}
    res["3g"] = reading_3g(data)
    res["verdict"] = verdict3(res["3a"], res["3b"], res["3c"], few)
    return res


def fmt_step3(res):
    """Markdown summary of step 3."""
    L = ["**PLAN 4.18 step 3** (seeds %s; arms C, D, Dk at 100 mT)" % ", ".join(str(s) for s in S3_SEEDS), ""]
    L += ["| arm | seed | relax iterations | samples with a Bloch point, cycles >= 1 | cycles >= 2 | "
          "mean dt (ps), cycles >= 1 | w |", "|---|---|---|---|---|---|---|"]
    for a in S3_ARMS:
        for s in sorted(res["relax"][a]):
            b = res["bloch_samples"][a][s]
            w = res["3f"]["w"][a].get(s)
            L.append("| %s | %d | %s | %d of %d | %d | %.2f | %s |" % (
                a, s, res["relax"][a][s]["iterations"], b["cycles_ge1"], b["samples_ge1"], b["cycles_ge2"],
                1e12 * res["3d"][a][s], "%.4g" % w if w else "–"))
    a3 = res["3a"]
    L += ["", "3a (D): %d samples with a Bloch point in the cycles >= 1; share with a cut cell in the cube %s; share "
          "with a large-angle pair in the cube %s; holds: %s" % (
              a3["samples"], "%.3f" % a3["pooled_share_cut"] if a3["samples"] else "–",
              "%.3f" % a3["pooled_share_large"] if a3["samples"] else "–", a3["holds"])]
    L.append("3b: seeds with a Bloch point of D in the cycles >= 2: %s; fails %d; holds: %s" % (
        res["3b"]["seeds"], res["3b"]["fails"], res["3b"]["holds"]))
    L.append("3c: seeds that take part: %s; holds: %s" % (
        {s: v["pass"] for s, v in res["3c"]["seeds"].items()}, res["3c"]["holds"]))
    for k in ("D/C", "Dk/C", "Dk/D"):
        f = res["3f"][k]
        if "factor" in f:
            L.append("3f %s: factor %.3f [%.3f, %.3f] (exp of the mean ln, SE of ln %.3f, 95 %% t-interval, n %d); "
                     "per seed %s" % (k, f["factor"], *f["ci95_factor"], f["se"], f["n"],
                                      {s: round(v, 3) for s, v in f["per_seed"].items()}))
    g = res["3g"]
    for a in S3_ARMS:
        L.append("3g %s: pooled %s; per seed %s" % (a, g[a]["pooled"], g[a]["seeds"]))
    L += ["3g place of D: %s (seeds with at least 5 rows: %s)" % (g["place_D"], g["seeds_ge5_D"]), "",
          "Verdict of step 3: %s" % res["verdict"]]
    return "\n".join(L) + "\n"


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
    s3 = sub.add_parser("step3")
    s3.add_argument("runs")
    s3.add_argument("--slope_a", type=float, default=None)
    s3.add_argument("--ref416", default=None)
    s3.add_argument("--no_fields", action="store_true")
    s3.add_argument("--out", default=None)
    a = p.parse_args(argv)
    if a.cmd == "step3":
        res = step3(a.runs, slope_a=a.slope_a, ref416=a.ref416, fields=not a.no_fields)
        txt = json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
        if a.out:
            out = pathlib.Path(a.out)
            out.parent.mkdir(parents=True, exist_ok=True)
            out.with_suffix(".json").write_text(txt)
            out.with_suffix(".md").write_text(fmt_step3(res))
        print(fmt_step3(res))
        return
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
