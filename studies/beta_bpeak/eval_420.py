#!/usr/bin/env python
"""PLAN 4.20 evaluation: tied cut cells (E6, --tie_cells), stages A, B and C.

Terms (PLAN 4.20): a cut cell is a magnetic cell with 0 < f < 1; a thin cell is a cut cell with f below l_ex/dx; a tied
cell has no degree of freedom (its m comes from its masters); a free cell is a particle cell that is not tied. The tilt
of a cell is its angle to the sum of m over the full cells of its particle within 2 cells (eval_418.tilt_deg). A
large-angle pair is a pair of neighbour cells of one particle with an angle above 90 deg.

    eval_420.py stageA <runs dir A> --b8 <runs dir of 4.18 step 1> [--out results/plan420_stageA]
        readings of stage A and gate A (arms AV1, AV2 in <runs dir A>/AV<v>_s<seed>; B8 and C8t of 4.18)
    eval_420.py stageB <runs dir> --gateA V1,V2 [--out results/plan420_stageB]
        readings B-0, B-a, B-b, B-c and the selection for stage C (arms C8, D8, D8V1, D8V2, F8 in <runs dir>/<arm>_s<seed>);
        writes <out>.json, <out>.md and <runs dir>/SELECTED (V1, V2 or NONE)
    eval_420.py stageC <runs dir> --variant V1|V2 --ref416 results/hf416.json [--out results/plan420_stageC]
        readings C-a, C-b and the verdict (Dv jobs D<v>_s<seed>_b<mT>, q-16 jobs D16_s<seed>, D<v>16_s<seed>, spare
        seed 925; R0_s921 and HFP_s922 checks)
    eval_420.py r0 <runs dir> --ref416 results/hf416.json
        the control R0: exit code 0 if w of R0_s921 equals w of HFD_s921_b100 of 4.16 within 1e-6 relative, else 1

Bloch flags of a job (rule of eval_rough.py, all cycles from cycle 0, as in 7.142): (its relaxed state has a Bloch
point, a time sample of the cycles has a Bloch point). Bloch samples: the samples of all cycles with a Bloch point.
A job "ran its cycles" if its count file has at least one sample. P = drive frequency x K_d x w (W/cm^3), w as in 4.16
(eval_418.w_tail).
"""
import argparse
import json
import math
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from eval_418 import geometry, load_run, tilt_deg, overlap, read_csv_rows, w_tail, log_factor  # noqa: E402

SEEDS_A = (941, 942, 943, 944)
VARIANTS = {"V1": "thin_face", "V2": "ghost_mirror"}
GATE_A_MIN_SEEDS = 3


def tie_of(cfg, ids, frac):
    """The TieMap of a run (field_terms_extra; V1 needs the sub-cells of the run), or None without --tie_cells."""
    from field_terms_extra import tie_thin_face, tie_ghost_mirror
    from geometry import fcc_subcells
    kind = cfg.get("tie_cells")
    if kind is None:
        return None
    if kind == "thin_face":
        q = cfg["staircase_q"]
        _, _, occ = fcc_subcells(cfg["N_used"], cfg["d"], cfg["dx"], q)
        return tie_thin_face(ids, frac, cfg["dx_lex"], occ, q)
    return tie_ghost_mirror(ids, frac, cfg["N_used"], cfg["d"], cfg["dx"])


def pair_types(m, ids, frac, dx_lex, tied=None):
    """Large-angle pairs (above 90 deg) by type: full-full, cut-full, cut-cut; each also with a thin cell; and the pairs
    with a tied cell. Returns {type: number} and the largest pair angle (deg)."""
    cut = (ids >= 0) & (frac < 1.0)
    thin = cut & (frac < 1.0 / dx_lex)
    t = np.zeros(ids.shape, bool) if tied is None else tied
    out = {"full_full": 0, "cut_full": 0, "cut_cut": 0, "with_thin": 0, "with_tied": 0, "all": 0}
    amax = 0.0
    for ax in range(3):
        same = (ids >= 0) & (ids == np.roll(ids, 1, axis=ax))
        c = (m * np.roll(m, 1, axis=ax)).sum(-1)
        ang = np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))
        if same.any():
            amax = max(amax, float(ang[same].max()))
        big = same & (ang > 90.0)
        nc = cut.astype(int) + np.roll(cut, 1, axis=ax).astype(int)
        out["all"] += int(big.sum())
        out["full_full"] += int((big & (nc == 0)).sum())
        out["cut_full"] += int((big & (nc == 1)).sum())
        out["cut_cut"] += int((big & (nc == 2)).sum())
        out["with_thin"] += int((big & (thin | np.roll(thin, 1, axis=ax))).sum())
        out["with_tied"] += int((big & (t | np.roll(t, 1, axis=ax))).sum())
    return out, amax


def readings_A(run, m_b8):
    """Stage-A readings of one relaxed state (init.pt and summary.json of the run)."""
    cfg, m, init = load_run(run)
    ids, frac = geometry(cfg)
    tie = tie_of(cfg, ids, frac)
    tied = np.zeros(ids.size, bool)
    if tie is not None:
        tied[tie.tied] = True
    tied = tied.reshape(ids.shape)
    th, ref = tilt_deg(m, ids, frac)
    free_cut = (ids >= 0) & (frac < 1.0) & ~tied & ref
    pt, amax = pair_types(m, ids, frac, cfg["dx_lex"], tied)
    s = json.loads((pathlib.Path(run) / "summary.json").read_text())["init"]
    nb = s.get("bloch_virgin", [None, None])
    return {"relax_converged": bool(s.get("relax_converged")), "relax_iterations": s.get("relax_iterations"),
            "n_tied": int(tied.sum()), "free_cut_tilt_max": float(th[free_cut].max()) if free_cut.any() else None,
            "free_cut_gt45": int((th[free_cut] > 45.0).sum()), "pairs_gt90": pt, "pair_angle_max": amax,
            "bloch": nb, "overlap_B8": overlap(m, m_b8, ids) if m_b8 is not None else None,
            "tie": cfg.get("diag", {}).get("tie")}


def gate_A(R):
    """Gate A per variant: passes if, in at least 3 of the 4 seeds, the relaxed state has a converged relaxation, no
    large-angle pair and no Bloch point. A seed without a run counts as failed."""
    out = {}
    for v in VARIANTS:
        ok = []
        for s in SEEDS_A:
            r = R.get(v, {}).get(s)
            ok.append(bool(r and r["relax_converged"] and r["pairs_gt90"]["all"] == 0
                           and r["bloch"][0] == 0 and r["bloch"][1] == 0))
        out[v] = {"seeds_ok": [s for s, o in zip(SEEDS_A, ok) if o], "passes": sum(ok) >= GATE_A_MIN_SEEDS}
    return out


def stageA(runs, b8_runs):
    R = {}
    for arm, d in (("B8", b8_runs), ("C8t", b8_runs)):
        R[arm] = {}
        for s in SEEDS_A:
            run = pathlib.Path(d) / ("%s_s%d" % (arm, s))
            if (run / "init.pt").exists():
                R[arm][s] = readings_A(run, None)
    m_b8 = {s: load_run(pathlib.Path(b8_runs) / ("B8_s%d" % s))[1] for s in SEEDS_A
            if (pathlib.Path(b8_runs) / ("B8_s%d" % s) / "init.pt").exists()}
    for v in VARIANTS:
        R[v] = {}
        for s in SEEDS_A:
            run = pathlib.Path(runs) / ("A%s_s%d" % (v, s))
            if (run / "DONE").exists() and (run / "init.pt").exists():
                R[v][s] = readings_A(run, m_b8.get(s))
    return {"readings": {a: {str(s): r for s, r in x.items()} for a, x in R.items()}, "gate_A": gate_A(R)}


def fmt_stageA(res):
    L = ["# PLAN 4.20 stage A: relaxed states of the tied variants (staircase geometry 48^3, q 8)\n",
         "Pairs above 90 deg by type (ff full-full, cf cut-full, cc cut-cut; thin: with a thin cell; tied: with a tied "
         "cell). Bloch points (+, -) of the relaxed state. Overlap: mean m . m_B8 per particle.\n",
         "| arm | seed | conv. | iterations | tied | free cut tilt max (deg) | free cut > 45 deg | pairs > 90 deg "
         "(all / ff / cf / cc / thin / tied) | largest pair angle (deg) | Bloch | overlap with B8 |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    for a, x in res["readings"].items():
        for s, r in x.items():
            p = r["pairs_gt90"]
            L.append("| %s | %s | %s | %s | %d | %s | %d | %d / %d / %d / %d / %d / %d | %.1f | %s | %s |" % (
                a, s, r["relax_converged"], r["relax_iterations"], r["n_tied"],
                "%.1f" % r["free_cut_tilt_max"] if r["free_cut_tilt_max"] is not None else "-", r["free_cut_gt45"],
                p["all"], p["full_full"], p["cut_full"], p["cut_cut"], p["with_thin"], p["with_tied"],
                r["pair_angle_max"], r["bloch"],
                " ".join("%.4f" % o for o in r["overlap_B8"]) if r["overlap_B8"] else "-"))
    L.append("\n## Gate A (a variant passes with no pair above 90 deg and no Bloch point in at least 3 of 4 seeds)\n")
    for v, g in res["gate_A"].items():
        L.append("- %s (%s): seeds without a pair above 90 deg and a Bloch point: %s; %s" % (
            v, VARIANTS[v], g["seeds_ok"] or "none", "passes" if g["passes"] else "fails"))
    return "\n".join(L) + "\n"


# --- jobs with cycles (stages B and C) -------------------------------------------------------------------------------

def load_job(run):
    """Data of one job with cycles (no DONE needed: the Bloch samples of a failed job count if it ran its cycles), or
    None without summary.json: Bloch flags and samples, w, P, the time step, the run time, exit state."""
    run = pathlib.Path(run)
    if not (run / "summary.json").exists():
        return None
    s = json.loads((run / "summary.json").read_text())
    cfg = json.loads((run / "config.json").read_text()) if (run / "config.json").exists() else {}
    init = s.get("init", {})
    b = read_csv_rows(run / "bloch_00_amp00.csv")
    st = s["stages"][0] if s.get("stages") else None
    w = w_tail(st) if st else None
    cyc = st["cycles"] if st else []
    P = None
    if w is not None and cfg:
        P = cfg["freq"] * cfg["units"]["Kd"] * w / 1e6
    n_tie = (cfg.get("diag") or {}).get("tie", {}).get("n_tied")
    return {"done": (run / "DONE").exists(), "relax_converged": bool(init.get("relax_converged")),
            "relax_iterations": init.get("relax_iterations"),
            "bloch_relaxed": bool(sum(init.get("bloch_virgin") or [0]) > 0),
            "ran_cycles": len(b) > 0, "bloch_samples": int(sum(1 for x in b if x[2] + x[3] > 0)),
            "bloch_cycles": sorted({int(x[0]) for x in b if x[2] + x[3] > 0}),
            "n_cycles": len(cyc), "complete": bool(st and st.get("complete")), "w": w, "P_W_cm3": P,
            "mean_dt": float(np.mean([c["mean_dt"] for c in cyc])) if cyc else None,
            "run_h": (init.get("wall_s", 0.0) + sum(c.get("wall_s", 0.0) for c in cyc)) / 3600.0,
            "n_tied": n_tie}


def has_bloch_job(j):
    """A job fails the Bloch rule if its relaxed state or a time sample of its cycles has a Bloch point."""
    return bool(j["bloch_relaxed"] or j["bloch_samples"] > 0)


# --- stage B ----------------------------------------------------------------------------------------------------------

SEEDS_B = (941, 942, 943, 944)
ARMS_B = ("C8", "D8", "D8V1", "D8V2", "F8")


def fine_reference(m_f, ids_f, frac_f, ids_c):
    """F8 onto the coarse mesh: per coarse cell the f-weighted sum of m of its fine cells (k^3 per coarse cell) of the
    same particle, normalised; cells without such a fine cell are NaN (left out)."""
    Nc, Nf = ids_c.shape[0], ids_f.shape[0]
    k = Nf // Nc
    m = (m_f * frac_f[..., None]).reshape(Nc, k, Nc, k, Nc, k, 3)
    same = ids_f.reshape(Nc, k, Nc, k, Nc, k) == ids_c[:, None, :, None, :, None]
    same &= ids_c[:, None, :, None, :, None] >= 0
    v = (m * same[..., None]).sum(axis=(1, 3, 5))
    n = np.linalg.norm(v, axis=-1)
    out = np.full(v.shape, np.nan)
    ok = n > 0
    out[ok] = v[ok] / n[ok][:, None]
    return out


def overlap_B(m, m_ref, ids, frac):
    """B-b of one arm and seed: per particle the mean of m . m_ref over the coarse cut cells; the cut cells and the full
    cells next to the cut layer (a face neighbour of the same particle is cut) with an angle above 45 deg."""
    cut = (ids >= 0) & (frac < 1.0)
    full = (ids >= 0) & (frac >= 1.0)
    near = np.zeros(ids.shape, bool)
    for ax in range(3):
        for sh in (1, -1):
            near |= np.roll(cut, sh, axis=ax) & (np.roll(ids, sh, axis=ax) == ids)
    near &= full
    ok = ~np.isnan(m_ref[..., 0])
    d = np.where(ok, (m * np.nan_to_num(m_ref)).sum(-1), np.nan)
    per = [float(np.nanmean(np.where(cut & (ids == p), d, np.nan))) for p in range(int(ids.max()) + 1)]
    return {"cut_overlap": per, "cut_overlap_mean": float(np.mean(per)),
            "cut_gt45": int((cut & ok & (d < math.cos(math.radians(45.0)))).sum()),
            "full_next_gt45": int((near & ok & (d < math.cos(math.radians(45.0)))).sum()),
            "cut_left_out": int((cut & ~ok).sum())}


def b_0(J):
    """B-0: F8 is a valid reference if it has values (w) in at least 3 seeds and no Bloch sample in any seed (and no
    Bloch point in a relaxed state)."""
    F = {s: J["F8"].get(s) for s in SEEDS_B}
    n_val = sum(1 for j in F.values() if j and j["w"] is not None)
    n_bloch = sum(1 for j in F.values() if j and has_bloch_job(j))
    return {"seeds_with_value": n_val, "seeds_with_bloch": n_bloch, "valid": bool(n_val >= 3 and n_bloch == 0)}


def b_a(J, variants):
    """B-a: a variant passes if no job of it that ran its cycles has a Bloch sample or a Bloch point in the relaxed
    state; valid only with at least 2 seeds in which D8 has Bloch samples."""
    d8 = [s for s in SEEDS_B if J["D8"].get(s) and J["D8"][s]["bloch_samples"] > 0]
    out = {"D8_seeds_with_bloch": d8, "valid": len(d8) >= 2, "variants": {}}
    for v in variants:
        jobs = {s: J["D8" + v].get(s) for s in SEEDS_B}
        ran = [s for s, j in jobs.items() if j and j["ran_cycles"]]
        bad = [s for s in ran if has_bloch_job(jobs[s])]
        out["variants"][v] = {"seeds_ran": ran, "seeds_with_bloch": bad, "passes": bool(ran and not bad)}
    return out


MIN_SEEDS_RAN = 2


def select(ba, bb, gate_variants):
    """The variant for stage C: of the variants that pass B-a (all of gate A if B-a is not valid), the one with the
    larger mean B-b cut-cell overlap; V1 if the absolute mean paired difference is smaller than its SE or B-b is not
    valid; NONE if no variant passes B-a. A variant whose jobs ran their cycles in fewer than 2 seeds is no candidate
    (stage B incomplete; the plan does not cover it, so it gets no stage C)."""
    ran = [v for v in gate_variants if len(ba["variants"][v]["seeds_ran"]) >= MIN_SEEDS_RAN]
    if not ran:
        return "NONE", "stage B incomplete: no variant ran its cycles in %d seeds" % MIN_SEEDS_RAN
    if ba["valid"]:
        cand = [v for v in ran if ba["variants"][v]["passes"]]
    else:
        cand = list(ran)
    if not cand:
        return "NONE", "no variant passes B-a"
    if len(cand) == 1:
        return cand[0], "only %s passes" % cand[0]
    if bb is None:
        return "V1", "B-b not valid: V1"
    seeds = [s for s in SEEDS_B if s in bb.get("D8V1", {}) and s in bb.get("D8V2", {})]
    if len(seeds) < 2:
        return "V1", "fewer than 2 paired seeds in B-b: V1"
    d = np.array([bb["D8V2"][s]["cut_overlap_mean"] - bb["D8V1"][s]["cut_overlap_mean"] for s in seeds])
    m, se = float(d.mean()), float(d.std(ddof=1) / math.sqrt(len(d)))
    if abs(m) < se:
        return "V1", "mean paired difference V2 - V1 %.4g smaller than its SE %.4g: V1" % (m, se)
    return ("V2" if m > 0 else "V1"), "mean paired difference V2 - V1 %.4g, SE %.4g" % (m, se)


def stageB(runs, gate_variants):
    runs = pathlib.Path(runs)
    J = {a: {s: load_job(runs / ("%s_s%d" % (a, s))) for s in SEEDS_B} for a in ARMS_B}
    res = {"jobs": {a: {str(s): j for s, j in x.items()} for a, x in J.items()}}
    res["B0"] = b_0(J)
    res["Ba"] = b_a(J, gate_variants)
    bb = None
    if res["B0"]["valid"]:
        bb = {}
        for s in SEEDS_B:
            fr = runs / ("F8_s%d" % s)
            if not (fr / "init.pt").exists():
                continue
            cf, mf, _ = load_run(fr)
            ids_f, frac_f = geometry(cf)
            for a in ("C8", "D8", "D8V1", "D8V2"):
                r = runs / ("%s_s%d" % (a, s))
                if not (r / "init.pt").exists():
                    continue
                c, m, _ = load_run(r)
                ids, frac = geometry(c)
                ref = fine_reference(mf, ids_f, frac_f, ids)
                bb.setdefault(a, {})[s] = overlap_B(m, ref, ids, frac)
    res["Bb"] = {a: {str(s): v for s, v in x.items()} for a, x in bb.items()} if bb is not None else None
    Bc = {}
    for a in ARMS_B:
        for ref in ("F8", "D8"):
            if a == ref:
                continue
            pairs = {s: (J[a][s]["w"], J[ref][s]["w"]) for s in SEEDS_B
                     if J[a][s] and J[ref][s] and J[a][s]["w"] and J[ref][s]["w"]}
            if pairs:
                Bc["%s/%s" % (a, ref)] = log_factor(pairs)
    res["Bc"] = Bc
    res["selected"], res["selection_reason"] = select(res["Ba"], bb, gate_variants)
    return res


def fmt_stageB(res):
    L = ["# PLAN 4.20 stage B: staircase geometry, 100 mT, seeds 941 ... 944\n",
         "| arm | seed | done | relaxed conv. | iterations | Bloch relaxed | Bloch samples (cycles) | cycles | complete "
         "| w | P (W/cm^3) | mean dt (s) | run time (h) |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for a, x in res["jobs"].items():
        for s, j in x.items():
            if j is None:
                L.append("| %s | %s | no run |  |  |  |  |  |  |  |  |  |  |" % (a, s))
                continue
            L.append("| %s | %s | %s | %s | %s | %s | %d %s | %d | %s | %s | %s | %s | %.2f |" % (
                a, s, j["done"], j["relax_converged"], j["relax_iterations"], j["bloch_relaxed"], j["bloch_samples"],
                j["bloch_cycles"], j["n_cycles"], j["complete"], "%.5g" % j["w"] if j["w"] else "-",
                "%.4g" % j["P_W_cm3"] if j["P_W_cm3"] else "-", "%.3g" % j["mean_dt"] if j["mean_dt"] else "-",
                j["run_h"]))
    b0, ba = res["B0"], res["Ba"]
    L.append("\n- B-0: F8 has values in %d seeds and Bloch points in %d seeds: %s." % (
        b0["seeds_with_value"], b0["seeds_with_bloch"], "valid reference" if b0["valid"] else "NOT a valid reference"))
    L.append("- B-a: D8 has Bloch samples in seeds %s: B-a %s." % (ba["D8_seeds_with_bloch"],
                                                                  "valid" if ba["valid"] else "NOT valid"))
    for v, x in ba["variants"].items():
        L.append("  - %s: ran its cycles in %s, Bloch in %s: %s" % (v, x["seeds_ran"], x["seeds_with_bloch"] or "none",
                                                                  "passes" if x["passes"] else "fails"))
    if res["Bb"]:
        L.append("\n## B-b: overlap with the f-weighted F8 reference\n")
        L.append("| arm | seed | cut-cell overlap per particle | mean | cut cells > 45 deg | full cells next to the cut "
                 "layer > 45 deg | cut cells left out |")
        L.append("|---|---|---|---|---|---|---|")
        for a, x in res["Bb"].items():
            for s, v in x.items():
                L.append("| %s | %s | %s | %.4f | %d | %d | %d |" % (a, s, " ".join("%.4f" % o for o in v["cut_overlap"]),
                         v["cut_overlap_mean"], v["cut_gt45"], v["full_next_gt45"], v["cut_left_out"]))
    L.append("\n## B-c: factors on P (exp of the mean of ln over the seeds, 95 % t-interval)\n")
    for k, f in res["Bc"].items():
        if "factor" in f:
            L.append("- %s: %.3f [%.3f, %.3f], n %d, SE of ln %.3f" % (k, f["factor"], *f["ci95_factor"], f["n"], f["se"]))
        else:
            L.append("- %s: n %d (too few for an interval)" % (k, f["n"]))
    L.append("\n**Selection for stage C: %s** (%s)" % (res["selected"], res["selection_reason"]))
    return "\n".join(L) + "\n"


# --- stage C ----------------------------------------------------------------------------------------------------------

SEEDS_C = (921, 922, 923, 924)
COUNTED = (921, 923, 924)
AMPS = (9, 50, 100)
# D at 100 mT of the spare seed 925: S3D_s925 of PLAN 4.18 step 3 (summary.json: bloch_virgin [0, 0], 31000 relaxation
# iterations; Bloch samples in 13 cycles >= 1, results/plan418_step3.json)
D_S925_B100 = {"bloch_relaxed": False}


def c_a(D, Dv, D16=None, Dv16=None, spare=None):
    """C-a: counted are the Dv jobs of seeds 921, 923, 924 whose D job (4.16 D or D16) has a relaxed state without a
    Bloch point. A counted job fails if its relaxed state or a time sample of its cycles has a Bloch point (a Bloch
    point in the relaxed state of a seed fails all its jobs). C-a holds if no counted job fails and at least 7 of the
    9 q-8 jobs (2 of the 3 Dv16 jobs) ran their cycles; with fewer, the spare seed 925 (spare: {(seed, mT): job})
    counts in addition; else "not determined". D: {(seed, mT): {"bloch_relaxed": bool}} of 4.16."""
    def count(arm, ref, keys):
        cnt = [k for k in keys if ref.get(k) is not None and not ref[k]["bloch_relaxed"]]
        ran = [k for k in cnt if arm.get(k) and arm[k]["ran_cycles"]]
        relaxed_bad = {k[0] for k in cnt if arm.get(k) and arm[k]["bloch_relaxed"]}
        fail = [k for k in cnt if arm.get(k) and (has_bloch_job(arm[k]) or k[0] in relaxed_bad)]
        return {"counted": [list(k) for k in cnt], "ran": [list(k) for k in ran], "fail": [list(k) for k in fail]}
    out = {"q8": count(Dv, D, [(s, a) for s in COUNTED for a in AMPS])}
    need = 7
    if len(out["q8"]["ran"]) < need and spare:
        sp = count(spare.get("Dv", {}), spare.get("D", {}), [(925, a) for a in AMPS])
        out["spare"] = sp
        out["q8"] = {k: out["q8"][k] + sp[k] for k in ("counted", "ran", "fail")}
    q8_ok = len(out["q8"]["ran"]) >= need
    out["q8"]["status"] = ("fails" if out["q8"]["fail"] else ("holds" if q8_ok else "not determined"))
    if Dv16 is not None:
        r16 = count(Dv16, D16, [(s, 100) for s in COUNTED])
        d16_bloch = [s for s in COUNTED if D16.get((s, 100)) and D16[(s, 100)]["bloch_samples"] > 0]
        r16["D16_seeds_with_bloch"] = d16_bloch
        r16["valid"] = len(d16_bloch) >= 2
        ok16 = len(r16["ran"]) >= 2
        r16["status"] = ("not valid" if not r16["valid"] else
                         "fails" if r16["fail"] else ("holds" if ok16 else "not determined"))
        out["q16"] = r16
    return out


def verdict_C(ca):
    s = [ca["q8"]["status"]] + ([ca["q16"]["status"]] if "q16" in ca else [])
    if "fails" in s:
        return "remedy fails"
    if ca["q8"]["status"] == "holds":
        return "no Bloch points in these runs"
    return "not determined"


def stageC(runs, variant, ref416):
    """Readings of stage C from the box runs and the 4.16 record (results/hf416.json: per job w and the Bloch flags)."""
    runs = pathlib.Path(runs)
    R = json.loads(pathlib.Path(ref416).read_text())["jobs"]
    D = {(s, a): {"bloch_relaxed": bool(R["HFD_s%d_b%d" % (s, a)]["bloch"][0]),
                  "bloch_cycles": R["HFD_s%d_b%d" % (s, a)]["bloch"][1], "w": R["HFD_s%d_b%d" % (s, a)].get("w")}
         for s in SEEDS_C for a in AMPS if "HFD_s%d_b%d" % (s, a) in R}
    Cr = {(s, a): R["HFC_s%d_b%d" % (s, a)].get("w") for s in SEEDS_C for a in AMPS if "HFC_s%d_b%d" % (s, a) in R}
    Dv = {(s, a): load_job(runs / ("D%s_s%d_b%d" % (variant, s, a))) for s in SEEDS_C for a in AMPS}
    D16 = {(s, 100): load_job(runs / ("D16_s%d" % s)) for s in COUNTED}
    Dv16 = {(s, 100): load_job(runs / ("D%s16_s%d" % (variant, s))) for s in COUNTED}
    have16 = any(v is not None for v in Dv16.values())
    sp = {"Dv": {(925, a): load_job(runs / ("D%s_s925_b%d" % (variant, a))) for a in AMPS},
          "D": {(925, a): load_job(runs / ("D_s925_b%d" % a)) for a in (9, 50)}}
    sp["D"] = {k: v for k, v in sp["D"].items() if v is not None}
    sp["D"][(925, 100)] = D_S925_B100
    ca = c_a(D, Dv, D16 if have16 else None, Dv16 if have16 else None, sp)
    res = {"variant": variant, "Ca": ca, "verdict": verdict_C(ca)}
    p = runs / "HFP_s922" / "summary.json"
    if p.exists():
        i = json.loads(p.read_text())["init"]
        res["HFP_s922"] = {"relax_iterations": i.get("relax_iterations"), "E_over_Kd": i.get("E_over_Kd"),
                           "equal_416": bool(i.get("relax_iterations") == 25500
                                             and abs(i.get("E_over_Kd", 0) / 1.6218e-4 - 1.0) < 1e-3)}
    # seed 922 is paired with 4.16 only if its start state was made again as in 4.16
    seeds_b = [s for s in SEEDS_C if s != 922 or res.get("HFP_s922", {}).get("equal_416")]
    # C-b: Bloch samples, ln w and beta against D and C of 4.16
    res["bloch"] = {"Dv": {"%d_%d" % k: (j["bloch_samples"] if j else None) for k, j in Dv.items()}}
    lw = {}
    for ref, W in (("D", {k: v["w"] for k, v in D.items()}), ("C", Cr)):
        for a in AMPS:
            pairs = {s: (Dv[(s, a)]["w"], W[(s, a)]) for s in seeds_b
                     if Dv.get((s, a)) and Dv[(s, a)]["w"] and W.get((s, a))}
            if pairs:
                lw["Dv/%s_%dmT" % (ref, a)] = log_factor(pairs)
    res["ln_w"] = lw
    beta = {}
    for ref, W in (("D", {k: v["w"] for k, v in D.items()}), ("C", Cr)):
        for a0, a1 in ((9, 50), (50, 100), (9, 100)):
            x = []
            for s in seeds_b:
                v = [Dv.get((s, a)) and Dv[(s, a)]["w"] for a in (a0, a1)]
                r = [W.get((s, a)) for a in (a0, a1)]
                if all(v) and all(r):
                    x.append((math.log(v[1] / v[0]) - math.log(r[1] / r[0])) / math.log(a1 / a0))
            if len(x) >= 2:
                from scipy.stats import t as tdist
                m, se = float(np.mean(x)), float(np.std(x, ddof=1) / math.sqrt(len(x)))
                h = float(tdist.ppf(0.975, len(x) - 1)) * se
                beta["Dv-%s_%d_%d" % (ref, a0, a1)] = {"n": len(x), "mean": m, "se": se, "ci95": [m - h, m + h]}
    res["d_beta"] = beta
    if have16:
        res["q16"] = {"D16": {str(k[0]): j for k, j in D16.items()}, "Dv16": {str(k[0]): j for k, j in Dv16.items()}}
    if (runs / "R0_s921" / "summary.json").exists():
        res["R0"] = r0_check(runs, ref416)
    res["jobs"] = {"D%s_s%d_b%d" % (variant, s, a): j for (s, a), j in Dv.items()}
    return res


def r0_check(runs, ref416):
    """The control R0 (PLAN 4.20): w of R0_s921 (the 4.16 arguments of HFD_s921_b100 on the commit of this plan)
    against w of HFD_s921_b100 of 4.16, within 1e-6 relative. A job without w fails."""
    w0 = json.loads(pathlib.Path(ref416).read_text())["jobs"]["HFD_s921_b100"]["w"]
    j = load_job(pathlib.Path(runs) / "R0_s921")
    w = j["w"] if j else None
    rel = abs(w / w0 - 1.0) if w else None
    return {"w": w, "w_416": w0, "rel": rel, "passes": bool(rel is not None and rel <= 1e-6)}


def fmt_stageC(res):
    v = res["variant"]
    L = ["# PLAN 4.20 stage C: variant %s in the 4.16 geometry\n" % v,
         "| job | done | relaxed conv. | Bloch relaxed | Bloch samples (cycles) | complete | w | P (W/cm^3) | run time (h) |",
         "|---|---|---|---|---|---|---|---|---|"]
    for n, j in res["jobs"].items():
        if j is None:
            L.append("| %s | no run |  |  |  |  |  |  |  |" % n)
            continue
        L.append("| %s | %s | %s | %s | %d %s | %s | %s | %s | %.2f |" % (
            n, j["done"], j["relax_converged"], j["bloch_relaxed"], j["bloch_samples"], j["bloch_cycles"],
            j["complete"], "%.5g" % j["w"] if j["w"] else "-", "%.4g" % j["P_W_cm3"] if j["P_W_cm3"] else "-",
            j["run_h"]))
    ca = res["Ca"]
    L.append("\n- C-a (q 8): counted %d jobs, ran their cycles %d, fail %s: %s" % (
        len(ca["q8"]["counted"]), len(ca["q8"]["ran"]), ca["q8"]["fail"] or "none", ca["q8"]["status"]))
    if "q16" in ca:
        L.append("- C-a (q 16): D16 has Bloch samples in seeds %s; counted %d, ran %d, fail %s: %s" % (
            ca["q16"]["D16_seeds_with_bloch"], len(ca["q16"]["counted"]), len(ca["q16"]["ran"]),
            ca["q16"]["fail"] or "none", ca["q16"]["status"]))
    L.append("\n**Verdict: %s**\n" % res["verdict"])
    L.append("## C-b: Delta ln w and Delta beta (paired by seed, mean +- SE, 95 % t-interval)\n")
    for k, f in res["ln_w"].items():
        if "factor" in f:
            L.append("- w %s: factor %.3f [%.3f, %.3f], n %d" % (k, f["factor"], *f["ci95_factor"], f["n"]))
    for k, b in res["d_beta"].items():
        L.append("- beta %s: %+.3f +- %.3f [%+.3f, %+.3f], n %d" % (k, b["mean"], b["se"], *b["ci95"], b["n"]))
    if "R0" in res:
        L.append("\n- R0: w %.8g against %.8g of 4.16, relative %.2g: %s" % (
            res["R0"]["w"] or float("nan"), res["R0"]["w_416"], res["R0"]["rel"] or float("nan"),
            "passes" if res["R0"]["passes"] else "FAILS"))
    if "HFP_s922" in res:
        L.append("- HFP_s922: %s iterations, E/K_d %.5g: %s" % (
            res["HFP_s922"]["relax_iterations"], res["HFP_s922"]["E_over_Kd"],
            "equal to 4.16" if res["HFP_s922"]["equal_416"] else "DIFFERENT from 4.16"))
    return "\n".join(L) + "\n"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("stageA")
    a.add_argument("runs")
    a.add_argument("--b8", required=True)
    a.add_argument("--out", default=str(HERE / "results" / "plan420_stageA"))
    b = sub.add_parser("stageB")
    b.add_argument("runs")
    b.add_argument("--gateA", required=True, help="the variants that passed gate A, for example V1,V2")
    b.add_argument("--out", default=str(HERE / "results" / "plan420_stageB"))
    c = sub.add_parser("stageC")
    c.add_argument("runs")
    c.add_argument("--variant", required=True, choices=list(VARIANTS))
    c.add_argument("--ref416", default=str(HERE / "results" / "hf416.json"))
    c.add_argument("--out", default=str(HERE / "results" / "plan420_stageC"))
    r = sub.add_parser("r0")
    r.add_argument("runs")
    r.add_argument("--ref416", default=str(HERE / "results" / "hf416.json"))
    args = p.parse_args(argv)
    if args.cmd == "r0":
        res = r0_check(args.runs, args.ref416)
        print("R0: w %s against %.8g of 4.16, relative %s: %s" % (res["w"], res["w_416"], res["rel"],
                                                               "passes" if res["passes"] else "FAILS"))
        sys.exit(0 if res["passes"] else 1)
    if args.cmd == "stageB":
        res = stageB(args.runs, [v for v in args.gateA.split(",") if v])
        pathlib.Path(args.out + ".json").write_text(json.dumps(res, indent=1))
        txt = fmt_stageB(res)
        pathlib.Path(args.out + ".md").write_text(txt)
        (pathlib.Path(args.runs) / "SELECTED").write_text(res["selected"] + "\n")
        print(txt)
    if args.cmd == "stageC":
        res = stageC(args.runs, args.variant, args.ref416)
        pathlib.Path(args.out + ".json").write_text(json.dumps(res, indent=1))
        txt = fmt_stageC(res)
        pathlib.Path(args.out + ".md").write_text(txt)
        print(txt)
    if args.cmd == "stageA":
        res = stageA(args.runs, args.b8)
        pathlib.Path(args.out + ".json").write_text(json.dumps(res, indent=1))
        txt = fmt_stageA(res)
        pathlib.Path(args.out + ".md").write_text(txt)
        print(txt)


if __name__ == "__main__":
    main()
