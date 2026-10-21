#!/usr/bin/env python
"""
Why some seeds respond to the stray-field changes and others not (PROTOKOLL 3.30): per run and particle,
the structure of the prepared state (init.pt) and the share of the loop energy.

Per particle k (cells of the arm's geometry, as in check_events.py):
- u: |mean of m| over the cells of the particle in the prepared state (1: uniform; small: vortex or domains);
- f60: the fraction of the cells whose m makes more than 60° with the mean direction of the particle;
- rR60: the mean r/R of these cells (r: distance from the sphere centre, FCC site with minimum image);
- w: loop energy of the particle in the last cycle, mu0 ∮ H dM_k (M_k: mean M along the field over the particle,
  samples_00_amp00.csv), in J/m³; share: w_k · n_k / sum(w · n) with n_k the cells of the particle;
- swing: (max − min) of M_k / Ms in the last cycle.

    python check_seed_structure.py --runs ROOT1 ... --arms T3 T3Pc D15 --events results/events_Q.json \
        --json results/seed_structure.json > results/seed_structure.md
"""
import argparse
import json
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "beta_bpeak"))
sys.path.insert(0, str(HERE / "fe"))
import eval_conv as EC  # noqa: E402
from anisotropy_model import box_edge  # noqa: E402
from geometry import FCC_SITES  # noqa: E402
from units import units  # noqa: E402
from check_events import magnetic_ids, load_m  # noqa: E402

MU0 = 4e-7 * np.pi
EXTRA = 1.3                                     # w / w_base above this: a particle with extra loss
_trapz = getattr(np, "trapezoid", None) or np.trapz


def folder(roots, name):
    return next((r / name for r in roots if (r / name / "DONE").exists()), None)


def last_cycle(csv):
    a = np.loadtxt(csv, delimiter=",", comments="#")
    return a[a[:, 0] == a[:, 0].max()]


def run_structure(d, arm, cache):
    cfg = json.loads((d / "config.json").read_text())
    un = units(cfg["Js"], cfg["A"])
    lex, Ms = un["l_ex"], cfg["Js"] / MU0
    a_lex, d_lex = box_edge(cfg["d_lex"], cfg["phi"], cfg["a_quant"])
    N = int(round(a_lex / cfg["dx_lex"]))
    dx, D, box = cfg["dx_lex"] * lex, d_lex * lex, a_lex * lex
    if (arm, N) not in cache:
        cache[(arm, N)] = magnetic_ids(arm, N, D, dx)
    ids = cache[(arm, N)]
    m = load_m(d / "init.pt")
    xyz = (np.stack(np.indices((N,) * 3), -1) + 0.5) * dx
    c = last_cycle(d / "samples_00_amp00.csv")
    H = c[:, 2]
    parts = []
    for k in range(4):
        sel = ids == k
        mk = m[sel]
        mean = mk.mean(0)
        u = float(np.linalg.norm(mean))
        cosang = mk @ (mean / max(u, 1e-12))
        off = cosang < 0.5
        r = xyz[sel] - FCC_SITES[k] * box
        r -= box * np.round(r / box)
        rR = np.linalg.norm(r, axis=-1) / (0.5 * D)
        Mk = c[:, 9 + k]
        parts.append({"cells": int(sel.sum()), "u": u, "f60": float(off.mean()),
                      "rR60": float(rR[off].mean()) if off.any() else None,
                      "w": float(MU0 * _trapz(H, Mk)), "swing": float((Mk.max() - Mk.min()) / Ms)})
    wn = np.array([p["w"] * p["cells"] for p in parts])
    for p, x in zip(parts, wn):
        p["share"] = float(x / wn.sum()) if wn.sum() != 0 else None
    return parts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True, help="run roots, searched in this order")
    ap.add_argument("--seeds", type=int, nargs="+", default=list(range(941, 953)))
    ap.add_argument("--arms", nargs="+", default=["T3", "T3Pc", "D15"])
    ap.add_argument("--base", type=int, nargs="+", default=[943, 945, 948],
                    help="seeds whose loop energy is the same in the 4 particles (the reference of w / w_base)")
    ap.add_argument("--events", nargs="+", default=[], help="JSON files of check_events.py")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    roots = [pathlib.Path(r) for r in a.runs]
    ev = {}
    for f in a.events:
        ev.update(json.loads(pathlib.Path(f).read_text()))
    cache, out = {}, {}
    for arm in a.arms:
        for s in a.seeds:
            name = "%s_s%d" % (arm, s)
            d = folder(roots, name)
            if d is None:
                continue
            e = ev.get(name)
            out[name] = {"arm": arm, "seed": s, "lnP": EC.run_record(d)["lnP"],
                         "event_particles": sorted(int(k) for k in e["changed"]) if e else None,
                         "particles": run_structure(d, arm, cache)}
    print("Reference of each arm: the seeds %s (mean of the loop energy w of their 4 particles: w_base; mean of their "
          "ln P: ln P_base).\n" % " ".join(map(str, a.base)))
    print("| arm | seed | ln P − ln P_base | w / w_base per particle 0 / 1 / 2 / 3 | event in particles | "
          "u per particle (prepared state) |")
    print("|---|---|---|---|---|---|")
    for arm in a.arms:
        rs = [r for r in out.values() if r["arm"] == arm]
        base = [r for r in rs if r["seed"] in a.base]
        if not base:
            continue
        wb = float(np.mean([p["w"] for r in base for p in r["particles"]]))
        lb = float(np.mean([r["lnP"] for r in base]))
        for r in rs:
            r["w_rel"] = [p["w"] / wb for p in r["particles"]]
            r["lnP_excess"] = r["lnP"] - lb
            evp = r["event_particles"]
            print("| %s | %d | %+.3f | %s | %s | %s |" % (
                arm, r["seed"], r["lnP_excess"], " / ".join("%.2f" % x for x in r["w_rel"]),
                "–" if evp is None else (", ".join(map(str, evp)) or "none"),
                " / ".join("%.2f" % p["u"] for p in r["particles"])))
    ev_all = [(r, k) for r in out.values() if "w_rel" in r for k in (r["event_particles"] or [])]
    print("\n## Events and extra loss\n")
    print("- Particles with an event: %d; with w / w_base > %.1f in the same run: %d; the particle with the largest "
          "w / w_base of the run: %d." % (len(ev_all), EXTRA, sum(r["w_rel"][k] > EXTRA for r, k in ev_all),
                                         sum(r["w_rel"][k] == max(r["w_rel"]) for r, k in ev_all)))
    # the post-event loop of a switched particle is large anyway: compare with an event-free run of the same seed
    # (T3 if it has no event, else T15g), whose loss shows where the extra loss is before any change
    n_ref, n_ext, n_top = 0, 0, 0
    for r, k in ev_all:
        ref = next((out.get("%s_s%d" % (x, r["seed"])) for x in ("T3", "T15g")
                    if out.get("%s_s%d" % (x, r["seed"])) and not out["%s_s%d" % (x, r["seed"])]["event_particles"]
                    and "w_rel" in out["%s_s%d" % (x, r["seed"])]), None)
        if ref is None or ref is r:
            continue
        n_ref += 1
        n_ext += ref["w_rel"][k] > EXTRA
        n_top += ref["w_rel"][k] == max(ref["w_rel"])
    print("- The same particles in an event-free run of the same seed (T3, else T15g): %d event particles; with "
          "w / w_base > %.1f there: %d; the particle with the largest w / w_base there: %d." % (
              n_ref, EXTRA, n_ext, n_top))
    for arm in a.arms:
        sens = [r["seed"] for r in out.values() if r["arm"] == arm and "w_rel" in r and max(r["w_rel"]) > EXTRA]
        print("- %s: seeds with a particle with w / w_base > %.1f: %s (%d of %d)." % (
            arm, EXTRA, " ".join(map(str, sens)), len(sens), sum(r["arm"] == arm for r in out.values())))
    print("\n## Pairs of particles without an event in that particle (count only, not evaluated)\n")
    for x, y in (("T3Pc4", "T3Pc"), ("T3Pc4", "X3q4"), ("T3Pc", "T3")):
        n_all, n_ext, seeds = 0, 0, set()
        for s in a.seeds:
            ra, rb = out.get("%s_s%d" % (x, s)), out.get("%s_s%d" % (y, s))
            if not (ra and rb and "w_rel" in rb):
                continue
            for k in range(4):
                if k in (ra["event_particles"] or []) + (rb["event_particles"] or []):
                    continue
                n_all += 1
                if rb["w_rel"][k] > EXTRA:
                    n_ext += 1
                    seeds.add(s)
        print("- %s − %s: %d particle pairs; %d of them with w / w_base > %.1f in %s (seeds %s)." % (
            x, y, n_all, n_ext, EXTRA, y, " ".join(map(str, sorted(seeds)))))
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
