#!/usr/bin/env python
"""
History test (PROTOKOLL 7.123, 7.125; Chris 2026-10-19): what moves inside the particles between the two zero
crossings of the drive field. The snapshots of the history runs are at 180° (H = 0, falling) and 360° (H = 0, rising)
of each cycle (drive H0 sin(ωt)). At H = 0 the difference d = m(360°) − m(180°) is the out-of-phase part of the
response, thus the lossy part. Per cycle and run, over the particle cells, with x = drive direction:
  net       mean of d_x (= 2 M_r / M_s of the loop at H = 0; the part that reaches the box average)
  with      mean of max(d_x, 0); against: mean of max(−d_x, 0) (cells that move with / against the net change)
  cancel    against / with (0: all cells move with the net change; 1: the motion cancels)
  perp      mean of |d_⊥| (motion across the drive direction)
  rot       mean of arccos(m(180°) · m(360°)) / π (sign-blind rotation; 1 = full reversal)
  rot/net   rotation per net change (how much rotation one unit of net change costs)
  particles net per particle (do particles move against each other?)
Report per run: the mean over the cycles from --first on (default 3) and the cycles before; ad against nx per
microstructure.

    python eval_hist_snapshots.py --runs <runs folder> --json out.json K4_0_ad_h16 K4_0_nx_h16 ... > report.md
"""
import argparse
import json
import pathlib
import re

import numpy as np

import domains as DM

SNAP = re.compile(r"m_amp00_c(\d+)_ph(180|360)\.vti$")


def cycle_measures(s180, s360, e_x):
    m1, m2, mask, ids = s180["m"], s360["m"], s360["mask"], s360["ids"]
    d = (m2 - m1)[mask]
    dx = d @ e_x
    perp = np.linalg.norm(d - dx[:, None] * e_x[None, :], axis=1)
    rot = np.arccos(np.clip((m1[mask] * m2[mask]).sum(-1), -1.0, 1.0)) / np.pi
    out = {"net": float(dx.mean()), "with": float(np.clip(dx, 0, None).mean()),
           "against": float(np.clip(-dx, 0, None).mean()), "perp": float(perp.mean()), "rot": float(rot.mean())}
    out["cancel"] = out["against"] / out["with"] if out["with"] > 0 else float("nan")
    out["rot_per_net"] = out["rot"] / abs(out["net"]) if out["net"] != 0 else float("nan")
    pid = ids[mask]
    out["particles"] = {int(p): float(dx[pid == p].mean()) for p in np.unique(pid)}
    return out


def run_measures(run_dir):
    cfg = json.loads((run_dir / "config.json").read_text())
    e_x = np.asarray(cfg.get("direction", [1.0, 0.0, 0.0]), float)
    e_x /= np.linalg.norm(e_x)
    files = {}
    for f in run_dir.glob("m_amp00_c*_ph*.vti"):
        mm = SNAP.search(f.name)
        if mm:
            files.setdefault(int(mm.group(1)), {})[mm.group(2)] = f
    res = {}
    for c in sorted(files):
        if set(files[c]) != {"180", "360"}:
            continue
        res[c] = cycle_measures(DM.load_state(files[c]["180"]), DM.load_state(files[c]["360"]), e_x)
    return res


KEYS = ("net", "with", "against", "cancel", "perp", "rot", "rot_per_net")


def mean_over(res, cycles):
    v = [res[c] for c in cycles if c in res]
    return {k: float(np.mean([x[k] for x in v])) for k in KEYS} if v else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True)
    ap.add_argument("--json", default=None)
    ap.add_argument("--first", type=int, default=3)
    ap.add_argument("names", nargs="+")
    a = ap.parse_args()
    root = pathlib.Path(a.runs)
    allres = {n: run_measures(root / n) for n in a.names}
    L = ["# History test: motion between the zero crossings of the drive (eval_hist_snapshots.py)", "",
         "d = m(360°) − m(180°) over the particle cells; x = drive direction. Means over the cycles from %d on "
         "(late) and before (early)." % a.first, "",
         "| run | cycles | net | with | against | cancel | perp | rot | rot/net | net per particle (late) |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for n, res in allres.items():
        cyc = sorted(res)
        for lab, sel in (("early", [c for c in cyc if c < a.first]), ("late", [c for c in cyc if c >= a.first])):
            mo = mean_over(res, sel)
            if mo is None:
                continue
            parts = ""
            if lab == "late":
                pv = {}
                for c in sel:
                    for p, x in res[c]["particles"].items():
                        pv.setdefault(p, []).append(x)
                parts = " / ".join("%d: %+.4f" % (p, np.mean(x)) for p, x in sorted(pv.items()))
            L.append("| %s | %s %d…%d (%d) | %+.4f | %.4f | %.4f | %.2f | %.4f | %.4f | %.1f | %s |" % (
                n, lab, sel[0], sel[-1], len(sel), mo["net"], mo["with"], mo["against"], mo["cancel"], mo["perp"],
                mo["rot"], mo["rot_per_net"], parts))
    print("\n".join(L))
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps({n: {str(c): r for c, r in res.items()} for n, res in allres.items()},
                                                   indent=1))


if __name__ == "__main__":
    main()
