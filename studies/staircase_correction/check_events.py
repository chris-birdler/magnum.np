#!/usr/bin/env python
"""
Design B of the dynamic test (PROTOKOLL 3.26): local events in the particles during the cycles.

For each run: the cells inside the particles whose m changed between the start of the cycles (init.pt) and the end
(checkpoint.pt), m · m0 < 0.9, per particle with the range of r/R (distance from the sphere centre over the radius; the
centre is the FCC site with the minimum image, not the mean of the cell positions, which is wrong for a particle that
crosses the box edge). Also the loop energy W_loop of each cycle (summary.json), which shows when a change occurs.
Magnetic cells: the voxel geometry (T3, D15), the cells of the q-2 staircase (T3Pc, X3), the locked dx-3 geometry (T15g).

    python check_events.py --runs ROOT1 ROOT2 ... --seeds 941 ... 952 --out results/events_B   (writes .md and .json)
"""
import argparse
import json
import pathlib
import sys

import numpy as np
import torch

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "beta_bpeak"))
sys.path.insert(0, str(HERE / "fe"))
from anisotropy_model import box_edge  # noqa: E402
from geometry import fcc_box, FCC_SITES  # noqa: E402
from units import units  # noqa: E402
from fe_states import fine_ids  # noqa: E402
from magnumnp.field_terms.demag_staircase import staircase_cells  # noqa: E402

ARMS = ("T3", "T3Pc", "X3", "T15g", "D15")


def load_m(path):
    return torch.load(str(path), map_location="cpu", weights_only=False)["m"].double().numpy()


def magnetic_ids(arm, N, d, dx):
    """Particle index of the magnetic cells of the arm, -1: void."""
    if arm in ("T3", "D15"):
        return fcc_box(N, d, dx)["ids"]
    if arm in ("T3Pc", "X3"):
        n2, pid2, _ = staircase_cells(fine_ids(N, 2, d, dx), 2)
        return np.where(n2 > 0, pid2, -1)
    ids = fcc_box(N // 2, d, 2 * dx)["ids"]                     # T15g: the voxel geometry of dx 3 on the dx-1.5 mesh
    return ids.repeat(2, 0).repeat(2, 1).repeat(2, 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True, help="run roots, searched in this order")
    ap.add_argument("--seeds", type=int, nargs="+", default=list(range(941, 953)))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    roots = [pathlib.Path(r) for r in a.runs]
    cache, out, rows = {}, {}, []
    for arm in ARMS:
        for s in a.seeds:
            name = "%s_s%d" % (arm, s)
            d = next((r / name for r in roots if (r / name / "DONE").exists()), None)
            ck = next((r / name / "checkpoint.pt" for r in roots if (r / name / "checkpoint.pt").exists()), None)
            if d is None or ck is None:
                rows.append("| %s | %d | missing | | | |" % (arm, s))
                continue
            cfg = json.loads((d / "config.json").read_text())
            lex = units(cfg["Js"], cfg["A"])["l_ex"]
            a_lex, d_lex = box_edge(cfg["d_lex"], cfg["phi"], cfg["a_quant"])
            N = int(round(a_lex / cfg["dx_lex"]))
            dx, D, box = cfg["dx_lex"] * lex, d_lex * lex, a_lex * lex
            if (arm, N) not in cache:
                cache[(arm, N)] = magnetic_ids(arm, N, D, dx)
            ids = cache[(arm, N)]
            dot = (load_m(d / "init.pt") * load_m(ck)).sum(-1)
            xyz = (np.stack(np.indices((N,) * 3), -1) + 0.5) * dx
            W = [c["W_loop"] for c in json.loads((d / "summary.json").read_text())["stages"][0]["cycles"]]
            parts = {}
            for k in range(4):
                sel = (ids == k) & (dot < 0.9)
                if sel.any():
                    r = xyz[sel] - FCC_SITES[k] * box
                    r -= box * np.round(r / box)
                    rr = np.linalg.norm(r, axis=-1) / (0.5 * D)
                    parts[k] = {"cells": int(sel.sum()), "rR_min": float(rr.min()), "rR_max": float(rr.max()),
                                "dot_min": float(dot[sel].min())}
            out[name] = {"W_loop": W, "changed": parts}
            txt = "; ".join("particle %d: %d cells, r/R %.2f … %.2f, min m·m0 %.2f" % (
                k, p["cells"], p["rR_min"], p["rR_max"], p["dot_min"]) for k, p in parts.items()) or "none"
            rows.append("| %s | %d | %s | %s | %.3f | %s |" % (arm, s, txt, " / ".join("%.3g" % w for w in W),
                                                              W[0] / W[-1], "yes" if parts else "no"))
            print(rows[-1], flush=True)
    md = ["# Design B: local events in the particles during the cycles (PROTOKOLL 3.26)\n",
          "Changed cells: cells inside the particles with m · m0 < 0.9 between the start of the cycles and the end. "
          "W_loop: loop energy of each cycle in J/m³ (cycle 0 starts from rest, so in a run without a change it is "
          "about 0.45 … 0.57 of the last cycle).\n",
          "| arm | seed | changed cells per particle | W_loop per cycle (J/m³) | W_loop cycle 0 / last | event |",
          "|---|---|---|---|---|---|"] + rows
    pathlib.Path(a.out + ".md").write_text("\n".join(md) + "\n")
    pathlib.Path(a.out + ".json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
