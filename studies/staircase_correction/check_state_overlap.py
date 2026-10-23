#!/usr/bin/env python
"""
Design B of the dynamic test (PROTOKOLL 3.11, 3.15): did a T arm stay in the state that it got from D15?

For each seed and T arm: the reference is the start state, i.e. the prepared state of D15 (init.pt) block averaged onto
the mesh and geometry of the arm (the same transfer_m as run_loops --init_transfer). Overlap = mean of m · m_ref over
the magnetic cells of the arm, for the relaxed state (init.pt) and for the state at the end of the cycles
(checkpoint.pt). Also: the overlap of T3Pc with T3 on the common voxel cells, and the mean m of each particle.

    python check_state_overlap.py --runs /root/runs > results/state_overlap_B.md
"""
import argparse
import pathlib
import sys

import numpy as np
import torch

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "beta_bpeak"))
from anisotropy_model import box_edge  # noqa: E402
from run_loops import transfer_m, locked_geometry  # noqa: E402

SEEDS = (941, 942, 943, 944)


def load(path):
    return torch.load(str(path), map_location="cpu", weights_only=False)


def geometry(cfg, N):
    """ids of the arm, as run_loops builds them (lex units; only the ratios matter). The diameter comes from box_edge
    (d_lex of the config is the nominal value; box_edge sets d for φ at the quantised box edge)."""
    a_lex, d = box_edge(cfg["d_lex"], cfg["phi"], cfg["a_quant"])
    dx = float(cfg["dx_lex"])
    if int(round(a_lex / dx)) != N:
        raise ValueError("box edge %g l_ex / dx %g does not give the mesh of the state (%d cells)" % (a_lex, dx, N))
    g = cfg.get("geom_dx_lex")
    rep = int(round(g / dx)) if g else 1
    sub = (cfg.get("staircase_q", 2) if cfg.get("staircase", "none") != "none" else 8) \
        if cfg.get("surface") == "fraction" else 0
    return locked_geometry(N, d, dx, rep, sub=sub)["ids"]


def overlap(m, m_ref, mag):
    dot = (m[mag] * m_ref[mag]).sum(-1)
    return float(dot.mean()), float((dot < 0).mean())


def particle_means(m, ids):
    return [m[ids == p].mean(0) for p in range(4)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="/root/runs")
    a = ap.parse_args()
    root = pathlib.Path(a.runs)
    print("# Design B: overlap of the T arms with their start state from D15 (PROTOKOLL 3.15)\n")
    print("Overlap: mean of m · m_ref over the magnetic cells (1: same state). Reversed: share of cells with "
          "m · m_ref < 0. m_ref: the prepared state of D15, block averaged onto the arm (the start of the arm).\n")
    print("| seed | arm | overlap after relaxation | overlap at the end of the cycles | reversed at the end |")
    print("|---|---|---|---|---|")
    ends = {}
    for s in SEEDS:
        src = load(root / ("D15_s%d" % s) / "init.pt")
        mf = src["m"].double().numpy()
        ids_f = geometry(src["config"], mf.shape[0])
        for arm in ("T3", "T3Pc", "T15g"):
            d = root / ("%s_s%d" % (arm, s))
            ini = load(d / "init.pt")
            m0 = ini["m"].double().numpy()
            m1 = load(d / "checkpoint.pt")["m"].double().numpy()
            ids = geometry(ini["config"], m0.shape[0])
            ref, _ = transfer_m(mf, ids_f, ids)
            mag = ids >= 0
            o0, _ = overlap(m0, ref, mag)
            o1, r1 = overlap(m1, ref, mag)
            ends[(arm, s)] = (m1, ids)
            print("| %d | %s | %.4f | %.4f | %.4f |" % (s, arm, o0, o1, r1))
    print("\nT3Pc against T3 at the end of the cycles (common voxel cells), and the mean m of each particle "
          "(x, y, z):\n")
    print("| seed | overlap T3Pc · T3 | reversed | particle means T3 | particle means T3Pc |")
    print("|---|---|---|---|---|")
    for s in SEEDS:
        (m3, ids3), (mc, idsc) = ends[("T3", s)], ends[("T3Pc", s)]
        mag = ids3 >= 0
        o, r = overlap(mc, m3, mag)
        f = lambda pm: "; ".join("(%+.2f, %+.2f, %+.2f)" % tuple(v) for v in pm)  # noqa: E731
        print("| %d | %.4f | %.4f | %s | %s |" % (s, o, r, f(particle_means(m3, ids3)), f(particle_means(mc, idsc))))


if __name__ == "__main__":
    main()
