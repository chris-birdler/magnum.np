#!/usr/bin/env python
"""
Paired field images for the A effect (decision Chris 2026-10-12, option 3): the corner K1 (A 20 pJ/m at d = 1 µm, i.e.
d/l_ex 212; L_eff 30 l_ex = K_eff 1.0 kJ/m³; K_σ 1.5 K_eff; φ 0.62) again with d/l_ex 300 (A 10 pJ/m), the same seeds
2011 … 2013 and all other arguments equal (snapshots at 9 / 50 / 150 mT, accumulators). In this model A enters only
as d/l_ex, so the pair differs only in A. The Herzer-cell grid scales with the box, so the random axes are not the
same cells: the pair shares seed, stress axes, packing and history, not the cell-by-cell microstructure.

    python pair_jobs.py > jobs/pair_K1A10.txt
    python pair_jobs.py --freq > jobs/freq_K1pair.txt      # option C′: f/2 and 2 f
"""
import shlex

from rerun_jobs import job_lines, set_arg

F_REL_30MHZ = 0.0007136497332286743            # f / f_M at 30 MHz, Js = 1.5 T (run_loops default)
FREQS = (("fhalf", 0.5), ("f1", 1.0), ("f2", 2.0))     # f1: only K1 (K1A10 at f = the pair runs), 7.73


def pair_lines():
    """K1A10_k: K1_k with d/l_ex 300; steady rule (--min_calm 2, at most 12 cycles per stage) as the later runs."""
    src = dict(job_lines("jobs/sobol_corners.txt"))
    out = []
    for k in range(3):
        tok = set_arg(shlex.split(src["K1_%d" % k]), "d_lex", ["300.0"])
        tok = set_arg(set_arg(tok, "max_cycles_per_amp", ["12"]), "min_calm", ["2"])
        tok = set_arg(tok, "acc_harm", ["7"])             # harmonics 2f … 7f per cell (option O2, PROTOKOLL 7.72)
        tok = set_arg(tok, "samples", ["512"])            # the C′ sample interval T / 512 at f (decision Chris, 7.76)
        out.append(("K1A10_%d" % k, tok))
    return out


def freq_lines():
    """Option C′ (Chris 2026-10-12): the pair K1 (A 20 pJ/m) and K1A10 (A 10 pJ/m) at f/2 and 2 f (15 and 60 MHz), the
    same seeds, no snapshots, accumulators with the harmonics k ≤ 7 (O2), steady rule. Longest jobs first."""
    from rerun_jobs import strip_snap
    src = dict(job_lines("jobs/sobol_corners.txt"))
    base = [("K1A10_%d" % k, tok) for k, tok in [(k, t) for k, (_, t) in enumerate(pair_lines())]]
    base += [("K1_%d" % k, set_arg(set_arg(shlex.split(src["K1_%d" % k]), "max_cycles_per_amp", ["12"]), "min_calm", ["2"]))
             for k in range(3)]
    out = []
    for tag, fac in FREQS:
        for name, tok in base:
            if fac == 1.0 and name.startswith("K1A10"):
                continue                                    # K1A10 at f: the pair runs (pair_lines, with O2)
            t = strip_snap(list(tok))
            if fac != 1.0:                                  # f1 keeps the default f of the K1 runs
                t = set_arg(t, "f_rel", ["%.10g" % (F_REL_30MHZ * fac)])
            t = set_arg(t, "acc_harm", ["7"])               # harmonics 2f … 7f per cell (option O2, PROTOKOLL 7.72)
            t = set_arg(t, "samples", [str(int(round(512 / fac)))])   # the same sample interval T/512 at f (7.75)
            out.append((fac, "%s_%s" % (name, tag), t))
    out.sort(key=lambda x: (x[0] != 0.5, x[0], not x[1].startswith("K1A10")))   # f/2 of the large box first
    return [(n, t) for _, n, t in out]


if __name__ == "__main__":
    import sys
    lines = freq_lines() if "--freq" in sys.argv else pair_lines()
    for name, tok in lines:
        print(name, " ".join(tok))
