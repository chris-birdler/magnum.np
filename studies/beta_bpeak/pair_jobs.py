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


SWEEP_MHZ = (21.2132, 42.4264, 84.8528, 120.0, 240.0)     # 15·√2, 30·√2, 60·√2, 120, 240 MHz (15 / 30 / 60: C′)


def sweep_lines():
    """Frequency sweep at 9 mT (decision Chris 2026-10-14, option B; PROTOKOLL 7.89): the C′ line of each run (K1A10_k:
    the pair line without snapshots, K1_k: the K1 line with the steady rule), only the 9 mT stage (n_amp 1, b_list
    0.006), started from the relaxed state of the same seed (--init_from <run>/init.pt: the start state of C′),
    harmonics k ≤ 7, samples max(256, 512 · 30 MHz / f). Control: one run per A at 30 MHz (must equal the 9 mT stage
    of C′). Longest jobs first."""
    from rerun_jobs import strip_snap
    src = dict(job_lines("jobs/sobol_corners.txt"))
    base = [("K1A10_%d" % k, tok) for k, (_, tok) in enumerate(pair_lines())]
    base += [("K1_%d" % k, set_arg(set_arg(shlex.split(src["K1_%d" % k]), "max_cycles_per_amp", ["12"]), "min_calm", ["2"]))
             for k in range(3)]
    out = []
    for f in SWEEP_MHZ + (30.0,):
        for name, tok in base:
            if f == 30.0 and not name.endswith("_0"):
                continue                                    # one control run per A
            fac = f / 30.0
            t = strip_snap(list(tok))
            t = set_arg(t, "n_amp", ["1"])
            t = set_arg(t, "b_list", ["0.006"])
            if f != 30.0:                                   # the control keeps the default f, as the C′ runs at f
                t = set_arg(t, "f_rel", ["%.10g" % (F_REL_30MHZ * fac)])
            t = set_arg(t, "acc_harm", ["7"])
            t = set_arg(t, "samples", [str(max(256, int(round(512 / fac))))])
            t = set_arg(t, "init_from", ["%s/init.pt" % name])
            cost = (2.8 if name.startswith("K1A10") else 1.0) / fac
            out.append((cost, "%s_sw%03d" % (name, round(f)), t))
    out.sort(key=lambda x: -x[0])
    return [(n, t) for _, n, t in out]


if __name__ == "__main__":
    import sys
    lines = freq_lines() if "--freq" in sys.argv else (sweep_lines() if "--sweep" in sys.argv else pair_lines())
    for name, tok in lines:
        print(name, " ".join(tok))
