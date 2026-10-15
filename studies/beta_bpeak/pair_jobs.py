#!/usr/bin/env python
"""
Paired field images for the A effect (decision Chris 2026-10-12, option 3): the corner K1 (A 20 pJ/m at d = 1 µm, i.e.
d/l_ex 212; L_eff 30 l_ex = K_eff 1.0 kJ/m³; K_σ 1.5 K_eff; φ 0.62) again with d/l_ex 300 (A 10 pJ/m), the same seeds
2011 … 2013 and all other arguments equal (snapshots at 9 / 50 / 150 mT, accumulators). In this model A enters only
as d/l_ex, so the pair differs only in A. The Herzer-cell grid scales with the box, so the random axes are not the
same cells: the pair shares seed, stress axes, packing and history, not the cell-by-cell microstructure.

    python pair_jobs.py > jobs/pair_K1A10.txt
    python pair_jobs.py --freq > jobs/freq_K1pair.txt      # option C′: f/2 and 2 f
    python pair_jobs.py --atol60 > jobs/atol60.txt         # tolerance test at 60 MHz (7.100, 7.101)
    python pair_jobs.py --atol30 > jobs/atol30.txt         # tolerance twins of the main run at 30 MHz (7.101)
    python pair_jobs.py --long60 > jobs/long60.txt         # long runs at 60 MHz after the C′ history (7.103)
"""
import shlex

from rerun_jobs import job_lines, set_arg

F_REL_30MHZ = 0.0007136497332286743            # f / f_M at 30 MHz, Js = 1.5 T (run_loops default)
FREQS = (("fhalf", 0.5), ("f1", 1.0), ("f2", 2.0))     # f1: only K1 (K1A10 at f = the pair runs), 7.73
INIT_ATOL = 1e-5                                # atol of the relaxed states K1_k/init.pt (run_loops default)


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


def atol60_lines():
    """Long tolerance test at 60 MHz (decision Chris 2026-10-14, option 2, PROTOKOLL 7.100): the sweep line of the
    small particle (K1_k, A 20 runs) at 60 MHz, 9 mT only, exactly 40 cycles, with --atol 1e-5 (L5) and 1e-6 (L6).
    Reference: the 9 mT stage of K1_k_f2 (C′); the sweep protocol equals the C′ first stage (control at 30 MHz,
    7.91). atol is an INIT_KEY (7.88): K1_k/init.pt was relaxed with atol 1e-5, so only L5 starts from it; L6 relaxes
    again with atol 1e-6 from the same seed (PROTOKOLL 7.101). Longest (1e-6) first."""
    out = []
    for tol, tag in (("1e-06", "L6"), ("1e-05", "L5")):
        for k in range(3):
            _, tok = next((n, t) for n, t in sweep_lines() if n == "K1_%d_sw085" % k)
            t = set_arg(list(tok), "f_rel", ["%.10g" % (F_REL_30MHZ * 2.0)])
            t = set_arg(t, "samples", ["256"])
            t = set_arg(t, "atol", [tol])
            t = set_arg(t, "cycles_per_amp", ["40"])
            t = set_arg(t, "max_cycles_per_amp", ["40"])
            if float(tol) != INIT_ATOL:
                t = set_arg(t, "init_from", None)
            out.append(("K1_%d_sw060%s" % (k, tag), t))
    return out


def atol30_lines():
    """Tolerance twins of the main run at 30 MHz (decision Chris 2026-10-14, option 1 with patch, PROTOKOLL 7.101): the
    production line of K1_k (jobs/sobol_corners.txt: 7 amplitudes, 7 cycles each, virgin relaxation) with --atol 1e-6.
    atol is an INIT_KEY (7.88), so each twin relaxes again from its seed with atol 1e-6. Only the snapshots are removed
    (disk, no physics). Reference: the main runs K1_k (atol 1e-5, the same line)."""
    from rerun_jobs import strip_snap
    src = dict(job_lines("jobs/sobol_corners.txt"))
    return [("K1_%d_L6" % k, set_arg(strip_snap(shlex.split(src["K1_%d" % k])), "atol", ["1e-06"])) for k in range(3)]


def long60_lines():
    """Long runs at 60 MHz after the C′ history (decision Chris 2026-10-15, point 5 of the open questions; PROTOKOLL
    7.103): the large particle (K1A10_k_f2) at 9 mT and the small particle (K1_k_f2) at 20 and 35 mT, each with the C′
    history up to the tested amplitude (earlier stages with the C′ stop rule, ≤ 12 cycles, ≥ 2 calm) and then exactly
    40 cycles at it (rerun_jobs.cond). No --init_from: C′ relaxed in the same process, so the history repeats bit for bit
    (7.59). Longest jobs (large particle: relaxation 39 min, 122 s per cycle) first."""
    from rerun_jobs import cond
    big = cond(["K1A10_%d_f2" % k for k in range(3)], "jobs/freq_K1pair.txt", amps=(9,), max_prev=12, last=40)
    small = cond(["K1_%d_f2" % k for k in range(3)], "jobs/freq_K1pair.txt", amps=(35, 20), max_prev=12, last=40)
    return [tuple(x.split(" ", 1)) for x in big + small]


if __name__ == "__main__":
    import sys
    if "--long60" in sys.argv:
        for name, args in long60_lines():
            print(name, args)
        sys.exit(0)
    if "--atol60" in sys.argv or "--atol30" in sys.argv:
        for name, tok in (atol60_lines() if "--atol60" in sys.argv else atol30_lines()):
            print(name, " ".join(tok))
        sys.exit(0)
    lines = freq_lines() if "--freq" in sys.argv else (sweep_lines() if "--sweep" in sys.argv else pair_lines())
    for name, tok in lines:
        print(name, " ".join(tok))
