#!/usr/bin/env python
"""
Job lists of the beta(B_peak) study (fictitious nanocrystalline powder,
reduced units) and a cost estimate.

    python jobs.py                  # writes jobs/*.txt, prints the cost table
    python jobs.py --speedup 2.2    # cost for a GPU 2.2x faster than a V100

Job file format: one job per line, "<name> <run_loops.py arguments>", sorted
by estimated cost (largest first).

Groups:
  bench     1 cycle of a d/l_ex = 300 run (timing -> calibrate the cost model)
  meshtest  d/l_ex = 96, dx = 3 / 2 / 1.5 l_ex, floor (Q_eff = 0) and L_eff = 12:
            does the grid pin the unresolved vortex cores? (ran 2026-09-30: cycles
            not steady at b <= 0.25 within 3-4 cycles -> steady test)
  steady    same geometry, 2 amplitudes (b ~ 0.25 and 0.06), 12 cycles each,
            dx = 3 vs 1.5, floor and L_eff = 12, alpha = 0.02 and 0.1 (L12):
            mesh test on steady cycles. Runs BEFORE the pilot.
  numfloor  same geometry, dx = 3, L_eff = 12, virgin state, B_peak = 9 / 50 / 150 mT
            ascending, 7 cycles: numerical floor of the small loops (atol 1e-5 vs
            1e-6, fp32 vs fp64).
  initproto same geometry, L_eff = 12, B_peak = 9 / 50 / 150 mT, 7 cycles: initial
            state and order of the amplitudes. A = virgin state (random m, full
            relaxation, like the anneal without field), ascending, dx = 1.5 (dx = 3:
            F_si_atol1e-05);
            B = AC demagnetization (decaying saturating cycles), ascending;
            C = present protocol (1 saturating cycle, descending).
  conv      same geometry, alpha = 0.1, virgin state (the same start on every mesh:
            random m on 6 l_ex blocks), B_peak = 9 / 20 / 50 mT: mesh convergence
            dx = 2 / 1.5 / 1; at dx = 1.5 also a second initial state (init_seed 2)
            and 50 mT directly from the virgin state (is the ascending history needed?).
  scatter   as conv at dx = 1.5, B_peak = 9 / 50 mT, seeds 3, 4, 5 (other cubes and
            other initial states): realisation scatter (uses the free GPU during C_dx1).
  freqtest  PLAN.md step 1: d/l_ex 300, dx 3, virgin state, alpha 0.1, 7 amplitudes
            9 ... 150 mT; f and 2 f; base point (L_eff 12, r_p 0, seeds 901-904)
            and corner (L_eff 30, r_p 3, seeds 905-908).
  relaxtest virgin state relaxed to convergence (--relax_mode converge): d/l_ex 300,
            dx 3, base point, alpha 0.1 and 0.01, seeds 901-903, 9 and 20 mT; compare
            the 9 mT energy balance and loss with F_base_f_* and A010_* (fixed relax).
  smoothtest surface roughness: dx 3 with a volume-fraction surface (--surface fraction)
            at d/l_ex 212 (seeds 911-914, as R212_*) and 300 (seeds 1-3, one job per
            amplitude, as P1_d300_dx3_*); compare with the staircase dx 3 and dx 1.5.
  alphatest is alpha a real factor at 30 MHz? Base point as freqtest (d/l_ex 300,
            dx 3, L_eff 12, r_p 0, phi 0.65), alpha 0.01 and 0.03, seeds 901-903,
            9 / 20 / 50 / 100 / 150 mT. alpha 0.1: F_base_f_s901-903 (freqtest).
  alphatest2 alpha 0.01, seeds 904-906, otherwise as alphatest, but the virgin state relaxed
            to convergence (hard gate; PROTOKOLL 7.13).
  snaptest  domain structure: base point, alpha 0.1, relaxed to convergence, 20 and 150 mT,
            seeds 901 and 902, m as .vti at 4 phases of the last cycle (--snap_phases 4).
  baseline  repeat after the audit (PROTOKOLL 7.13): realisation scatter at 1 um with a
            converged virgin state: PROD (d/l_ex 300, dx 3, 7 amplitudes), base point,
            alpha 0.1, seeds 921-928; m as .vti at 4 phases of the last cycle of each amplitude.
  mesh300   repeat after the audit: dx 1.5 vs dx 3 at d/l_ex 300 with a converged virgin
            state, seeds 921-924, 9 / 50 / 100 mT, 5 cycles; the 9 mT job relaxes and
            writes init.pt, the 50 and 100 mT jobs start from it (--init_from, --init_wait_h).
  freq9     check before the Sobol design: f/2 and 2 f at 9 and 20 mT, seeds 921-923, from the relaxed
            virgin state of B300 (f is in the baseline): is the loss at small B dynamic or hysteretic?
  phase16   check before the Sobol design: seeds 921, 922, 9 / 20 / 150 mT, 16 snapshots per cycle in the
            last 3 cycles, and a ring-down test (H held at H max for 2 periods, p_dis at 512 per period).
  rough300  surface roughness with the converged state: as the dx 3 jobs of mesh300 (seeds 921-924,
            9 / 50 / 100 mT, 5 cycles) but --surface fraction; compared with M300_dx3_* (staircase).
  roughtest PLAN.md step 1b: mesh error at the small end of the d/l_ex range:
            d/l_ex 212, dx 3 vs dx 1.5, 4 seeds each (911-914), 9 and 50 mT,
            base point (L_eff 12, r_p 0, phi 0.65). Compared with d/l_ex 300
            (planning pilot): same mesh error -> d/l_ex effect clean.
  planpilot virgin state, alpha = 0.1, B_peak = 9 / 50 mT (PROTOKOLL 7.8):
            block 1 scatter vs size (d/l_ex 150 and 300) and dx 3 vs 1.5 at 300,
            block 2 effect sizes (L_eff 30, r_p 3) at d/l_ex 96,
            block 3 mesh pairs dx 1 / 2 / 3 at d/l_ex 96 (seeds 3-5).
  demagtest same geometry, alpha = 0.1, B_peak = 9 / 50 / 150 mT: AC demagnetization at
            f and at 4 f (only the end state counts), dx = 3 and 1.5, and the
            virgin state at dx = 1.5 as the reference of a clean virgin state.

Amplitudes: the study range is B_peak = 10 ... 150 mT (Chris, 2026-09-30):
b = B_peak/Js with Js = 1.5 T. run_loops.py --b_list sets the drive per
stage from the measured b/h of the stage before.
  pilot   d/l_ex = 150, 4 amplitudes: pinning signal against the floor, mesh,
          damping, drive direction, frequency, realisation scatter. The
          production matrix is decided after the pilot (PROTOKOLL.md).

Cost model (calibrated on V100 2026-10-02 with the planning pilot, alpha = 0.1):
  steps/cycle = 2 pi / (f_rel * dt_tau),  dt_tau = dt * gamma Ms
                alpha = 0.1:  dt_tau = 1.62 (dx/3)^1.5  (measured 1.61 / 0.91 / 0.56 / 0.26
                              at dx = 3 / 2 / 1.5 / 1)
                alpha = 0.02: dt_tau = 0.75 (dx/3)^0.5  (older runs, 2026-09-30)
  s/step      = max(s_min_step, s_per_step_Mcell * n_cells / 1e6) / speedup
                (0.012 s/Mcell up to 10 M cells, 0.0164 above (26 M cells measured);
                small grids are limited by the overhead s_min_step = 0.008 s)
  relaxation of the virgin state: about 2 cycles
  cycles      = 1 + n_stages * cycles_per_amp   (fixed number of cycles), else
                1 + n_stages * cycles_avg
                cycles_avg = 3.5: with 10 % scatter per cycle (target, dW_tol = 0.10)
                about half of the stages need the 4th cycle
The unit "V100-h" is only a cost unit. Any fp32 GPU can run the jobs.
"""
import argparse
import math
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from anisotropy_model import box_edge   # noqa: E402
from units import F_REL_DEFAULT         # noqa: E402

JS_T = 1.5                                         # Js of the unit system (units.py)


def b_of_mT(*mT):
    return [round(x / 1000.0 / JS_T, 6) for x in mT]


B_STUDY = b_of_mT(150, 100, 70, 50, 35, 20, 9)     # production amplitudes (9 mT: 10 mT inside the range)
B_PILOT = b_of_mT(150, 100, 50, 25, 9)

BASE = dict(d_lex=300.0, dx_lex=3.0, phi=0.65, Leff_lex=18.0, r_p=0.0,
            alpha=0.02, b_list=B_STUDY, n_amp=len(B_STUDY), cycles_per_amp=7, max_cycles_per_amp=7,
            samples=256, precision="single", seed=1)
PILOT = dict(BASE, d_lex=150.0, b_list=B_PILOT, n_amp=len(B_PILOT))
# frozen protocol of the mesh test and the steady test (they ran with it)
_OLD = dict(d_lex=150.0, dx_lex=3.0, phi=0.65, Leff_lex=18.0, r_p=0.0,
            alpha=0.02, n_amp=4, h_factor=0.25, cycles_per_amp=3, max_cycles_per_amp=4, dW_tol=0.10,
            samples=256, precision="single", seed=1)
MESHTEST = dict(_OLD, d_lex=96.0)       # a = 144 l_ex: a, a/2 and L_eff on cell faces for dx = 3, 2, 1.5
# steady test: b ~ 0.25 and 0.06 (from the meshtest), 12 cycles each, no early stop
STEADY = dict(MESHTEST, h_list=[0.02917, 0.007292], n_amp=2, cycles_per_amp=12, max_cycles_per_amp=12)


def job(base, name, **kw):
    a = dict(base)
    a.update(kw)
    return name, a


def matrix():
    bench = [job(BASE, "bench_base", max_cycles=1, relax_tau=100.0)]

    meshtest = []
    for dx in (3.0, 2.0, 1.5):
        tag = ("%g" % dx).replace(".", "")
        meshtest.append(job(MESHTEST, "M_floor_dx%s" % tag, Q_eff=0.0, dx_lex=dx))
        meshtest.append(job(MESHTEST, "M_L12_dx%s" % tag, Leff_lex=12.0, dx_lex=dx))

    steady = []
    for dx in (3.0, 1.5):
        tag = ("%g" % dx).replace(".", "")
        steady.append(job(STEADY, "S_floor_dx%s" % tag, Q_eff=0.0, dx_lex=dx))
        steady.append(job(STEADY, "S_L12_dx%s" % tag, Leff_lex=12.0, dx_lex=dx))
        steady.append(job(STEADY, "S_L12_a010_dx%s" % tag, Leff_lex=12.0, dx_lex=dx, alpha=0.1))

    SMALL = dict(MESHTEST, Leff_lex=12.0, b_list=b_of_mT(9, 50, 150), n_amp=3,
                 cycles_per_amp=7, max_cycles_per_amp=7)
    del SMALL["h_factor"], SMALL["dW_tol"]
    numfloor = []
    for prec in ("single", "double"):
        for atol in (1e-5, 1e-6):
            numfloor.append(job(SMALL, "F_%s_atol%.0e" % (prec[:2], atol), protocol="virgin_asc",
                                precision=prec, atol=atol))
    # A at dx = 3 is F_si_atol1e-05 (same arguments)
    initproto = [job(SMALL, "I_A_virgin_dx15", protocol="virgin_asc", dx_lex=1.5),
                 job(SMALL, "I_B_acdemag_dx3", protocol="acdemag_asc"),
                 job(SMALL, "I_C_reset_dx3", protocol="reset_desc")]

    pilot = []
    for L in (12.0, 18.0, 24.0, 30.0):                              # Q_eff = 6.9e-3 ... 1.1e-3
        pilot.append(job(PILOT, "T_L%02.0f" % L, Leff_lex=L))
    pilot.append(job(PILOT, "T_floor", Q_eff=0.0))                   # no anisotropy: floor
    for rp in (1.0, 3.0):                                           # particle-scale stress
        pilot.append(job(PILOT, "T_L18_rp%.0f" % rp, r_p=rp))
    pilot.append(job(PILOT, "T_L18_alpha010", alpha=0.1))           # damping loss
    pilot.append(job(PILOT, "T_L18_d111", direction=[1.0, 1.0, 1.0]))  # cube lattice vs drive
    pilot.append(job(PILOT, "T_L18_f2", f_rel=2.0 * F_REL_DEFAULT))  # dynamic share (pair with T_L18)
    pilot.append(job(PILOT, "T_L30_dx2", Leff_lex=30.0, dx_lex=2.0))  # mesh, same cubes as T_L30
    pilot.append(job(PILOT, "T_L30_seed2", Leff_lex=30.0, seed=2))    # realisation scatter
    D = dict(SMALL, alpha=0.1)
    demagtest = []
    for dx in (3.0, 1.5):
        tag = ("%g" % dx).replace(".", "")
        demagtest.append(job(D, "D_B_dx%s" % tag, protocol="acdemag_asc", dx_lex=dx))
        demagtest.append(job(D, "D_B4f_dx%s" % tag, protocol="acdemag_asc", demag_f_factor=4.0, dx_lex=dx))
    demagtest.append(job(D, "D_A_dx15", protocol="virgin_asc", dx_lex=1.5))
    CV = dict(D, protocol="virgin_asc", init_block=6.0, b_list=b_of_mT(9, 20, 50), n_amp=3)
    conv = []
    for dx in (2.0, 1.5, 1.0):
        conv.append(job(CV, "C_dx%s" % ("%g" % dx).replace(".", ""), dx_lex=dx))
    conv.append(job(CV, "C_dx15_is2", dx_lex=1.5, init_seed=2))
    conv.append(job(CV, "C_dx15_direct50", dx_lex=1.5, b_list=b_of_mT(50), n_amp=1))
    scatter = [job(CV, "R_dx15_s%d" % sd, dx_lex=1.5, seed=sd, b_list=b_of_mT(9, 50), n_amp=2) for sd in (3, 4, 5)]
    PP = dict(CV, dx_lex=1.5, b_list=b_of_mT(9, 50), n_amp=2)
    pp = []
    for sd in (1, 2, 3, 4):                                            # block 1
        pp.append(job(PP, "P1_d150_s%d" % sd, d_lex=150.0, seed=sd))
    for sd in (1, 2, 3):
        for mT in (9, 50):                                             # one job per amplitude
            pp.append(job(PP, "P1_d300_s%d_b%d" % (sd, mT), d_lex=300.0, seed=sd, b_list=b_of_mT(mT), n_amp=1))
            pp.append(job(PP, "P1_d300_dx3_s%d_b%d" % (sd, mT), d_lex=300.0, dx_lex=3.0, seed=sd,
                          b_list=b_of_mT(mT), n_amp=1))
    for sd in (6, 7, 8):                                               # block 2
        pp.append(job(PP, "P2_L12_s%d" % sd, seed=sd))
    for sd in range(3, 9):
        pp.append(job(PP, "P2_L30_s%d" % sd, seed=sd, Leff_lex=30.0))
        pp.append(job(PP, "P2_rp3_s%d" % sd, seed=sd, r_p=3.0))
    for sd in (3, 4, 5):                                               # block 3 (pairs with R_dx15_s3..5)
        for dx in (1.0, 2.0, 3.0):
            pp.append(job(PP, "P3_dx%g_s%d" % (dx, sd), seed=sd, dx_lex=dx))
    B7 = b_of_mT(9, 20, 35, 50, 70, 100, 150)
    PROD = dict(CV, d_lex=300.0, dx_lex=3.0, phi=0.65, b_list=B7, n_amp=len(B7))   # PLAN.md fixed settings
    freqtest = []
    for tag, kw, seeds in (("base", dict(Leff_lex=12.0, r_p=0.0), (901, 902, 903, 904)),
                           ("corner", dict(Leff_lex=30.0, r_p=3.0), (905, 906, 907, 908))):
        for sd in seeds:
            freqtest.append(job(PROD, "F_%s_f_s%d" % (tag, sd), seed=sd, **kw))
            freqtest.append(job(PROD, "F_%s_2f_s%d" % (tag, sd), seed=sd, f_rel=2.0 * F_REL_DEFAULT, **kw))
    roughtest = []
    RT = dict(PROD, d_lex=212.0, Leff_lex=12.0, r_p=0.0, b_list=b_of_mT(9, 50), n_amp=2)
    for sd in (911, 912, 913, 914):
        roughtest.append(job(RT, "R212_dx15_s%d" % sd, seed=sd, dx_lex=1.5))
        roughtest.append(job(RT, "R212_dx3_s%d" % sd, seed=sd, dx_lex=3.0))
    alphatest = []
    AT = dict(PROD, Leff_lex=12.0, r_p=0.0, b_list=b_of_mT(9, 20, 50, 100, 150), n_amp=5)
    for al in (0.01, 0.03):
        for sd in (901, 902, 903):
            alphatest.append(job(AT, "A%03.0f_s%d" % (al * 1000, sd), seed=sd, alpha=al))
    relaxtest = []
    RX = dict(PROD, Leff_lex=12.0, r_p=0.0, b_list=b_of_mT(9, 20), n_amp=2,
              relax_mode="converge", relax_maxiter=50000, relax_chunk=500, relax_tol_E=1e-6)
    for al in (0.1, 0.01):
        for sd in (901, 902, 903):
            relaxtest.append(job(RX, "X%03.0f_s%d" % (al * 1000, sd), seed=sd, alpha=al))
    smoothtest = []
    for sd in (911, 912, 913, 914):
        smoothtest.append(job(RT, "S212_dx3f_s%d" % sd, seed=sd, dx_lex=3.0, surface="fraction"))
    for sd in (1, 2, 3):
        for mT in (9, 50):
            smoothtest.append(job(PP, "S300_dx3f_s%d_b%d" % (sd, mT), d_lex=300.0, dx_lex=3.0, seed=sd,
                                  b_list=b_of_mT(mT), n_amp=1, surface="fraction"))
    # repeats after the audit (PROTOKOLL 7.13): the virgin state relaxed to convergence (hard gate)
    CONV = dict(relax_mode="converge", relax_maxiter=80000, relax_chunk=500, relax_tol_E=1e-6, relax_confirm=2)
    B0 = dict(PROD, Leff_lex=12.0, r_p=0.0, alpha=0.1, **CONV)
    baseline = [job(B0, "B300_s%d" % sd, seed=sd, snap_phases=4) for sd in range(921, 929)]   # m at 4 phases: wall measure vs B (domains.py)
    mesh300 = []
    for sd in (921, 922, 923, 924):          # same seeds as baseline: the dx 3 side also checks B300
        for dx in (1.5, 3.0):
            tag = ("%g" % dx).replace(".", "")
            M = dict(B0, dx_lex=dx, n_amp=1, cycles_per_amp=5, max_cycles_per_amp=5,
                     relax_maxiter=80000)
            first = "M300_dx%s_s%d_b9" % (tag, sd)
            mesh300.append(job(M, first, seed=sd, b_list=b_of_mT(9)))
            for mT in (50, 100):             # directly from the relaxed virgin state of the 9 mT job
                mesh300.append(job(M, "M300_dx%s_s%d_b%d" % (tag, sd, mT), seed=sd, b_list=b_of_mT(mT),
                                   init_from="%s/init.pt" % first, init_wait_h=30.0))
    # surface roughness with the converged state (PROTOKOLL 7.11 withdrawn; Chris 2026-10-06): identical to
    # the dx 3 jobs of mesh300 except the surface (volume-fraction surface instead of the staircase)
    rough300 = []
    for sd in (921, 922, 923, 924):
        M3 = dict(B0, dx_lex=3.0, n_amp=1, cycles_per_amp=5, max_cycles_per_amp=5, relax_maxiter=80000,
                  surface="fraction")
        first = "R300f_s%d_b9" % sd
        rough300.append(job(M3, first, seed=sd, b_list=b_of_mT(9)))
        for mT in (50, 100):
            rough300.append(job(M3, "R300f_s%d_b%d" % (sd, mT), seed=sd, b_list=b_of_mT(mT),
                                init_from="%s/init.pt" % first, init_wait_h=5.0))
    # checks before the Sobol design (PROTOKOLL 7.24, Chris 2026-10-04): start from the relaxed virgin
    # state of the baseline runs (the relaxation does not depend on f, alpha or the amplitudes)
    freq9 = []
    for sd in (921, 922, 923):
        for tag, fac in (("half", 0.5), ("dbl", 2.0)):
            freq9.append(job(B0, "F9_%s_s%d" % (tag, sd), seed=sd, f_rel=fac * F_REL_DEFAULT, b_list=b_of_mT(9, 20),
                             n_amp=2, init_from="B300_s%d/init.pt" % sd))
    phase16 = [job(B0, "PH_s%d" % sd, seed=sd, b_list=b_of_mT(9, 20, 150), n_amp=3, snap_phases=16, snap_cycles=3,
                   ringdown_periods=2.0, ringdown_samples=512, init_from="B300_s%d/init.pt" % sd)
               for sd in (921, 922)]
    return {"bench": bench, "meshtest": meshtest, "steady": steady, "numfloor": numfloor,
            "initproto": initproto, "demagtest": demagtest, "conv": conv, "scatter": scatter,
            "planpilot": pp, "freqtest": freqtest, "roughtest": roughtest, "alphatest": alphatest,
            "smoothtest": smoothtest, "relaxtest": relaxtest,
            "alphatest2": [job(AT, "A010c_s%d" % sd, seed=sd, alpha=0.01, **CONV) for sd in (904, 905, 906)],
            "snaptest": [job(RX, "V_s%d" % sd, seed=sd, alpha=0.1, b_list=b_of_mT(20, 150), snap_phases=4)
                         for sd in (901, 902)],
            "baseline": baseline, "mesh300": mesh300, "freq9": freq9, "phase16": phase16, "rough300": rough300,
            "pilot": pilot}


def n_cells(a):
    a_lex, _ = box_edge(a["d_lex"], a["phi"])
    N = int(round(a_lex / a["dx_lex"]))
    return N, N**3


def cost_h(a, s_per_step_Mcell, dt_tau, cycles_avg, speedup, s_min_step=0.008):
    N, n = n_cells(a)
    f_rel = a.get("f_rel", F_REL_DEFAULT)
    al = a.get("alpha", 0.02)
    if al >= 0.05:
        dtt = 1.62 * (a["dx_lex"] / 3.0) ** 1.5
    elif abs(al - 0.02) < 1e-9:
        dtt = dt_tau * (a["dx_lex"] / 3.0) ** 0.5
    else:   # estimate: dt ~ alpha^0.47 between the measured 0.02 and 0.1 (dx = 3)
        dtt = 1.62 * (al / 0.1) ** 0.47 * (a["dx_lex"] / 3.0) ** 1.5
    steps = 2.0 * math.pi / (f_rel * dtt)
    if n > 10e6:
        s_per_step_Mcell = max(s_per_step_Mcell, 0.0164)
    steps *= (1e-5 / a.get("atol", 1e-5)) ** 0.2                      # RKF45: dt ~ atol^(1/5)
    t_step = max(s_min_step, s_per_step_Mcell * n / 1e6) / speedup * (2.0 if a["precision"] == "double" else 1.0)
    n_st = len(a.get("h_list") or a.get("b_list") or []) or a["n_amp"]
    per_amp = a["cycles_per_amp"] if a["cycles_per_amp"] == a["max_cycles_per_amp"] else cycles_avg
    pre = 2 if a.get("protocol") == "virgin_asc" else 1               # relaxation / reset cycle
    if a.get("relax_mode") == "converge":
        # measured d/l_ex 300, dx 3 (relaxtest): 26 000 iterations of 1e-11 s, 0.065 s each = 1.7 RKF45
        # steps; the number of steps per iteration scales as dt (dx^1.5); 30 000 iterations assumed
        pre = 0
        if not a.get("init_from"):
            t_relax = 30000 * 1.7 * (3.0 / a["dx_lex"]) ** 1.5 * max(s_min_step, s_per_step_Mcell * n / 1e6)
            pre = t_relax / (steps * max(s_min_step, s_per_step_Mcell * n / 1e6) / speedup)
    if a.get("protocol") == "acdemag_asc":
        pre = 17 / a.get("demag_f_factor", 1.0)                       # decaying cycles (q = 0.8)
    cycles = a.get("max_cycles", pre + n_st * per_amp)
    return cycles * steps * t_step / 3600.0


def to_cli(a):
    parts = []
    for k, v in a.items():
        if v is True:
            parts.append("--%s" % k)
        elif isinstance(v, (list, tuple)):
            parts.append("--%s %s" % (k, " ".join(repr(float(x)) for x in v)))
        else:
            parts.append("--%s %s" % (k, repr(v) if isinstance(v, float) else v))
    return " ".join(parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--s_per_step_Mcell", type=float, default=0.012)
    ap.add_argument("--dt_tau", type=float, default=0.75)
    ap.add_argument("--cycles_avg", type=float, default=3.5)
    ap.add_argument("--speedup", type=float, default=1.0)
    ap.add_argument("--outdir", default=str(HERE / "jobs"))
    a = ap.parse_args()

    outdir = pathlib.Path(a.outdir)
    outdir.mkdir(exist_ok=True)
    for f in outdir.glob("*.txt"):
        f.unlink()
    total = 0.0
    for group, jobs in matrix().items():
        costed = sorted(((cost_h(j, a.s_per_step_Mcell, a.dt_tau, a.cycles_avg, a.speedup), n, j) for n, j in jobs),
                        key=lambda x: -x[0])
        with open(outdir / ("%s.txt" % group), "w") as f:
            for c, n, j in costed:
                f.write("%s %s\n" % (n, to_cli(j)))
        g = sum(c for c, _, _ in costed)
        total += g
        print("\n%-11s %2d jobs  %6.1f V100-h" % (group, len(costed), g))
        for c, n, j in costed:
            print("   %-20s N=%3d  %6.2f h" % (n, n_cells(j)[0], c))
    print("\nTOTAL %.1f V100-h (speedup %.2g)" % (total, a.speedup))


if __name__ == "__main__":
    main()
