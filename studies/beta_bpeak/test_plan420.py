"""
CPU tests of the tied cut cells of PLAN 4.20 (E6: --tie_cells, field_terms_extra.TiedLLGSolver, tie_thin_face,
tie_ghost_mirror). The numbers of the tests are those of PLAN 4.20, section Code.

    CUDA_DEVICE=-1 python -m pytest studies/beta_bpeak/test_plan420.py -q
"""
import json
import math
import os
import pathlib
import subprocess
import sys

import numpy as np
import pytest

os.environ.setdefault("CUDA_DEVICE", "-1")
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent))

import torch  # noqa: E402
from magnumnp import (Mesh, State, ExchangeField, DemagField, UniaxialAnisotropyField, LLGSolver,  # noqa: E402
                      DemagStaircaseCorrection, constants)

from anisotropy_model import box_edge  # noqa: E402
from geometry import fcc_box, fcc_subcells, FCC_SITES  # noqa: E402
from units import units  # noqa: E402
from field_terms_extra import SinusoidalDrive, TiedLLGSolver, tie_thin_face, tie_ghost_mirror  # noqa: E402

BASE_COMMIT = "8c92122"          # the diag branch before E6: the code path of D without --tie_cells


@pytest.fixture
def fp64():
    old = torch.get_default_dtype()
    torch.set_default_dtype(torch.float64)
    yield
    torch.set_default_dtype(old)


def _geo(d_lex=24.0, phi=0.3, dx_lex=3.0, q=8):
    """The fraction geometry of run_loops with the correction (sub-cells q) for a synthetic sphere packing."""
    lex = units(1.5, 1e-11)["l_ex"]
    a_lex, d_lex = box_edge(d_lex, phi, 12.0)
    N = int(round(a_lex / dx_lex))
    d, dx = d_lex * lex, dx_lex * lex
    g = fcc_box(N, d, dx, sub=q)
    n_s, pid_s, occ_s = fcc_subcells(N, d, dx, q)
    assert np.array_equal(pid_s, g["ids"].astype(np.int64)) and np.allclose(n_s / q**3, g["frac"], atol=1e-6)
    return dict(N=N, d=d, dx=dx, dx_lex=dx_lex, q=q, ids=g["ids"], frac=g["frac"], n=n_s, occ=occ_s)


def _tie(g, kind):
    if kind == "thin_face":
        return tie_thin_face(g["ids"], g["frac"], g["dx_lex"], g["occ"], g["q"])
    return tie_ghost_mirror(g["ids"], g["frac"], g["N"], g["d"], g["dx"])


def _state_terms(g, seed=1, h_hold=0.05):
    """State and the terms of a D run (Newell stray field with true_pbc, the pair correction q 8 p 3 with the centroid,
    exchange, cube anisotropy with random axes, a constant drive), random m."""
    U = units(1.5, 1e-11)
    Ms, N = U["Ms"], g["N"]
    rng = np.random.default_rng(seed)
    mag4 = (g["ids"] >= 0)[..., None]
    f4 = g["frac"].astype(np.float64)[..., None]
    st = State(Mesh((N, N, N), (g["dx"],) * 3, pbc=(4, 4, 4)))
    T = lambda a: torch.tensor(a, dtype=st.dtype)
    st.material = {"alpha": 1.0}
    st.material["Ms"] = T(np.where(mag4, np.maximum(f4 * Ms, 1e-8 * Ms), 1e-8 * Ms))
    st.material["A"] = T(np.where(mag4, f4 * 1e-11, 0.0))
    ax = rng.standard_normal((N, N, N, 3))
    st.material["Ku"] = T(np.where(mag4, f4 * U["Kd"] / 16.0, 0.0))
    st.material["Ku_axis"] = T(ax / np.linalg.norm(ax, axis=-1, keepdims=True))
    demag = DemagField(true_pbc=True, extrapolate=True)
    corr = DemagStaircaseCorrection(g["n"], g["occ"], g["q"], 3, Ms, form="pair", centroid=True, demag=demag)
    drive = SinusoidalDrive(st, (1.0, 0.0, 0.0))
    drive.H_amp, drive.hold = h_hold * Ms, 1.0
    terms = [demag, corr, ExchangeField(), UniaxialAnisotropyField(), drive]
    m = rng.standard_normal((N, N, N, 3))
    st.m = T(m / np.linalg.norm(m, axis=-1, keepdims=True))
    return st, terms


# --- test 5: the tie maps ---------------------------------------------------------------------------------------------

def test_tie_maps_synthetic():
    """(5) V1: each thin cell is tied to the face neighbour of its particle with the largest shared occupied face (a
    direct count over the sub-cells), or to the end of its chain, or (fall-back) to a free cell; the master is never
    thin. V2: the tied cells are the cut cells with the centre outside (direct count); the mirror point lies inside;
    the masters have the centre inside and the same particle; the weights sum to 1. V1 is a subset of V2."""
    g = _geo()
    ids, frac, occ, q, N = g["ids"], g["frac"], g["occ"], g["q"], g["N"]
    v1, v2 = _tie(g, "thin_face"), _tie(g, "ghost_mirror")
    cut = (ids >= 0) & (frac < 1.0)
    thin = cut & (frac < 1.0 / g["dx_lex"])
    assert v1.info["n_thin"] == int(thin.sum()) > 0 and v1.info["n_cut"] == int(cut.sum())
    o = occ.reshape(N, N, N, q, q, q)
    unflat = lambda j: (j // (N * N), (j // N) % N, j % N)
    nxt = {}
    for j in v1.tied:
        c = unflat(int(j))
        best, bn = None, 0
        for ax, side in ((0, 1), (0, -1), (1, 1), (1, -1), (2, 1), (2, -1)):
            nb = list(c)
            nb[ax] = (nb[ax] + side) % N
            nb = tuple(nb)
            if ids[nb] != ids[c]:
                continue
            n = 0
            for u in range(q):
                for v in range(q):
                    s_own = [u, v]
                    s_own.insert(ax, q - 1 if side > 0 else 0)
                    s_nb = [u, v]
                    s_nb.insert(ax, 0 if side > 0 else q - 1)
                    n += bool(o[c][tuple(s_own)] and o[nb][tuple(s_nb)])
            if n > bn:
                best, bn = nb, n
        nxt[c] = best
    n_direct = 0
    for j, k in zip(v1.tied, v1.src[:, 0]):
        c, mk = unflat(int(j)), unflat(int(k))
        assert not thin[mk] and ids[mk] == ids[c]
        e, chain = nxt[c], 1
        while e is not None and thin[e] and chain <= 3:
            e, chain = nxt[e], chain + 1
        if e is not None and not thin[e] and chain <= 3:
            assert mk == e
            n_direct += 1
    assert n_direct == v1.info["n_thin"] - v1.info["n_fallback_no_face"] - v1.info["n_fallback_chain"]
    # V2
    a, R = N * g["dx"], 0.5 * g["d"]
    x = (np.arange(N) + 0.5) * g["dx"]
    outside = np.zeros_like(cut)
    for i, j, k in np.argwhere(ids >= 0):
        s = FCC_SITES[ids[i, j, k]] * a
        dd = np.array([x[i], x[j], x[k]]) - s
        dd -= a * np.round(dd / a)
        outside[i, j, k] = (dd**2).sum() > R**2
    assert v2.info["n_centre_outside"] == int((cut & outside).sum()) == len(v2.tied)
    assert set(v1.tied.tolist()) <= set(v2.tied.tolist())
    for j, src, w in zip(v2.tied, v2.src, v2.w):
        c = unflat(int(j))
        assert w.sum() == pytest.approx(1.0, abs=1e-12) and (w >= 0).all()
        for k in src[w > 0]:
            mk = unflat(int(k))
            assert ids[mk] == ids[c] and not outside[mk]
        s = FCC_SITES[ids[c]] * a
        dd = np.array([x[c[0]], x[c[1]], x[c[2]]]) - s
        dd -= a * np.round(dd / a)
        r = math.sqrt((dd**2).sum())
        assert 2.0 * R - r < R


def test_tie_counts_staircase_geometry():
    """(5) The counts of PLAN 4.20 at 48^3 (staircase geometry, q 8): 15584 cut cells, 7104 thin (V1), 8256 with the
    centre outside (V2)."""
    g = _geo(d_lex=96.0, phi=0.55, dx_lex=3.0, q=8)
    assert g["N"] == 48
    v1, v2 = _tie(g, "thin_face"), _tie(g, "ghost_mirror")
    assert (v1.info["n_cut"], v1.info["n_thin"], v2.info["n_centre_outside"]) == (15584, 7104, 8256)


# --- tests 2 ... 4: the solver ----------------------------------------------------------------------------------------

@pytest.mark.parametrize("kind", ["thin_face", "ghost_mirror"])
def test_tie_fixed_point_after_steps(fp64, kind):
    """(2) After 100 steps, each tied cell equals its tie value (bitwise: the step normalises, then ties); dm of the
    tied cells is 0."""
    g = _geo()
    st, terms = _state_terms(g)
    tie = _tie(g, kind)
    llg = TiedLLGSolver(terms, tie, st.material["Ms"], atol=1e-6)
    llg.tie_state(st)
    for _ in range(100):
        llg.step(st, 2e-13)
    m = st.m.clone()
    assert torch.equal(llg.tie(m)[0], m)
    d = llg.dm(st.t, m, state=st, alpha=1.0).reshape(-1, 3)
    assert torch.equal(d[torch.as_tensor(tie.tied)], torch.zeros(len(tie.tied), 3))
    assert float(d.abs().max()) > 0


@pytest.mark.parametrize("kind", ["thin_face", "ghost_mirror"])
def test_master_field_is_energy_gradient(fp64, kind):
    """(3) D_k h_k^* equals -(1/mu0 V) dE/dm_k of the tied energy (tangential part, relative 1e-6, fp64): the numerical
    gradient perturbs the free m_k, ties again and evaluates E. Control: a free cell that is no master (M_s h)."""
    g = _geo()
    st, terms = _state_terms(g, seed=2)
    tie = _tie(g, kind)
    llg = TiedLLGSolver(terms, tie, st.material["Ms"], atol=1e-6)
    llg.tie_state(st)
    x0 = st.m.clone()
    h, w = llg.field_p(st)
    V = g["dx"] ** 3
    masters, rows, wts = tie.master_table()
    pick = list(masters[np.argsort(-(wts > 0).sum(1), kind="stable")[:4]])      # masters with the most tied cells
    free = np.setdiff1d(np.flatnonzero((g["ids"] >= 0).reshape(-1)), np.concatenate([masters, tie.tied]))
    pick.append(int(free[len(free) // 2]))
    eps = 1e-5

    def E(x):
        st.m = llg.tie(x)[0]
        return float(llg.E(st))

    for k in pick:
        grad = np.zeros(3)
        for c in range(3):
            xp, xm = x0.clone(), x0.clone()
            xp.view(-1, 3)[k, c] += eps
            xm.view(-1, 3)[k, c] -= eps
            grad[c] = (E(xp) - E(xm)) / (2 * eps)
        mk = x0.view(-1, 3)[k].numpy()
        P = np.eye(3) - np.outer(mk, mk)
        num = P @ (-grad / (constants.mu_0 * V))
        ana = P @ (float(w.view(-1)[k]) * h.view(-1, 3)[k].numpy())
        assert np.linalg.norm(num - ana) <= 1e-6 * np.linalg.norm(ana), (k, num, ana)


def _energies(llg, st, n, dt):
    E = [float(llg.E(st))]
    for _ in range(n):
        llg.step(st, dt)
        E.append(float(llg.E(st)))
    return np.array(E)


def _p_box(llg, st, tied, a, V):
    """Dissipated power of the box (W) as run_loops.diagnostics: p_coef w |m x H|^2, summed over the cells times V."""
    p_coef = a * constants.gamma * constants.mu_0 / (1.0 + a**2)
    if tied:
        H, w = llg.field_p(st)
    else:
        H, w = sum(t.h(st) for t in llg._terms), st.material["Ms"]
    mxH = torch.linalg.cross(st.m, H)
    return float((p_coef * w * (mxH * mxH).sum(-1, keepdim=True)).sum()) * V


@pytest.mark.parametrize("kind", [None, "thin_face", "ghost_mirror"])
def test_dissipation_equals_energy_loss(fp64, kind):
    """(6a) The dissipation p_dis of run_loops (with --tie_cells: h^* and the group moment D_k, 0 at the tied cells)
    equals the loss of energy per time, -dE/dt, at alpha 0.1 with precession and a constant drive field: 5 steps of
    1e-15 s, relative difference at most 1e-3 (fp64; the same check without the tie)."""
    g = _geo()
    st, terms = _state_terms(g, seed=4)
    st.material["alpha"] = 0.1
    if kind is None:
        llg = LLGSolver(terms, atol=1e-6)
    else:
        llg = TiedLLGSolver(terms, _tie(g, kind), st.material["Ms"], atol=1e-6)
        llg.tie_state(st)
    for _ in range(5):
        E0, P0, t0 = float(llg.E(st)), _p_box(llg, st, kind is not None, 0.1, g["dx"] ** 3), float(st.t)
        llg.step(st, 1e-15)
        E1, P1, t1 = float(llg.E(st)), _p_box(llg, st, kind is not None, 0.1, g["dx"] ** 3), float(st.t)
        rate = -(E1 - E0) / (t1 - t0)
        assert abs(rate - 0.5 * (P0 + P1)) <= 1e-3 * 0.5 * (P0 + P1), (kind, rate, P0, P1)


@pytest.mark.parametrize("kind", ["thin_face", "ghost_mirror"])
def test_energy_decreases(fp64, kind):
    """(4) At alpha 1, fp64, from a random start: the energy decreases in each of 200 steps. The same 200 steps without
    the tie are the reference for the size of a rise (the normalisation after each step is a projection)."""
    g = _geo()
    st, terms = _state_terms(g, seed=3)
    E0 = _energies(LLGSolver(terms, atol=1e-6), st, 200, 1e-13)
    up0 = max(0.0, float(np.max(np.diff(E0))))
    st, terms = _state_terms(g, seed=3)
    llg = TiedLLGSolver(terms, _tie(g, kind), st.material["Ms"], atol=1e-6)
    llg.tie_state(st)
    E1 = _energies(llg, st, 200, 1e-13)
    assert E1[-1] < E1[0]
    assert float(np.max(np.diff(E1))) <= max(up0, 1e-12 * abs(E1[0])), (np.max(np.diff(E1)), up0)


# --- tests 1 and 6: run_loops -----------------------------------------------------------------------------------------

def _run(script_dir, args, out, timeout=1800, device="-1"):
    env = dict(os.environ, CUDA_DEVICE=device, PYTHONPATH=str(HERE.parent.parent))
    return subprocess.run([sys.executable, str(script_dir / "run_loops.py"), "--out", str(out)] + args,
                          cwd=script_dir, env=env, capture_output=True, text=True, timeout=timeout)


BASE = ("--d_lex 24 --phi 0.3 --Leff_lex 12 --dx_lex 3 --protocol virgin_asc --relax_mode converge "
        "--relax_chunk 100 --relax_tol_E 1e-4 --relax_confirm 1 --seed 3 --b_list 0.066667 --n_amp 1 --samples 64 "
        "--alpha 0.1 --strayfield newell_ext --strayfield_shells 4 --atol 1e-06 --atol_cycles 1e-07 "
        "--rtol_cycles 1e-06").split()
D_ARGS = "--surface fraction --staircase pair --staircase_q 8 --staircase_p 3 --staircase_centroid".split()


def _rows(path):
    return [x for x in path.read_text().splitlines() if not x.startswith("#")]


def test_untied_path_equals_base_commit(tmp_path):
    """(1) Without --tie_cells, run_loops gives the outputs of the commit before E6 bitwise: a D run (fp32, the
    correction q 8) from a transferred voxel state, its relaxation and 1 cycle (samples with p_dis, Bloch points)."""
    old = tmp_path / "old"
    old.mkdir()
    repo = HERE.parent.parent
    tar = subprocess.run(["git", "-C", str(repo), "archive", BASE_COMMIT, "studies/beta_bpeak"], capture_output=True)
    assert tar.returncode == 0
    subprocess.run(["tar", "-x", "-C", str(old)], input=tar.stdout, check=True)
    old_dir = old / "studies" / "beta_bpeak"
    src = tmp_path / "vox"
    r = _run(HERE, BASE + ["--precision", "single", "--surface", "voxel", "--relax_maxiter", "5000", "--prep_only"], src)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    args = BASE + D_ARGS + ["--precision", "single", "--relax_maxiter", "3000", "--init_transfer", str(src / "init.pt"),
                            "--max_cycles", "1", "--diag_samples"]
    outs = {}
    for tag, sd in (("new", HERE), ("old", old_dir)):
        outs[tag] = tmp_path / tag / "run"
        r = _run(sd, args, outs[tag])
        assert r.returncode == 0, tag + r.stdout[-2000:] + r.stderr[-2000:]
    a, b = (torch.load(str(outs[t] / "init.pt"), weights_only=False) for t in ("new", "old"))
    assert torch.equal(a["m"], b["m"])
    for f in ("samples_00_amp00.csv", "bloch_00_amp00.csv", "diag_00_amp00.csv", "bloch_pos_00_amp00.csv"):
        assert _rows(outs["new"] / f) == _rows(outs["old"] / f), f


@pytest.mark.parametrize("device", ["-1", pytest.param("0", marks=pytest.mark.skipif(
    not torch.cuda.is_available(), reason="no GPU"))])
def test_run_loops_tie_end_to_end(tmp_path, device):
    """(6) run_loops with --tie_cells (V1, V2) on the CPU (device -1) and on a GPU (device 0, box test B5): a D run from
    a transferred voxel state of the same mesh (as F8 and stage C), relaxation and 1 cycle with --diag_samples. The
    record has the tie counts; the relaxed state and the state after the cycle are fixed points of the tie (bitwise);
    on the GPU, a second run of V2 gives the same samples bitwise. The energy balance of the cycle is not a criterion
    here: the first cycle after the relaxation is not closed (closure 0.35 for D, |w_dis/w_loop - 1| 1.7 for D on
    the CPU); test 6a checks p_dis against -dE/dt directly. The balances are printed."""
    src = tmp_path / "vox"
    r = _run(HERE, BASE + ["--precision", "single", "--surface", "voxel", "--relax_maxiter", "5000", "--prep_only"],
             src, device=device)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    common = BASE + D_ARGS + ["--precision", "single", "--relax_maxiter", "5000", "--init_transfer",
                              str(src / "init.pt"), "--max_cycles", "1", "--diag_samples"]
    g = _geo()
    bal = {}
    for kind in (None, "thin_face", "ghost_mirror"):
        out = tmp_path / (kind or "D")
        r = _run(HERE, common + (["--tie_cells", kind] if kind else []), out, device=device)
        assert r.returncode == 0, str(kind) + r.stdout[-2000:] + r.stderr[-2000:]
        s = json.loads((out / "summary.json").read_text())
        c = s["stages"][0]["cycles"][0]
        bal[kind] = abs(c["w_dis"] / c["w_loop"] - 1.0)
        if kind is None:
            continue
        info = json.loads((out / "config.json").read_text())["diag"]["tie"]
        tie = _tie(g, kind)
        assert info["kind"] == kind and info["n_tied"] == len(tie.tied) > 0
        tied = torch.as_tensor(tie.tied)
        for f in ("init.pt", "checkpoint.pt"):
            if not (out / f).exists():
                continue
            mf = torch.load(str(out / f), weights_only=False)["m"].float()
            llg = TiedLLGSolver([], tie, torch.ones(g["N"], g["N"], g["N"], 1, dtype=torch.float32), atol=1e-6)
            mt = llg.tie(mf)[0].reshape(-1, 3)[tied]
            if device == "-1":
                assert torch.equal(mt, mf.reshape(-1, 3)[tied]), (kind, f)
            else:                                   # the tie on the CPU; the run tied on the GPU
                assert float((mt - mf.reshape(-1, 3)[tied]).abs().max()) < 1e-6, (kind, f)
    print("energy balance |w_dis/w_loop - 1| of the first cycle:", bal)
    if device != "-1":
        out2 = tmp_path / "ghost_mirror_2"
        r = _run(HERE, common + ["--tie_cells", "ghost_mirror"], out2, device=device)
        assert r.returncode == 0
        assert _rows(out2 / "samples_00_amp00.csv") == _rows(tmp_path / "ghost_mirror" / "samples_00_amp00.csv")


# --- eval_420: the rules on synthetic values ---------------------------------------------------------------------------

def test_eval_420_gate_A_rules():
    """Gate A: a variant passes with a converged relaxation, no pair above 90 deg and no Bloch point in at least 3 of
    the 4 seeds; a missing seed counts as failed."""
    from eval_420 import gate_A

    def r(n90=0, bloch=(0, 0), conv=True):
        return {"relax_converged": conv, "pairs_gt90": {"all": n90}, "bloch": list(bloch)}
    R = {"V1": {941: r(), 942: r(), 943: r(), 944: r(n90=1)},
         "V2": {941: r(), 942: r(bloch=(1, 0)), 943: r(conv=False)}}
    g = gate_A(R)
    assert g["V1"] == {"seeds_ok": [941, 942, 943], "passes": True}
    assert g["V2"] == {"seeds_ok": [941], "passes": False}
    R["V2"] = {941: r(), 942: r(), 943: r()}                       # 944 missing
    assert gate_A(R)["V2"]["passes"] is True
    R["V2"] = {941: r(), 942: r(), 943: r(bloch=(0, 1)), 944: r(n90=3)}
    assert gate_A(R)["V2"]["passes"] is False


def _job(bloch_relaxed=False, samples=0, ran=True, w=1e-5):
    return {"bloch_relaxed": bloch_relaxed, "bloch_samples": samples, "ran_cycles": ran, "w": w}


def test_eval_420_stage_B_rules():
    """B-0, B-a and the selection rule on synthetic values."""
    from eval_420 import b_0, b_a, select
    J = {"F8": {941: _job(), 942: _job(), 943: _job(w=None), 944: None},
         "D8": {941: _job(samples=3), 942: _job(samples=1), 943: _job(), 944: _job()},
         "D8V1": {941: _job(), 942: _job(), 943: _job(ran=False), 944: _job()},
         "D8V2": {941: _job(), 942: _job(samples=1), 943: _job(), 944: _job()}}
    assert b_0(J)["valid"] is False                                 # values in 2 seeds only
    J["F8"][944] = _job()
    assert b_0(J)["valid"] is True
    J["F8"][944] = _job(bloch_relaxed=True)
    assert b_0(J)["valid"] is False                                 # a Bloch point in F8
    ba = b_a(J, ["V1", "V2"])
    assert ba["valid"] and ba["variants"]["V1"]["passes"] and not ba["variants"]["V2"]["passes"]
    assert ba["variants"]["V1"]["seeds_ran"] == [941, 942, 944]
    assert select(ba, None, ["V1", "V2"])[0] == "V1"
    J["D8V1"][941] = _job(samples=2)
    assert select(b_a(J, ["V1", "V2"]), None, ["V1", "V2"])[0] == "NONE"
    J["D8"] = {s: _job() for s in (941, 942, 943, 944)}             # no Bloch in D8: B-a not valid
    ba = b_a(J, ["V1", "V2"])
    assert not ba["valid"]
    ov = lambda xs: {s: {"cut_overlap_mean": x} for s, x in zip((941, 942, 943, 944), xs)}
    bb = {"D8V1": ov([0.90, 0.91, 0.92, 0.93]), "D8V2": ov([0.95, 0.96, 0.97, 0.98])}
    assert select(ba, bb, ["V1", "V2"])[0] == "V2"                  # a clear paired difference
    bb["D8V2"] = ov([0.91, 0.90, 0.93, 0.92])                       # difference smaller than its SE: V1
    assert select(ba, bb, ["V1", "V2"])[0] == "V1"
    assert select(ba, bb, ["V2"])[0] == "V2"                        # only V2 passed gate A


def test_eval_420_stage_C_rules():
    """C-a and the verdict on synthetic values: the counted jobs, a Bloch point in a relaxed state fails the seed, the
    7-of-9 rule, the spare seed, and the q-16 positive control."""
    from eval_420 import c_a, verdict_C
    D = {(s, a): {"bloch_relaxed": s == 922} for s in (921, 922, 923, 924) for a in (9, 50, 100)}
    Dv = {(s, a): _job() for s in (921, 922, 923, 924) for a in (9, 50, 100)}
    ca = c_a(D, Dv)
    assert len(ca["q8"]["counted"]) == 9 and ca["q8"]["status"] == "holds"
    assert verdict_C(ca) == "no Bloch points in these runs"
    Dv[(922, 100)] = _job(samples=5)                                # seed 922 is not counted
    assert c_a(D, Dv)["q8"]["status"] == "holds"
    Dv[(923, 9)] = _job(bloch_relaxed=True)                         # fails all jobs of seed 923
    ca = c_a(D, Dv)
    assert sorted(map(tuple, ca["q8"]["fail"])) == [(923, 9), (923, 50), (923, 100)]
    assert verdict_C(ca) == "remedy fails"
    Dv[(923, 9)] = _job()
    for k in ((921, 9), (921, 50), (921, 100)):
        Dv[k] = _job(ran=False)
    ca = c_a(D, Dv)
    assert ca["q8"]["status"] == "not determined" and verdict_C(ca) == "not determined"
    sp = {"Dv": {(925, a): _job() for a in (9, 50, 100)},
          "D": {(925, 9): {"bloch_relaxed": False}, (925, 50): {"bloch_relaxed": False}}}
    ca = c_a(D, Dv, spare=sp)                                       # 925 at 100 mT: no D record, not counted
    assert len(ca["q8"]["ran"]) == 8 and ca["q8"]["status"] == "holds"
    D16 = {(s, 100): _job(samples=2 if s != 924 else 0) for s in (921, 923, 924)}
    Dv16 = {(s, 100): _job() for s in (921, 923, 924)}
    ca = c_a(D, Dv, D16, Dv16, sp)
    assert ca["q16"]["valid"] and ca["q16"]["status"] == "holds"
    Dv16[(921, 100)] = _job(samples=1)
    assert verdict_C(c_a(D, Dv, D16, Dv16, sp)) == "remedy fails"
    D16 = {(s, 100): _job() for s in (921, 923, 924)}
    assert c_a(D, Dv, D16, Dv16, sp)["q16"]["status"] == "not valid"
