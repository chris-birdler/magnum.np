"""
CPU tests of the PLAN 4.18 diagnostics (E1 ... E4 of run_loops, field_terms_extra.CutPairExchange, eval_418).

    CUDA_DEVICE=-1 python -m pytest studies/beta_bpeak/test_plan418.py -q
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
from magnumnp import Mesh, State, ExchangeField, DemagFieldPBC, LLGSolver, constants  # noqa: E402

STAIRCASE_RUNS = pathlib.Path("/home/christoph/Documents/code/magnum.np-staircase/studies/staircase_correction/"
                              "remote_20261021_box6/runs")


@pytest.fixture
def fp64():
    old = torch.get_default_dtype()
    torch.set_default_dtype(torch.float64)
    yield
    torch.set_default_dtype(old)


def _surface_state(n=8, seed=1):
    """A periodic state with a ball of full cells, a shell of cut cells (0 < f < 1) and void, random m."""
    g = np.random.default_rng(seed)
    x = np.arange(n) - (n - 1) / 2.0
    r = np.sqrt(sum(c ** 2 for c in np.meshgrid(x, x, x, indexing="ij")))
    f = np.clip(3.0 - r, 0.0, 1.0)
    f[(f > 0) & (f < 1)] = g.uniform(0.05, 0.95, int(((f > 0) & (f < 1)).sum()))
    Ms0, A0 = 1.0 / constants.mu_0, 1e-11
    st = State(Mesh((n, n, n), (1e-9,) * 3, pbc=(1, 1, 1)))
    st.material = {"Ms": torch.tensor(np.maximum(f, 1e-8)[..., None] * Ms0),
                   "A": torch.tensor(f[..., None] * A0), "alpha": 1.0}
    m = g.standard_normal((n, n, n, 3))
    st.m = torch.tensor(m / np.linalg.norm(m, axis=-1, keepdims=True))
    return st, f


def test_cut_pair_exchange_k1_equals_exchange_field(fp64):
    """E1: with K = 1 the field equals that of ExchangeField of magnumnp (relative difference below 1e-6, fp64)."""
    from field_terms_extra import CutPairExchange
    st, f = _surface_state()
    h0 = ExchangeField().h(st)
    h1 = CutPairExchange((f > 0) & (f < 1), 1.0).h(st)
    assert float((h1 - h0).norm() / h0.norm()) < 1e-6


def test_cut_pair_exchange_count_energy_and_locality(fp64):
    """E1: the number of changed pairs equals a direct count; the energy is linear in K (a sum over pairs: E(4) - E(1)
    = 3 (E(2) - E(1))); h(K) - h(1) is zero away from the cut cells and their neighbours; a mesh that is not periodic
    stops."""
    from field_terms_extra import CutPairExchange
    st, f = _surface_state()
    cut = (f > 0) & (f < 1)
    n = f.shape[0]
    direct = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
                    q = ((i + d[0]) % n, (j + d[1]) % n, (k + d[2]) % n)
                    if f[i, j, k] > 0 and f[q] > 0 and (cut[i, j, k] or cut[q]):
                        direct += 1
    t = {K: CutPairExchange(cut, K) for K in (1.0, 2.0, 4.0)}
    assert t[4.0].n_pairs_changed(st.material["A"]) == direct > 0
    E = {K: float(tm.E(st)) for K, tm in t.items()}
    assert E[4.0] - E[1.0] == pytest.approx(3.0 * (E[2.0] - E[1.0]), rel=1e-9)
    dh = (t[4.0].h(st) - t[1.0].h(st)).norm(dim=-1).numpy()
    near = cut.copy()
    for ax in range(3):
        for sh in (1, -1):
            near |= np.roll(cut, sh, axis=ax)
    assert dh[~near].max() == 0.0 and dh[cut].max() > 0.0
    st2 = State(Mesh((4, 4, 4), (1e-9,) * 3))
    st2.material = {"Ms": 1.0, "A": 1e-11}
    st2.m = torch.ones(4, 4, 4, 3) / math.sqrt(3)
    with pytest.raises(ValueError):
        CutPairExchange(np.zeros((4, 4, 4), bool), 2.0).h(st2)


def test_term_fields_sum_is_effective_field(fp64):
    """E2: the sum of the term fields of run_loops.term_fields gives the right-hand side of the solver (alpha 1)."""
    from field_terms_extra import CutPairExchange
    from run_loops import term_fields
    st, f = _surface_state()
    terms = [DemagFieldPBC(), CutPairExchange((f > 0) & (f < 1), 4.0)]
    llg = LLGSolver(terms)
    h = sum(term_fields(st, terms, ["demag", "exchange"]).values())
    m = st.m
    gp = constants.gamma / 2.0
    rhs = -gp * torch.linalg.cross(m, torch.linalg.cross(m, h)) - gp * torch.linalg.cross(m, h)
    dm = llg.dm(st.t, m, state=st, alpha=1.0)
    assert float((dm - rhs).norm() / rhs.norm()) < 1e-10


def test_aligned_m_sets_cut_cells_along_full_cell_sum():
    """E2: aligned_m turns each cut cell to the sum of the full cells of its particle within 2 cells and keeps the
    full cells; a cut cell without a full cell in the block keeps its m and is counted."""
    from run_loops import aligned_m, full_cell_sum
    n = 10
    ids = -np.ones((n, n, n), np.int64)
    f = np.zeros((n, n, n))
    ids[2:7, 2:7, 2:7] = 0
    f[2:7, 2:7, 2:7] = 1.0
    f[2, 2:7, 2:7] = 0.4                                  # a cut layer
    ids[9, 9, 9], f[9, 9, 9] = 0, 0.3                     # a cut cell far from all full cells
    m = np.zeros((n, n, n, 3))
    m[..., 2] = 1.0
    m[2, 2:7, 2:7] = [1.0, 0.0, 0.0]                      # cut layer tilted by 90 deg
    m[9, 9, 9] = [0.0, 1.0, 0.0]
    ma, n_noref = aligned_m(m, ids, f)
    assert n_noref == 1
    assert np.allclose(ma[2, 2:7, 2:7], [0.0, 0.0, 1.0])
    assert np.allclose(ma[9, 9, 9], [0.0, 1.0, 0.0]) and np.allclose(ma[3:7, 2:7, 2:7], m[3:7, 2:7, 2:7])
    s = full_cell_sum(m, ids, f)
    assert s[2, 4, 4, 2] == 2 * 5 * 5                     # 2 full layers x 5 x 5 cells in the block


def test_bloch_rows_synthetic_hedgehog():
    """E3: a hedgehog in a ball gives one row at its cube, with the charge +1, its particle and the largest angle in
    the cube (acos(1/3) = 70.53 deg between neighbour corners)."""
    from run_loops import bloch_rows_t, bloch_mask_t
    n = 12
    x = np.arange(n) - 5.5
    X = np.stack(np.meshgrid(x, x, x, indexing="ij"), -1)
    r = np.linalg.norm(X, axis=-1)
    ids = np.where(r < 4.5, 2, -1)
    m = torch.tensor(X / r[..., None])
    rows = bloch_rows_t(m, torch.tensor(ids), bloch_mask_t(torch.tensor(ids)))
    assert len(rows) == 1
    i, j, k, Q, p, ang = rows[0]
    assert (i, j, k, Q, p) == (5, 5, 5, 1, 2)
    assert ang == pytest.approx(math.degrees(math.acos(1.0 / 3.0)), abs=1e-6)


def test_large_pair_120_counted_once():
    """E3: a pair of one particle at 120 deg is counted once; a pair at 60 deg and a pair of two particles are not."""
    from run_loops import large_pairs_t
    n = 6
    ids = -np.ones((n, n, n), np.int64)
    ids[1, 1, 1] = ids[2, 1, 1] = 0                       # one particle of two cells
    ids[1, 4, 1], ids[2, 4, 1] = 0, 1                     # two particles side by side
    ids[4, 1, 1] = ids[4, 2, 1] = 3
    m = np.zeros((n, n, n, 3))
    m[..., 0] = 1.0
    c, s = math.cos(math.radians(120)), math.sin(math.radians(120))
    m[2, 1, 1] = [c, s, 0.0]
    m[2, 4, 1] = [-1.0, 0.0, 0.0]
    m[4, 2, 1] = [0.5, math.sqrt(0.75), 0.0]              # 60 deg
    idt = torch.tensor(ids)
    masks = [(idt >= 0) & (idt == torch.roll(idt, 1, dims=ax)) for ax in range(3)]
    n90, cmin = large_pairs_t(torch.tensor(m), masks)
    assert n90 == 1 and cmin == pytest.approx(c)


def test_diag_init_error_rules():
    """E4: a state of the same mesh and geometry is accepted (also with another correction on the same cells); another
    argument, another mesh, another geometry (other sub-cells) or a relaxation that did not converge stop the run."""
    from run_loops import diag_init_error, locked_geometry, INIT_KEYS
    N, d, dx = 12, 5.0e-8, 1.0e-8
    geo_of = lambda sub: locked_geometry(N, d, dx, 1, sub=sub)
    cfg = {k: None for k in INIT_KEYS}
    cfg.update(surface="fraction", staircase="none", staircase_q=2, diag_k_cut=4.0, relax_maxiter=80000)
    geo = geo_of(8)
    src = dict(cfg, staircase="pair", staircase_q=8, staircase_p=3, staircase_centroid=True, diag_k_cut=None,
               relax_maxiter=1000)
    c = {"m": torch.zeros(N, N, N, 3), "init": {"relax_converged": True}, "config": src}
    assert diag_init_error(c, cfg, N, geo, geo_of) is None
    assert "other geometry" in diag_init_error(dict(c, config=dict(src, staircase_q=4)), cfg, N, geo, geo_of)
    assert "other arguments" in diag_init_error(dict(c, config=dict(src, seed=5)), cfg, N, geo, geo_of)
    assert "shape" in diag_init_error(dict(c, m=torch.zeros(6, 6, 6, 3)), cfg, N, geo, geo_of)
    assert "converge" in diag_init_error(dict(c, init={"relax_converged": False}), cfg, N, geo, geo_of)


def _run(args, out, timeout=900):
    env = dict(os.environ, CUDA_DEVICE="-1")
    return subprocess.run([sys.executable, str(HERE / "run_loops.py"), "--out", str(out)] + args, cwd=HERE, env=env,
                          capture_output=True, text=True, timeout=timeout)


BASE = ("--d_lex 24 --phi 0.3 --Leff_lex 12 --dx_lex 3 --surface fraction --protocol virgin_asc --relax_mode converge "
        "--relax_chunk 100 --relax_tol_E 1e-4 --relax_confirm 1 --seed 3 --b_list 0.006 --n_amp 1 --samples 8 "
        "--precision double").split()


def test_run_loops_diag_end_to_end(tmp_path):
    """E1, E2, E4 and E3 in run_loops on a small fraction-surface mesh: (A) a relaxed state; (B) the start from it with
    --diag_init_state, --diag_k_cut 4 and --diag_fields and no relaxation step: the start m equals the file, the record
    names the source, diag_fields.npz holds both states and all term fields; (C) a source of another mesh stops the
    run; (D) one cycle with --diag_samples writes the per-sample file (one row per sample) and the Bloch-point file."""
    a, b, c, d = (tmp_path / x for x in "abcd")
    r = _run(BASE + ["--relax_maxiter", "5000", "--prep_only"], a)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    r = _run(BASE + ["--relax_maxiter", "0", "--allow_unconverged", "--prep_only", "--diag_k_cut", "4",
                     "--diag_fields", "--diag_init_state", str(a / "init.pt")], b)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    ca, cb = (torch.load(str(x / "init.pt"), weights_only=False) for x in (a, b))
    assert torch.equal(ca["m"], cb["m"])
    assert cb["init"]["diag_init_state"]["source"] == str(a / "init.pt")
    assert json.loads((b / "config.json").read_text())["diag"]["k_cut"]["n_pairs_changed"] > 0
    z = np.load(b / "diag_fields.npz")
    names = [str(x) for x in z["term_names"]]
    assert names[:2] == ["demag", "exchange"]
    for tag in ("relaxed", "aligned"):
        for nm in names:
            assert z["h_%s_%s" % (tag, nm)].shape == z["m_relaxed"].shape
    assert np.allclose(z["m_relaxed"], ca["m"].numpy(), atol=1e-6)
    other = ["--d_lex", "36"] + BASE[2:]                  # another box, thus another mesh
    r = _run(other + ["--relax_maxiter", "0", "--prep_only", "--diag_init_state", str(a / "init.pt")], c)
    assert r.returncode != 0 and "--diag_init_state" in (r.stdout + r.stderr)
    r = _run(BASE + ["--relax_maxiter", "5000", "--init_from", str(a / "init.pt"), "--max_cycles", "1",
                     "--diag_samples"], d)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    rows = [x for x in (d / "diag_00_amp00.csv").read_text().splitlines() if not x.startswith("#")]
    assert len(rows) == 9 and all(len(x.split(",")) == 4 for x in rows)
    assert (d / "bloch_pos_00_amp00.csv").read_text().startswith("# cycle,sample,i,j,k,Q,particle")


@pytest.mark.skipif(not (STAIRCASE_RUNS / "T3Pc8_s941" / "init.pt").exists(), reason="GPU states of the staircase study")
def test_eval_418_reproduces_basis_table():
    """eval_418 gives the tilt table of the PLAN 4.18 basis for the GPU states T3Pc8_s941 and s944 (plan418_basis
    tilt_ref.py and pairs_loc.py: largest cut-cell tilt 73.3 and 82.1 deg, 99th percentile 18.6 and 18.3 deg, 9 and 12
    cut cells above 45 deg, 5 and 6 large-angle pairs, full cells up to 30.6 deg)."""
    import eval_418 as E
    want = {941: (73.3, 18.6, 9, 5), 944: (82.1, 18.3, 12, 6)}
    for s, (mx, p99, n45, npair) in want.items():
        cfg, m, rec = E.load_run(STAIRCASE_RUNS / ("T3Pc8_s%d" % s))
        ids, f = E.geometry(cfg)
        r = E.state_readings(m, ids, f)
        assert (round(r["cut_max"], 1), round(r["cut_p99"], 1), r["cut_gt45"], r["pairs_gt90"]) == (mx, p99, n45, npair)
        assert round(r["full_max"], 1) == 30.6 and rec["relax_converged"]


def test_eval_418_verdict_rules():
    """eval_418.verdict: the rules of step 1 (PLAN 4.18) on synthetic readings."""
    import eval_418 as E
    rd = lambda n45, pairs=0: {"cut_gt45": n45, "pairs_gt90": pairs}
    B8 = {s: rd(10, 5) for s in E.SEEDS}
    R = {"B8": B8, "C8t": {s: rd(0) for s in E.SEEDS}, "C8h": {}, "T4": {s: rd(4) for s in E.SEEDS},
         "T16": {s: rd(1) for s in E.SEEDS}}
    v = E.verdict(R, True)
    assert v["step1"] == "confirmed" and v["k_star"] == 4
    assert E.verdict(R, False)["step1"] == "not determined"                      # 1.0 fails
    R2 = dict(R, C8t={s: rd(6 if s in (941, 942) else 0) for s in E.SEEDS})     # C8t keeps half in 2 seeds
    assert E.verdict(R2, True)["step1"] == "refuted"
    R3 = dict(R, C8t={s: rd(2) for s in E.SEEDS})                                # 1a fails, below half
    assert E.verdict(R3, True)["step1"] == "not determined"
    R4 = dict(R, T4={s: rd(6) for s in E.SEEDS})
    assert E.verdict(R4, True)["k_star"] == 16
    R5 = dict(R, B8={s: rd(2) for s in E.SEEDS})                                 # no seed takes part
    assert E.verdict(R5, True)["step1"] == "not determined"
    R6 = dict(R, T16={941: rd(0), 942: rd(0)})                                   # T16 with 2 values only
    assert E.verdict(R6, True)["step1"] == "not determined"
