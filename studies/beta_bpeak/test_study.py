"""
CPU tests for the beta(B_peak) study helpers.

    CUDA_DEVICE=-1 python -m pytest studies/beta_bpeak/test_study.py -q
"""
import math
import os
import pathlib
import sys

import numpy as np
import pytest

os.environ.setdefault("CUDA_DEVICE", "-1")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from geometry import fcc_box, phi_from, a_from
from anisotropy_model import box_edge, cube_axes, particle_axes


def test_phi_formula():
    d = 1e-6
    for phi in (0.47, 0.65, 0.71):
        assert abs(phi_from(d, a_from(d, phi)) - phi) < 1e-12
    assert abs(phi_from(1.0, math.sqrt(2.0)) - 0.7405) < 1e-4


@pytest.mark.parametrize("N,d", [(72, 0.5e-6), (108, 0.75e-6)])
def test_fcc_box_study_grids(N, d):
    dx = 10.2569e-9
    g = fcc_box(N, d, dx)
    assert abs(g["phi_vox"] - g["phi_nom"]) < 0.02
    assert set(np.unique(g["ids"])) == {-1, 0, 1, 2, 3}
    counts = [np.sum(g["ids"] == p) for p in range(4)]
    assert max(counts) - min(counts) <= 0.01 * max(counts)


def test_fcc_box_rejects_overlap_and_contact():
    from geometry import _check_no_contact
    with pytest.raises(ValueError):
        fcc_box(20, 145e-9, 10e-9)     # d > a/sqrt(2) = 141 nm: overlap
    ids = -np.ones((6, 6, 6), np.int8)
    ids[1, 1, 1], ids[2, 1, 1] = 0, 1  # face neighbours of different particles
    with pytest.raises(ValueError):
        _check_no_contact(ids)
    ids[2, 1, 1] = -1; ids[2, 2, 1] = 1  # edge contact: rejected as well
    with pytest.raises(ValueError):
        _check_no_contact(ids)
    ids[2, 2, 1] = -1; ids[2, 2, 2] = 1  # corner contact: rejected as well
    with pytest.raises(ValueError):
        _check_no_contact(ids)
    ids[2, 2, 2] = -1; ids[3, 1, 1] = 1  # one void cell between: accepted
    _check_no_contact(ids)


def test_drive_value_matches_field():
    import torch
    from magnumnp import Mesh, State
    from field_terms_extra import SinusoidalDrive
    mesh = Mesh((2, 2, 2), (1e-9, 1e-9, 1e-9))
    state = State(mesh)
    drv = SinusoidalDrive(state, (1.0, 1.0, 0.0))
    drv.H_amp, drv.freq, drv.t0 = 1e4, 30e6, 1e-9
    state.t = torch.tensor(9e-9, dtype=torch.float64)
    h = drv.h(state)
    v = drv.value(9e-9)
    assert h.shape == (2, 2, 2, 3)
    assert abs(float(torch.linalg.norm(h[0, 0, 0])) - abs(v)) < 1e-9 * 1e4


def test_units_defaults():
    from units import units, F_REL_DEFAULT
    U = units()
    assert abs(U["l_ex"] - 3.34e-9) < 0.01e-9
    assert abs(U["f_M"] / 1.5 - 28.0e9) < 0.1e9          # 28 GHz/T
    assert abs(F_REL_DEFAULT * U["f_M"] - 30e6) < 1.0
    assert abs(U["Kd"] - 0.5 * U["Js"] ** 2 / 1.2566370614e-6) < 1e-6 * U["Kd"]


@pytest.mark.parametrize("drive", [(1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.3, -0.5, 0.8)])
def test_particle_axes_isotropic(drive):
    """The stress axes are an isotropic set (PROTOKOLL 7.36 / 7.37): orientation tensor <u u^T> = I/3 exactly
    (not only along the drive), and along the drive <cos^2> = 1/3, <cos^4> = 1/5."""
    ax = particle_axes(drive)
    e = np.asarray(drive) / np.linalg.norm(drive)
    assert np.allclose(np.linalg.norm(ax, axis=1), 1.0)
    assert np.allclose(ax.T @ ax / 4.0, np.eye(3) / 3.0, atol=1e-12)
    c = ax @ e
    assert abs(np.mean(c ** 2) - 1 / 3) < 1e-12 and abs(np.mean(c ** 4) - 1 / 5) < 1e-12


@pytest.mark.parametrize("d,phi", [(150.0, 0.65), (300.0, 0.55), (300.0, 0.69)])
def test_box_edge_quantized_and_phi_exact(d, phi):
    a, d_used = box_edge(d, phi)
    assert abs(a / 12.0 - round(a / 12.0)) < 1e-12
    assert abs(phi_from(d_used, a) - phi) < 1e-12
    assert abs(d_used / d - 1.0) < 0.03


def _axes_on_fine_grid(a, d, dx, L, g, seed=5):
    """Cube axes and particle ids of the mesh dx, sampled on a fine grid of step g
    (g divides dx): every fine cell lies completely inside one mesh cell."""
    N = int(round(a / dx))
    geo = fcc_box(N, d, dx)
    ax = cube_axes(geo["ids"], dx, a, d, L, seed)[0]
    k = int(round(dx / g))
    return (np.repeat(np.repeat(np.repeat(geo["ids"], k, 0), k, 1), k, 2),
            np.repeat(np.repeat(np.repeat(ax, k, 0), k, 1), k, 2))


def _n_mismatch(a, d, L, dx1, dx2, g):
    i1, a1 = _axes_on_fine_grid(a, d, dx1, L, g)
    i2, a2 = _axes_on_fine_grid(a, d, dx2, L, g)
    both = (i1 >= 0) & (i1 == i2)
    assert both.sum() > 0.9 * (i1 >= 0).sum()
    return int(np.any(a1[both] != a2[both], axis=-1).sum())


@pytest.mark.parametrize("L", [12.0, 18.0])
def test_cube_axes_same_on_all_meshes(L):
    """The cubes are defined in physical coordinates. On a common fine grid, the
    meshes dx = 3, 2 and 1.5 l_ex give exactly the same axis everywhere."""
    a, d = box_edge(96.0, 0.65)                    # a = 144 l_ex (mesh test geometry)
    assert _n_mismatch(a, d, L, 3.0, 1.5, 1.5) == 0
    assert _n_mismatch(a, d, L, 3.0, 2.0, 1.0) == 0
    assert _n_mismatch(a, d, L, 1.5, 0.75, 0.75) == 0      # convergence series 3 / 1.5 / 0.75 (PLAN 4.10)


def test_cube_axes_mesh_check_detects_bad_L():
    """Negative control: L = 15 l_ex is not a multiple of 6, so the cube faces do
    not lie on the faces of both meshes. The check must see this."""
    a, d = box_edge(96.0, 0.65)
    assert _n_mismatch(a, d, 15.0, 3.0, 2.0, 1.0) > 0


def test_cube_axes_isotropic_and_cube_count():
    a, d = box_edge(150.0, 0.65)
    dx, L = 3.0, 12.0
    N = int(round(a / dx))
    g = fcc_box(N, d, dx)
    ax, n_cubes = cube_axes(g["ids"], dx, a, d, L, seed=1)
    m = g["ids"] >= 0
    assert np.allclose(np.linalg.norm(ax[m], axis=-1), 1.0)
    assert np.all(ax[~m] == 0.0)
    # expected number of cubes that touch a sphere: about (pi/6)(d/L)^3 (plus the surface layer)
    n0 = math.pi / 6.0 * (d / L) ** 3
    assert all(n0 < n < 3.0 * n0 for n in n_cubes)
    # one axis per cube: count distinct axes per particle
    for p in range(4):
        assert len(np.unique(ax[g["ids"] == p].round(12), axis=0)) == n_cubes[p]
    # isotropy of the cube axes (all particles): <(e.x)^2> = 1/3
    u = np.unique(ax[m].round(12), axis=0)
    assert abs(np.mean(u[:, 0] ** 2) - 1.0 / 3.0) < 0.05


def test_cube_axes_depend_on_seed():
    a, d = box_edge(150.0, 0.65)
    N = int(round(a / 3.0))
    g = fcc_box(N, d, 3.0)
    ax1 = cube_axes(g["ids"], 3.0, a, d, 18.0, seed=1)[0]
    ax2 = cube_axes(g["ids"], 3.0, a, d, 18.0, seed=2)[0]
    assert not np.array_equal(ax1, ax2)



def test_volume_fraction_surface():
    """surface = fraction: the packing fraction is exact, fractions lie in [0, 1],
    the interior is 1, and the staircase variant is unchanged (frac 0 or 1)."""
    a, d = box_edge(212.0, 0.65)            # the partial cells need a gap of about 3 cells
    N = int(round(a / 3.0))
    g0 = fcc_box(N, d, 3.0)
    g8 = fcc_box(N, d, 3.0, sub=8)
    assert set(np.unique(g0["frac"])) <= {0.0, 1.0}
    assert np.all((g8["frac"] >= 0.0) & (g8["frac"] <= 1.0))
    assert abs(g8["phi_vox"] - phi_from(d, a)) < 2e-3
    inner = (g0["ids"] >= 0) & (g8["frac"] < 1.0)
    assert inner.sum() < 0.2 * (g0["ids"] >= 0).sum()      # only the surface layer is partial
    assert np.all(g8["ids"][g8["frac"] > 0] >= 0)


# --- domains.py: wall measure and switched volume on synthetic states with known walls --------
# sphere d = 300 l_ex (1 um), L_eff = 12 l_ex, 20 deg ripple; full check: `python domains.py --selftest`

def _dom_state(kind, dx=3.0, ripple=20.0, shift=0.0, seed=0):
    import domains
    st, tr = domains.synthetic_state(kind, dx, d=300.0, ripple_deg=ripple, shift=shift, seed=seed)
    if dx < 3.0:
        m, ic = domains.coarsen_state(st["m"], st["ids"], int(round(3.0 / dx)))
        st = domains.make_state(m, ic, 3.0, st["Leff"], st["name"])
    return st, tr


@pytest.mark.parametrize("kind", ["slab180", "slab90", "bubble"])
def test_texture_gives_wall_area(kind):
    import domains
    st, tr = _dom_state(kind)
    r = domains.measures(st, widths=(domains.W_WALL,))
    assert 0.85 < r["T"][domains.W_WALL] / tr["A"] < 1.0      # selftest: 0.875 ... 0.94


def test_texture_suppresses_ripple_and_is_mesh_independent():
    import domains
    W = (domains.W_WALL,)
    rip = domains.measures(_dom_state("ripple")[0], widths=W)["T"][domains.W_WALL]
    wall, _ = _dom_state("slab180")
    w3 = domains.measures(wall, widths=W)["T"][domains.W_WALL]
    w15 = domains.measures(_dom_state("slab180", dx=1.5)[0], widths=W)["T"][domains.W_WALL]
    assert rip < 0.15 * w3
    assert abs(w15 / w3 - 1.0) < 0.02


@pytest.mark.parametrize("kind", ["slab180", "bubble"])
def test_switched_volume_linear_in_wall_shift(kind):
    import domains
    a, _ = _dom_state(kind)
    b, tr = _dom_state(kind, shift=6.0)
    assert abs(domains.switched_volume(a, b) / tr["V_rev"] - 1.0) < 0.1


def test_texture_no_leak_between_particles():
    """FCC cell of the study (d/l_ex 300, phi 0.65, dx 3): four uniform particles with different
    directions have no wall; the smoothing kernel is wider than the gaps (before the fix: T(2 L_eff)
    = 4.2e4 l_ex^2, about 25 % of a typical wall area)."""
    import domains
    from geometry import fcc_box
    from anisotropy_model import box_edge
    lex = 3.34e-9
    N = int(round(box_edge(300.0, 0.65)[0] / 3.0))
    ids = fcc_box(N, 300.0 * lex, 3.0 * lex)["ids"]
    dirs = np.random.default_rng(1).normal(size=(4, 3))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    m = np.zeros(ids.shape + (3,), np.float32)
    for p in range(4):
        m[ids == p] = dirs[p]
    st = domains.make_state(m, ids, 3.0, 12.0, "uniform particles")
    assert domains.measures(st, widths=(domains.W_WALL,))["T"][domains.W_WALL] < 1.0


def test_wall_path_and_correlation_length():
    """Half-cycle wall path s = V_sw / A_w for 3 walls moved by 6 l_ex (selftest: 5.4, the
    calibration CAL of A_w gives -9 %); l_C grows with the ripple correlation length."""
    import domains
    a, _ = _dom_state("slab180", ripple=13.0)
    b, _ = _dom_state("slab180", ripple=13.0, shift=6.0)
    s = domains.analyze([a, a, b, b])["half_cycle"]["wall_path"][0]
    assert 6.0 * 0.85 < s < 6.0 * 1.15
    st1, _ = domains.synthetic_state("ripple", 3.0, d=300.0, ripple_deg=13.0, ripple_corr=1.0)
    st8, _ = domains.synthetic_state("ripple", 3.0, d=300.0, ripple_deg=13.0, ripple_corr=8.0)
    assert domains.correlation_length(st8) > 3 * domains.correlation_length(st1)


def test_rotation_test_separates_180_walls_at_1um():
    """F90 at d = 300: one 180 deg wall (selftest 12.3 %) is well above a vortex (4.3 %) and a
    helix (0 %)."""
    import domains
    def f90(kind, **kw):
        st, _ = domains.synthetic_state(kind, 3.0, d=300.0, **kw)
        return (domains.rotation_angles(st, window_leff=domains.W_ROT, n_lines=800) > 90.0).mean()
    wall, vortex, helix = f90("single180", ripple_deg=13.0), f90("vortex"), f90("helix")
    assert wall > 2.0 * vortex and helix < 0.01 and wall > 0.08


def test_wall_area_from_rotation_test():
    """A_F90 = 2 F90 V / W for three 180 deg walls (selftest: 1.85e5 vs 1.86e5)."""
    import domains
    st, tr = _dom_state("slab180", ripple=13.0)
    r = domains.measures(st, widths=(domains.W_WALL,))
    assert abs(r["A_F90"] / tr["A"] - 1.0) < 0.15


def test_cube_axes_uniform_on_the_sphere():
    """The cube axes of the study geometry (d/l_ex 300, dx 3, L_eff 12) are uniform on the sphere:
    <n n> = I/3 within 4 SE (SE = sqrt(4/45/n) diagonal, sqrt(1/15/n) off-diagonal), |cos| to each
    axis uniform on [0, 1] (Kolmogorov-Smirnov), neighbour cubes independent (<(n_i.n_j)^2> = 1/3).
    Check over 20 seeds (2026-10-03): <n_x^2>, <n_y^2>, <n_z^2> = 0.3331 / 0.3333 / 0.3335 +- 0.0004."""
    from scipy import stats
    lex = 3.34e-9
    a = box_edge(300.0, 0.65)[0]
    N = int(round(a / 3.0))
    ids = fcc_box(N, 300.0 * lex, 3.0 * lex)["ids"]
    m = ids >= 0
    ax, _ = cube_axes(ids, 3.0, N * 3.0, 300.0, 12.0, seed=901)
    u = np.unique(ax[m].round(12), axis=0)
    n = len(u)
    T = u.T @ u / n
    assert np.all(np.abs(np.diag(T) - 1.0 / 3.0) < 4 * math.sqrt(4.0 / 45.0 / n))
    assert np.abs(T - np.diag(np.diag(T))).max() < 4 * math.sqrt(1.0 / 15.0 / n)
    assert all(stats.kstest(np.abs(u[:, i]), "uniform").pvalue > 1e-3 for i in range(3))
    nb = m & (np.roll(ids, -4, 0) == ids)
    dots = (ax * np.roll(ax, -4, 0)).sum(-1)[nb]
    dots = dots[np.abs(dots) < 1 - 1e-9]
    assert abs(np.mean(dots ** 2) - 1.0 / 3.0) < 0.01


def test_features_power_law_windows():
    """analyze.features: for w = c b^1.7 on the 7 study amplitudes all window betas and beta_pow are 1.7;
    a window with a missing amplitude gives nan; ln w at 10 mT is interpolated."""
    import analyze
    amps = [9, 20, 35, 50, 70, 100, 150]
    rows = [{"b_peak": mT / 1500.0, "w_loop": 2.0 * (mT / 1500.0) ** 1.7, "lnw_err": 0.1} for mT in amps]
    f, e = analyze.features(rows)
    for k in ("beta@9-35mT", "beta@35-70mT", "beta@70-150mT", "beta_pow"):
        assert abs(f[k] - 1.7) < 1e-9 and e[k] > 0
    assert abs(f["lnw@10mT"] - math.log(2.0 * (10 / 1500.0) ** 1.7)) < 1e-9
    f5, _ = analyze.features([r for r in rows if round(r["b_peak"] * 1500) in (9, 20, 50, 100, 150)])
    assert math.isnan(f5["beta@9-35mT"]) and abs(f5["beta_pow"] - 1.7) < 1e-9


def test_sobol_design_mapping():
    """sobol_design.point: fixed mapping (independent of N), ranges and the 4 L_eff levels."""
    import sobol_design as SD
    assert SD.point([0.0, 0.0, 0.0, 0.0]) == dict(Leff_lex=12.0, r_p=0.0, phi=0.55, d_lex=212.0)
    p = SD.point([0.999999, 1.0, 1.0, 1.0])
    assert p["Leff_lex"] == 30.0 and abs(p["r_p"] - 3.0) < 1e-6 and abs(p["phi"] - 0.69) < 1e-6 and abs(p["d_lex"] - 300.0) < 1e-3
    assert [SD.point([u, .5, .5, .5])["Leff_lex"] for u in (0.1, 0.3, 0.6, 0.9)] == [12.0, 18.0, 24.0, 30.0]


def test_domains_batch_selects_last_cycle_quarter_phases(tmp_path):
    """domains_batch.select: 16 phases x 2 cycles -> only the last cycle, phases 90/180/270/360; 4 phases x 1
    cycle (old runs) -> all 4."""
    import domains_batch
    for c in (6, 7):
        for k in range(16):
            (tmp_path / ("m_amp00_c%d_ph%03d.vti" % (c, round(360.0 * (k + 1) / 16)))).touch()
    for deg in (90, 180, 270, 360):
        (tmp_path / ("m_amp03_c4_ph%03d.vti" % deg)).touch()
    sel = domains_batch.select(tmp_path)
    assert [f.name for f in sel["amp00"]] == ["m_amp00_c7_ph%03d.vti" % x for x in (90, 180, 270, 360)]
    assert [f.name for f in sel["amp03"]] == ["m_amp03_c4_ph%03d.vti" % x for x in (90, 180, 270, 360)]


def test_run_queue_expected_bytes():
    """run_queue.expected_bytes: Sobol job = 16 x 2 x 3 snapshots + 2 state files, N from the box edge."""
    import run_queue
    from anisotropy_model import box_edge
    line = open(pathlib.Path(__file__).resolve().parent / "jobs" / "sobol_main.txt").readline()
    name, _, args = line.strip().partition(" ")
    d = run_queue.job_args(args)
    N = round(box_edge(float(d["d_lex"][0]), float(d["phi"][0]))[0] / float(d["dx_lex"][0]))
    assert run_queue.expected_bytes(args) == N ** 3 * (96 * run_queue.B_SNAP_CELL + 2 * run_queue.B_PT_CELL
                                                       + 7 * run_queue.B_ACC_CELL)
    assert run_queue.expected_bytes("--d_lex 300 --phi 0.62 --dx_lex 3") == pytest.approx(
        round(box_edge(300.0, 0.62)[0] / 3) ** 3 * 24.0)


def test_phase_angles_from_index_no_leak():
    """phasemap_batch.phase_angles: equidistant angles from the index (names are rounded: 22, 68, ...). With them
    a static part + harmonics 1 and 3 give exactly these harmonics; with the angles from the rounded names the
    static part leaks into k = 4 (the bug found by the audit 2026-10-06)."""
    import phasemap_batch
    N = 16
    names = [round(360.0 * (j + 1) / N) for j in range(N)]
    assert names[0] == 22 and names[2] == 68
    th = phasemap_batch.phase_angles(names)
    assert np.allclose(th, 2 * np.pi * np.arange(1, N + 1) / N)
    sig = lambda t: 0.9 + 0.01 * np.cos(t) + 0.002 * np.sin(3 * t)
    F = lambda ang, k: np.sum(np.exp(-1j * k * ang) * sig(2 * np.pi * np.arange(1, N + 1) / N)) * 2.0 / N
    assert abs(abs(F(th, 1)) - 0.01) < 1e-12 and abs(abs(F(th, 3)) - 0.002) < 1e-12 and abs(F(th, 4)) < 1e-12
    bad = 2 * np.pi * np.array(names, float) / 360.0
    assert abs(F(bad, 4)) > 0.01                            # the old way: the static 0.9 leaks into k = 4
    with pytest.raises(ValueError):
        phasemap_batch.phase_angles([22, 45, 68])


def test_cleanup_checkpoint_needs_end_state_analysis(tmp_path):
    """cleanup_instance.checkpoint_ok: a snapshot run's checkpoint is deletable only after POSTPROC_DONE, with
    "final" in domains.json and domains.json fetched (local md5); a run without snapshots keeps the old rule."""
    import json
    import cleanup_instance as CI

    def mk(name, snap, postproc=False, final=False):
        d = tmp_path / name
        d.mkdir()
        (d / "DONE").touch()
        (d / "config.json").write_text(json.dumps({"snap_phases": snap}))
        (d / "checkpoint.pt").touch()
        if postproc:
            (d / "POSTPROC_DONE").touch()
            (d / "domains.json").write_text(json.dumps({"stages": [], **({"final": {}} if final else {})}))
        return d
    a = mk("S001", 16)                                   # finished, not yet post-processed
    b = mk("S002", 16, postproc=True, final=True)        # analysed, domains.json not fetched
    c = mk("S003", 16, postproc=True, final=True)        # analysed and fetched
    e = mk("S004", 16, postproc=True, final=False)       # analysed without end state
    f = mk("OLD1", 0)                                    # no snapshots
    (c / "init.pt").touch()
    local = {"%s/%s" % (x, f) for x in ("S001", "S002", "S003", "S004", "OLD1") for f in ("summary.json", "checkpoint.pt")} | \
        {"S003/domains.json", "S004/domains.json", "S003/init.pt"}
    assert [CI.checkpoint_ok(x, local) for x in (a, b, c, e, f)] == [False, False, True, False, True]
    assert not CI.checkpoint_ok(c, local - {"S003/summary.json"})
    assert not CI.checkpoint_ok(c, local - {"S003/checkpoint.pt"})        # checkpoint itself not fetched (C1)
    assert not CI.checkpoint_ok(c, local - {"S003/init.pt"})              # init.pt not fetched (C1)


def _fake_command(args, out):
    """Fake run_loops for the queue tests: the behaviour comes from the job args."""
    code = {"ok": "open(r'%s/DONE','w').write('ok')" % out,
            "fail": "raise SystemExit(1)",
            "gate": "raise SystemExit(3)",
            "gate --relax_maxiter 200000 --no_resume": "open(r'%s/DONE','w').write('ok')" % out,
            "gate --relax_maxiter 200000": "open(r'%s/DONE','w').write('ok')" % out,
            "gate2": "raise SystemExit(3)",
            "gate2 --relax_maxiter 200000 --no_resume": "raise SystemExit(3)",
            "hang": "import time; time.sleep(30)",
            "flaky": ("import os; f=r'%s/tried'\nif os.path.exists(f): open(r'%s/DONE','w').write('ok')\n"
                      "else: open(f,'w').write('1'); raise SystemExit(1)") % (out, out)}[args]
    return [sys.executable, "-c", code]


def _cmd(args, out):
    if args == "boom":
        raise OSError("boom")
    return _fake_command(args, out)


def test_run_queue_retry_gate_timeout_exception(tmp_path):
    """run_queue.run_jobs: a failed job runs once more (flaky -> DONE); exit code 3 is not retried; a hanging
    job is killed by the timeout; an exception in a job does not stop the worker; status file is written."""
    import json
    import run_queue
    jobs = [("A", "ok"), ("B", "flaky"), ("C", "gate"), ("D", "hang"), ("E", "boom"), ("F", "ok")]
    res = run_queue.run_jobs(jobs, ["g0", "g1"], tmp_path / "runs", retries=1, job_timeout_h=2.0 / 3600,
                             status=tmp_path / "status.json", command=_cmd, burst_n=99, status_s=0.2, gate_maxiter=0)
    assert sorted(res["done"]) == ["A", "B", "F"]
    assert set(res["failed"]) == {"C", "D", "E"}
    assert res["failed"]["C"] == "relaxation gate: excluded" and res["failed"]["D"].startswith("timeout")
    assert res["failed"]["E"].startswith("exception")
    st = json.loads((tmp_path / "status.json").read_text())
    assert st["final"] and st["done"] == 3 and not st["running"]
    # gate failure: one retry with raised relax_maxiter and --no_resume; a second gate failure is excluded
    res = run_queue.run_jobs([("G", "gate"), ("H", "gate2")], ["g0"], tmp_path / "runs_g", command=_cmd, burst_n=99)
    assert res["done"] == ["G"] and res["failed"]["H"].endswith("excluded")
    assert run_queue.gate_args("--seed 3 --relax_maxiter 80000", 200000) == "--seed 3 --relax_maxiter 200000 --no_resume"


def test_run_queue_stop_file(tmp_path):
    """PLAN 4.5: when the stop file exists, no new job starts; the jobs left are reported as not run."""
    import run_queue
    stop = tmp_path / "STOP_NEW_JOBS"
    stop.write_text("check_reruns: mismatch\n")
    res = run_queue.run_jobs([("A", "ok"), ("B", "ok")], ["g0"], tmp_path / "runs", command=_cmd, burst_n=99,
                             stop_file=stop)
    assert res["done"] == [] and res["failed"] == {"A": "not run (stop file)", "B": "not run (stop file)"}
    stop.unlink()
    res = run_queue.run_jobs([("A", "ok")], ["g0"], tmp_path / "runs2", command=_cmd, burst_n=99, stop_file=stop)
    assert res["done"] == ["A"]


def test_run_queue_failure_burst_pauses_then_disables_gpu(tmp_path):
    """run_queue.run_jobs: failures of 2 different jobs within the window pause the GPU; the burst_max-th burst
    disables it (ALERT); with one GPU the rest is reported as not run. One job failing twice is no burst."""
    import run_queue
    jobs = [("A", "fail"), ("B", "fail"), ("C", "fail"), ("D", "fail"), ("E", "ok")]
    res = run_queue.run_jobs(jobs, ["g0"], tmp_path / "runs", retries=0, command=_cmd, burst_n=2, burst_s=600,
                             pause_s=0.1, burst_max=2)
    alert = (tmp_path / "runs" / "ALERT").read_text()
    assert "gpu g0 paused" in alert and "gpu g0 disabled" in alert
    assert res["failed"]["E"] == "not run (no GPU left)"
    res = run_queue.run_jobs([("X", "fail"), ("Y", "ok")], ["g0"], tmp_path / "runs2", retries=1, command=_cmd,
                             burst_n=2, burst_s=600, pause_s=0.1, burst_max=1)
    assert res["done"] == ["Y"] and not (tmp_path / "runs2" / "ALERT").exists()


def test_all_study_scripts_compile():
    """Every script of the study compiles (run_loops.py is not imported by the other tests)."""
    import py_compile
    for f in sorted(pathlib.Path(__file__).resolve().parent.glob("*.py")):
        py_compile.compile(str(f), doraise=True)


def test_run_queue_state_survives_restart(tmp_path):
    """run_queue: the retry state is kept in QUEUE_STATE.json (audit 3, C2): an excluded job is not run again
    after a restart; a job in its gate retry resumes with the raised relax_maxiter, without --no_resume."""
    import json
    import run_queue
    runs = tmp_path / "runs"
    res = run_queue.run_jobs([("X", "flaky")], ["g0"], runs, retries=0, command=_cmd, burst_n=99)
    assert res["failed"]["X"] == "exit code 1"
    res = run_queue.run_jobs([("X", "flaky")], ["g0"], runs, retries=0, command=_cmd, burst_n=99)
    assert res["failed"]["X"].startswith("excluded earlier") and not (runs / "X" / "DONE").exists()
    run_queue.save_queue_state(runs / "G", {"attempts": 0, "gate_retried": True, "excluded": None,
                                            "args": "gate --relax_maxiter 200000 --no_resume"})
    res = run_queue.run_jobs([("G", "gate")], ["g0"], runs, command=_cmd, burst_n=99)
    assert res["done"] == ["G"]
    assert json.loads((runs / "X" / "QUEUE_STATE.json").read_text())["excluded"] == "exit code 1"


def _synthetic_recs(shape_of_x, seed=1):
    """Synthetic runs on the real design (results/sobol_design.csv, alpha 0.1): ln E[w](B) = -9 + beta(x) ln(B/50)
    + 0.2 x_phi, beta(x) = 1.8 + 0.15 x_rp (coded x); w = E[w] * Gamma(k(x), 1/k(x)) (mean 1) per amplitude."""
    import csv
    import analyze_design as AD
    rows = [r for r in csv.DictReader(open(pathlib.Path(__file__).resolve().parent / "results" / "sobol_design.csv"))
            if float(r["alpha"]) == 0.1]
    rng = np.random.default_rng(seed)
    recs = []
    for r in rows:
        rec = {"name": r["name"], "lnQ": math.log(1.0 / float(r["Leff_lex"]) ** 2), "r_p": float(r["r_p"]),
               "phi": float(r["phi"]), "lnd": math.log(float(r["d_lex"])), "d_over_L": float(r["d_over_L"]),
               "flagged": False}
        x = AD.coded([rec])[0]
        beta = 1.8 + 0.15 * x[1]
        k = shape_of_x(x)
        lnw = np.array([-9.0 + beta * math.log(m / 50.0) + 0.2 * x[2] for m in AD.AMPS_MT])
        lnw = lnw + np.log(rng.gamma(k, 1.0 / k, len(lnw)))
        for m, v in zip(AD.AMPS_MT, lnw):
            rec["lnw@%gmT" % m] = float(v)
        for m, v in zip(AD.AMPS_MT, lnw):
            rec["lnwdis@%gmT" % m] = float(v)
        rec.update(AD.per_run_outputs(lnw))
        recs.append(rec)
    return recs


def test_analyze_design_recovers_known_effects():
    """analyze_design on synthetic data with a known answer (constant scatter): the true range effects lie in the
    95 % CI of the main fit; r_p is a lever for beta_all (true +0.30), phi for ln E[w] 50 mT (true +0.40); the
    other 6 primary tests are not levers."""
    import json
    import analyze_design as AD
    recs = _synthetic_recs(lambda x: 50.0)
    res = AD.analyse(recs, n_boot=200)
    truth = {("beta_all", f): (0.30 if f == "r_p" else 0.0) for f in AD.FACTORS}
    truth.update({("lnw@50mT", f): (0.40 if f == "phi" else 0.0) for f in AD.FACTORS})
    for k, t in truth.items():
        lo, hi = res["main"][k]["ci95"]
        assert lo - 0.02 <= t <= hi + 0.02, (k, t, lo, hi)
    assert res["classes"][("beta_all", "r_p")] == "lever" and res["classes"][("lnw@50mT", "phi")] == "lever"
    assert sum(v == "lever" for v in res["classes"].values()) == 2
    assert res["boot"]["n_fail"] == 0
    txt, js = AD.report(res), AD.to_json(res)              # the report and the JSON can be written
    assert "Primary outputs" in txt and "lever" in txt and json.loads(js)["n_runs"] == len(recs)


def test_analyze_design_mean_estimand_under_heteroscedastic_scatter():
    """The scatter grows with r_p (Gamma shape 64 -> 4), the mean E[w] does not depend on r_p. The main fit (mean)
    finds no r_p effect on ln E[w]; fit A (OLS of ln w, E[ln w]) is biased by -(1/8 - 1/128) ~ -0.12; the scatter
    model finds the r_p slope."""
    import analyze_design as AD
    recs = _synthetic_recs(lambda x: 16.0 * 4.0 ** (-x[1]), seed=2)
    res = AD.analyse(recs, n_boot=200)
    m, a = res["main"][("lnw@50mT", "r_p")], res["A"][("lnw@50mT", "r_p")]
    assert m["ci95"][0] - 0.02 <= 0.0 <= m["ci95"][1] + 0.02
    assert a["eff"] < -0.05 and abs(m["eff"]) < abs(a["eff"])
    assert res["scatter"]["lnw@50mT"]["significant"]["r_p"]


def test_steady_tail_rule():
    """Steady-tail rule (PROTOKOLL 7.54): only the cycles after the last event cycle (|w_dis/w_loop - 1| > tol or
    w_loop <= 0) are kept; an event in the last cycle leaves nothing. A missing amplitude gives nan at its nominal
    amplitude only; the slopes use the other amplitudes; the main fit and fit A run with the gaps and recover the
    known effects."""
    import analyze
    import analyze_design as AD
    c = lambda wl, wd: {"w_loop": wl, "w_dis": wd}
    cyc = [c(1, 1.5), c(1, 1.1), c(1, 2.0), c(1, 1.0), c(1, 0.9)]
    assert analyze.steady_tail(cyc, 0.2) == cyc[3:]
    assert analyze.steady_tail(cyc[:3], 0.2) == []
    assert analyze.steady_tail([c(-1, 1), c(1, 1)], 0.2) == [c(1, 1)]
    assert analyze.steady_tail(cyc[3:], 0.2) == cyc[3:]
    bm = np.array(AD.AMPS_MT) * 1.02
    lw = 2.0 * np.log(bm) - 9.0
    full = AD.to_nominal(bm, lw, np.zeros(7, bool))
    assert np.allclose(full, 2.0 * np.log(AD.AMPS_MT) - 9.0)
    miss = np.zeros(7, bool)
    miss[2] = True
    gap = AD.to_nominal(bm, lw, miss)
    assert np.isnan(gap[2]) and np.allclose(np.delete(gap, 2), np.delete(full, 2))
    o = AD.per_run_outputs(gap)
    assert abs(o["beta_all"] - 2.0) < 1e-9 and abs(o["beta@9-35mT"] - 2.0) < 1e-9
    assert math.isfinite(o["lnw@50mT"]) and math.isfinite(o["lnw@10mT"])
    o = AD.per_run_outputs(np.where(np.arange(7) == 1, np.nan, full))
    assert math.isnan(o["lnw@10mT"])                        # 10 mT needs 9 and 20 mT
    recs = _synthetic_recs(lambda x: 50.0, seed=3)
    rng = np.random.default_rng(4)
    for r in rng.choice(recs, 6, replace=False):              # 6 gaps at 20 / 35 mT, as in the real data
        k = int(rng.choice([1, 2]))
        r["lnw@%gmT" % AD.AMPS_MT[k]] = float("nan")
        lnw = np.array([r["lnw@%gmT" % m] for m in AD.AMPS_MT])
        r.update(AD.per_run_outputs(lnw))
    res = AD.analyse(recs, n_boot=100)
    lo, hi = res["main"][("beta_all", "r_p")]["ci95"]
    assert lo - 0.02 <= 0.30 <= hi + 0.02
    assert res["boot"]["n_fail"] == 0


def test_acc_reduce_relative_slice_is_finite():
    """PROTOKOLL 7.61: the loss-density slice is stored relative to its mean (float16 overflowed at 65504 W/m³)."""
    from acc_reduce import relative_slice
    plane = np.full((4, 4), 3e8)
    plane[1, 1] = 9e9
    inside = np.ones((4, 4), bool)
    inside[0, 0] = False
    plane[0, 0] = 0.0
    rel, mean, scale = relative_slice(plane, inside)
    assert rel.dtype == np.float16 and np.all(np.isfinite(rel))
    assert abs(mean - plane[inside].mean()) < 1e-6 * mean
    assert abs(scale - mean) < 1e-6 * mean                     # a quantity >= 0: scale = mean, rel = p / mean as before
    assert np.allclose(rel.astype(float) * mean, plane, rtol=2e-3)


def test_acc_reduce_signed_slice_round_trip(tmp_path):
    """PROTOKOLL 7.90: a signed plane (the rest) with a mean <= 0 or near 0 is stored, not zeroed; slice_values reads
    new files (rel x scale) and old files (rel x mean, no scale); plane_known marks an old zeroed plane."""
    from acc_reduce import relative_slice, slice_values, plane_known
    rng = np.random.default_rng(3)
    inside = np.ones((20, 20), bool)
    inside[:3] = False
    for shift in (-0.3, -1e-9, 0.0, 1e-9, 0.3):              # mean below 0, near 0 and above 0
        plane = rng.normal(shift, 1.0, (20, 20)) * 1e7
        plane[~inside] = 0.0
        plane[inside] -= plane[inside].mean() - shift * 1e7
        rel, mean, scale = relative_slice(plane, inside)
        assert scale > 0 and np.all(np.isfinite(rel.astype(float)))
        np.savez(tmp_path / "z.npz", a_rest_rel=rel, a_rest_mean=mean, a_rest_scale=scale)
        z = np.load(tmp_path / "z.npz")
        assert plane_known(z, "a", "rest")
        assert np.allclose(slice_values(z, "a", "rest"), plane, rtol=0, atol=2e-3 * np.abs(plane).max())
    one = np.zeros((4, 4))
    one[2, 2] = 1.0                                           # one cell holds everything: rel = 16 fits float16
    big = np.zeros((300, 300))
    big[5, 5] = 1.0                                           # rel = 9e4 > float16 maximum: stored as float32
    for pl, dt in ((one, np.float16), (big, np.float32)):
        rel, mean, scale = relative_slice(pl, np.ones(pl.shape, bool))
        assert rel.dtype == dt and np.isclose(float(rel.max()) * scale, 1.0, rtol=1e-3)
    np.savez(tmp_path / "old.npz", a_rest_rel=np.zeros((2, 2), np.float16), a_rest_mean=-2.0,
             b_rest_rel=np.full((2, 2), 1.0, np.float16), b_rest_mean=5.0)
    z = np.load(tmp_path / "old.npz")
    assert not plane_known(z, "a", "rest") and plane_known(z, "b", "rest")
    assert np.allclose(slice_values(z, "b", "rest"), 5.0)


def test_min_calm_stop_rule():
    """PLAN 4.5: the calm tail counts the calm cycles at the end of the measured cycles (from 2 on after a drive
    correction); it equals the number of cycles that analyze.steady_tail keeps. With --min_calm N > 0 a stage stops
    when its calm tail is >= N or at the maximum (also when the stage has more cycles than max_cycles_per_amp);
    with N = 0 the old closure / dW rule holds."""
    import types
    import analyze
    from run_loops import calm_tail, measured_cycles, stop_stage
    rng = np.random.default_rng(5)
    for _ in range(200):
        cyc = [{"w_loop": float(rng.choice([-1.0, 1.0, 1.0, 1.0])), "w_dis": float(rng.choice([0.9, 1.0, 1.5]))}
               for _ in range(int(rng.integers(0, 9)))]
        assert calm_tail(cyc, 0.2) == len(analyze.steady_tail(cyc, 0.2))
    c = lambda wd, cl=0.001, dw=0.01: {"w_loop": 1.0, "w_dis": wd, "closure": cl, "dW_rel": dw}
    args = types.SimpleNamespace(min_calm=2, calm_tol=0.2, closure_tol=0.02, dW_tol=0.1, max_cycles_per_amp=15)
    rec = {"h_amp_initial": 1.0, "cycles": [c(5.0), c(5.0), c(1.0), c(1.0), c(1.0), c(1.5), c(1.0)]}
    assert measured_cycles(rec) == rec["cycles"][2:]
    stop, info = stop_stage(rec, 7, 7, args)
    assert not stop and info["calm_tail"] == 1                   # event in cycle 5: one calm cycle
    rec["cycles"].append(c(1.05))
    stop, info = stop_stage(rec, 8, 7, args)
    assert stop and info["calm_ok"] and info["calm_tail"] == 2
    rec = {"cycles": [c(1.0)] + [c(2.0)] * 14}                    # never calm: stops at the maximum
    assert not stop_stage(rec, 14, 7, args)[0] and stop_stage(rec, 15, 7, args)[0]
    assert not stop_stage(rec, 15, 30, args)[0] and stop_stage(rec, 30, 30, args)[0]   # --cycles_last_amp 30
    old = types.SimpleNamespace(min_calm=0, calm_tol=0.2, closure_tol=0.02, dW_tol=0.1, max_cycles_per_amp=7)
    rec = {"cycles": [c(2.0, 0.001, 0.01)] * 7}                   # old rule: closed and in scatter -> stop
    assert stop_stage(rec, 7, 7, old)[0]
    rec = {"cycles": [c(1.0, 0.05, 0.5)] * 5}
    assert not stop_stage(rec, 5, 5, old)[0] and stop_stage(rec, 7, 5, old)[0]


def test_rerun_jobs_and_gate(tmp_path):
    """PLAN 4.5: rerun lines keep the original arguments, set --init_from, --min_calm 2 and the maximum, and drop
    the snapshots; the conditioning lines take the first 1 / 2 / 3 amplitudes. The gate compares cycle by cycle up
    to the first stage where the numbers of cycles differ."""
    import shlex
    import check_reruns as CR
    import rerun_jobs as RJ
    line = ("--d_lex 252.0 --cycles_per_amp 7 --max_cycles_per_amp 7 --b_list 0.006 0.013333 0.023333 0.1 --seed 9 "
            "--acc --snap_phases 16 --snap_cycles 2 --snap_mT 9.0 50.0 150.0")
    t = RJ.strip_snap(shlex.split(line))
    assert "--snap_mT" not in t and "--acc" in t and RJ.get_arg(t, "b_list") == ["0.006", "0.013333", "0.023333", "0.1"]
    t = RJ.set_arg(t, "max_cycles_per_amp", ["15"])
    assert RJ.get_arg(t, "max_cycles_per_amp") == ["15"] and t.count("--max_cycles_per_amp") == 1
    src = tmp_path / "jobs.txt"
    src.write_text("C0 " + line + "\n")
    cl = RJ.cond(["C0"], str(src))
    assert [x.split()[0] for x in cl] == ["cond_C0_b9", "cond_C0_b20", "cond_C0_b35"]
    t = shlex.split(cl[1].partition(" ")[2])
    assert RJ.get_arg(t, "b_list") == ["0.006", "0.013333"] and RJ.get_arg(t, "cycles_last_amp") == ["30"]
    assert RJ.get_arg(t, "init_from") is None and RJ.get_arg(t, "min_calm") == ["2"]   # same start path (7.59)
    runs = tmp_path / "runs"
    for n, cfg in (("S1_mc2", '{"init_from": "S1/init.pt"}'), ("S2_mc2", '{"init_from": null}'),
                   ("M1_mc2", '{"init_from": "M1_b9/init.pt"}')):
        (runs / n).mkdir(parents=True)
        (runs / n / "config.json").write_text(cfg)
    qf = tmp_path / "queue.txt"
    qf.write_text("S1_mc2 --seed 1 --min_calm 2\nS2_mc2 --seed 2\nM1_mc2 --init_from M1_b9/init.pt\n")
    assert RJ.retire(qf, runs) == ["S1_mc2 -> S1_mc2init"]
    assert (runs / "S1_mc2init").exists() and (runs / "S2_mc2").exists() and (runs / "M1_mc2").exists()
    assert CR.reference("S033_mc2") == "S033" and CR.reference("cond_C3_b20") == "C3" and CR.reference("S033") is None
    cy = lambda v: {"w_loop": v, "w_dis": v, "b_peak": 0.01}
    ref = [{"name": "amp00", "cycles": [cy(1.0)] * 7}, {"name": "amp01", "cycles": [cy(2.0)] * 7},
           {"name": "amp02", "cycles": [cy(3.0)] * 7}]
    new = [{"name": "amp00", "cycles": [cy(1.0)] * 7}, {"name": "amp01", "cycles": [cy(2.0)] * 7 + [cy(2.5)] * 3},
           {"name": "amp02", "cycles": [cy(9.0)] * 7}]
    assert CR.compare(ref, new) == (14, None)                    # amp02 starts from another state: not compared
    new[0]["cycles"] = [cy(1.0)] * 6 + [cy(1.0 + 1e-4)]
    n, why = CR.compare(ref, new)
    assert n == 6 and why is not None


def test_cond_jobs_alpha001_keep_design_history():
    """Decision Chris 2026-10-12: α 0.01 long runs at 9 and 20 mT; the earlier stage keeps the design length (9 cycles),
    the start state stays C<k>/init.pt, no snapshots; the gate maps cond_C0_a001_b20 to C0_a001."""
    import shlex
    import check_reruns as CR
    import rerun_jobs as RJ
    src = pathlib.Path(__file__).resolve().parent / "jobs" / "sobol_a001.txt"
    lines = RJ.cond(["C0_a001"], str(src), (9, 20), 9)
    assert [x.split()[0] for x in lines] == ["cond_C0_a001_b9", "cond_C0_a001_b20"]
    t = shlex.split(lines[1].partition(" ")[2])
    assert RJ.get_arg(t, "b_list") == ["0.006", "0.013333"] and RJ.get_arg(t, "max_cycles_per_amp") == ["9"]
    assert RJ.get_arg(t, "cycles_last_amp") == ["30"] and RJ.get_arg(t, "init_from") == ["C0/init.pt"]
    assert RJ.get_arg(t, "alpha") == ["0.01"] and "--snap_phases" not in t
    assert CR.reference("cond_C0_a001_b20") == "C0_a001"


def test_pair_and_freq_jobs():
    """Decisions Chris 2026-10-12: K1A10 = K1 with d/l_ex 300 (A 10 pJ/m), snapshots kept, steady rule; frequency test
    C′ = the pair at f/2 and 2 f, no snapshots, the rest equal."""
    import shlex
    import os
    import rerun_jobs as RJ
    cwd = os.getcwd()
    os.chdir(pathlib.Path(__file__).resolve().parent)
    try:
        import pair_jobs as PJ
        pl, fl = PJ.pair_lines(), PJ.freq_lines()
        src = dict(RJ.job_lines("jobs/sobol_corners.txt"))
    finally:
        os.chdir(cwd)
    assert [n for n, _ in pl] == ["K1A10_0", "K1A10_1", "K1A10_2"]
    for n, t in pl:
        r = shlex.split(src["K1_" + n[-1]])
        assert RJ.get_arg(t, "d_lex") == ["300.0"] and RJ.get_arg(t, "seed") == RJ.get_arg(r, "seed")
        assert RJ.get_arg(t, "snap_phases") == ["16"] and RJ.get_arg(t, "min_calm") == ["2"]
    for n, t in pl:
        assert RJ.get_arg(t, "acc_harm") == ["7"] and RJ.get_arg(t, "samples") == ["512"]   # O2 (7.73), T / 512 (7.76)
        assert 512 % int(RJ.get_arg(t, "snap_phases")[0]) == 0
    # C′: the pair at f/2 and 2 f, K1 also at f (7.73); harmonics; the same sample interval T(f) / 512 (7.75)
    assert len(fl) == 15 and fl[0][0] == "K1A10_0_fhalf"
    assert sorted(n for n, _ in fl if n.endswith("_f1")) == ["K1_0_f1", "K1_1_f1", "K1_2_f1"]
    for n, t in fl:
        fac = 0.5 if n.endswith("fhalf") else (1.0 if n.endswith("_f1") else 2.0)
        f = float(RJ.get_arg(t, "f_rel")[0]) / PJ.F_REL_30MHZ if "--f_rel" in t else 1.0
        assert abs(f - fac) < 1e-6 and "--snap_phases" not in t
        assert RJ.get_arg(t, "d_lex") == (["300.0"] if n.startswith("K1A10") else ["212.0"])
        assert RJ.get_arg(t, "acc_harm") == ["7"] and RJ.get_arg(t, "samples") == [str(int(512 / fac))]


def test_size_jobs_equal_mesh300_except_size_and_seed():
    """PLAN 4.6: each line of the d/l_ex 212 mesh check equals its mesh300 template of seed 921 except d_lex, seed, init_from and the
    steady-rule extension; 50 / 100 mT start from the 9 mT state of the same size, mesh and seed."""
    import shlex
    import rerun_jobs as RJ
    import size_jobs as SJ
    src = pathlib.Path(__file__).resolve().parent / "jobs" / "mesh300.txt"
    tmpl = dict(RJ.job_lines(src))
    out = SJ.lines(src)
    assert len(out) == 24 and len({n for n, _ in out}) == 24
    assert [n for n, _ in out[:4]] == ["M212_dx15_s%d_b9" % s for s in range(931, 935)]
    for name, a in out:
        size, mesh, s, b = name.split("_")
        t, r = shlex.split(a), shlex.split(tmpl["M300_%s_s921_%s" % (mesh, b)])
        assert RJ.get_arg(t, "d_lex") == ["%s.0" % size[1:]] and RJ.get_arg(t, "seed") == [s[1:]]
        assert RJ.get_arg(t, "min_calm") == ["2"] and RJ.get_arg(t, "max_cycles_per_amp") == ["12"]
        init = RJ.get_arg(t, "init_from")
        assert init == (None if b == "b9" else ["%s_%s_%s_b9/init.pt" % (size, mesh, s)])
        keys = {x[2:] for x in r if x.startswith("--")}
        for k in keys - set(SJ.CHANGED):
            assert RJ.get_arg(t, k) == RJ.get_arg(r, k), (name, k)
        assert {x[2:] for x in t if x.startswith("--")} <= keys | set(SJ.CHANGED)


def test_eval_size_rules():
    """PLAN 4.6: no statement from fewer than 3 seeds; Welch difference of two summaries; sign change only if both
    90 % CIs exclude 0; the mesh error of one seed from the stage dicts."""
    import eval_size as ES
    assert ES.summary([0.1, 0.2]) is None
    s = ES.summary([0.1, 0.2, 0.3])
    assert abs(s[0] - 0.2) < 1e-12 and s[2] == 3 and abs(s[3] - 2.919986 * 0.1 / math.sqrt(3)) < 1e-5
    a, b = (0.5, 0.1, 4, 0.2), (0.2, 0.1, 4, 0.2)
    d = ES.difference(a, b)
    assert abs(d[0] - 0.3) < 1e-12 and abs(d[1] - math.sqrt(0.02)) < 1e-12
    assert abs(d[2] - stats_t(0.95, 6.0) * math.sqrt(0.02)) < 1e-9          # equal SE and n: dof = 2 (n - 1)
    assert ES.sign_change((0.5, 0.1, 4, 0.2), (-0.5, 0.1, 4, 0.2))
    assert not ES.sign_change((0.5, 0.1, 4, 0.2), (-0.1, 0.1, 4, 0.2)) and ES.difference(a, None) is None
    fine = {9: (1.0, 9.0), 50: (40.0, 50.0), 100: (200.0, 100.0)}
    coarse = {9: (2.0, 9.0), 50: (40.0, 50.0), 100: (200.0, 100.0)}
    e = ES.mesh_error(fine, coarse)
    assert abs(e["lnP 9"] - math.log(2.0)) < 1e-12 and e["lnP 50"] == 0.0
    assert abs(e["beta 9-50"] + math.log(2.0) / math.log(50.0 / 9.0)) < 1e-12 and e["beta"] < 0
    assert "lnP 100" not in ES.mesh_error(fine, {9: (2.0, 9.0)}) and "beta" not in ES.mesh_error(fine, {9: (2.0, 9.0)})


def test_eval_cond_decision_rule():
    """PLAN 4.5 D: design value stands if the 90 % CI of the mean drift lies inside ±0.05; systematic error if the
    95 % CI excludes 0 and |mean| >= 0.05; else not determined; fewer than 3 seeds: not determined."""
    import eval_cond as EC
    assert EC.decide([0.01, -0.01, 0.0, 0.005, -0.005, 0.0])[-1] == "design value stands"
    assert EC.decide([0.10, 0.11, 0.09, 0.10, 0.12, 0.08])[-1] == "systematic error"
    assert EC.decide([0.10, -0.08, 0.15, -0.02, 0.0, 0.05])[-1] == "not determined"
    assert EC.decide([0.0, 0.0])[-1].startswith("not determined")
    cy = [{"cycle": i, "w_loop": 1.0, "w_dis": 1.0, "b_peak": 0.006} for i in range(30)]
    cy[3]["w_dis"] = 2.0
    cy[28]["w_loop"] = cy[28]["w_dis"] = 1.2
    design, late, B, ev = EC.values({"h_amp_initial": 0.1, "cycles": cy})
    assert design == 1.0 and abs(late - 1.1) < 1e-12 and ev == [3] and abs(B - 9.0) < 1e-9


def stats_t(q, dof):
    from scipy import stats
    return stats.t.ppf(q, dof)


def test_analyze_design_rerun_mapping(tmp_path):
    """PLAN 4.5: with reruns=True a finished X_mc2 replaces X (name = folder used, orig = X); without, the original
    is used; an unfinished rerun is not used."""
    import shutil
    import analyze_design as AD
    src = pathlib.Path(__file__).resolve().parent / "runs" / "S005"
    if not (src / "summary.json").exists():
        pytest.skip("run S005 not present")
    for name in ("S005", "S005_mc2", "S006", "S006_mc2"):
        (tmp_path / name).mkdir()
        for f in ("config.json", "summary.json"):
            shutil.copy(src / f, tmp_path / name / f)
    for name in ("S005", "S005_mc2", "S006"):
        (tmp_path / name / "DONE").write_text("ok\n")
    recs, _ = AD.load_runs(tmp_path, reruns=True)
    assert sorted((r["name"], r["orig"]) for r in recs) == [("S005_mc2", "S005"), ("S006", "S006")]
    recs, _ = AD.load_runs(tmp_path)
    assert sorted(r["name"] for r in recs) == ["S005", "S006"]


def test_analyze_design_flag_rule_on_archived_runs():
    """Flag rule (PROTOKOLL 7.39): a run is flagged if |w_dis/w_loop - 1| > 0.20 at any amplitude. On the archived
    runs of the stopped design (if present) this flags a minority, not all runs as the old rule did."""
    import analyze_design as AD
    arch = pathlib.Path(__file__).resolve().parent / "runs" / "_archive_sobol_v2_oldaxes"
    if not arch.exists():
        pytest.skip("archive not present")
    recs, _ = AD.load_runs(arch)
    n_flag = sum(r["flagged"] for r in recs)
    assert len(recs) == 20 and 0 < n_flag < len(recs) // 2
    assert all(all(math.isfinite(r["lnwdis@%gmT" % m]) for m in AD.AMPS_MT) for r in recs)


def test_virgin_m0_isotropic_and_seeded():
    """Physics invariant (PROTOKOLL 7.40): the random initial state has no preferred direction: the mean over the
    blocks is ~0 within its statistical error and <m m^T> = I/3 within the error; it depends on the seed; it is
    constant on the blocks (init_block 6 l_ex at dx 3 = 2 cells)."""
    from run_loops import virgin_m0
    N, a, dx, blk = 96, 288.0, 3.0, 6.0
    m = virgin_m0(N, a, dx, blk, 1001)
    nb = int(round(a / blk))
    blocks = m[::2, ::2, ::2].reshape(-1, 3)                 # one value per block
    assert blocks.shape[0] == nb ** 3
    assert np.allclose(m[0, 0, 0], m[1, 1, 1])               # constant on a block
    se = 1.0 / math.sqrt(3 * nb ** 3)
    assert np.all(np.abs(blocks.mean(0)) < 5 * se)
    assert np.allclose(blocks.T @ blocks / len(blocks), np.eye(3) / 3, atol=5 * 2 / math.sqrt(45 * nb ** 3) + 1e-3)
    assert not np.allclose(m, virgin_m0(N, a, dx, blk, 1002))


def test_drive_zero_mean_and_phase():
    """Physics invariant (PROTOKOLL 7.40): the drive is a pure sine starting at phase 0: zero mean over the 256
    samples of a cycle, H(t0 + T/4) = +H_amp, H(t0 + 3T/4) = -H_amp."""
    from magnumnp import Mesh, State
    from field_terms_extra import SinusoidalDrive
    drv = SinusoidalDrive(State(Mesh((2, 2, 2), (1e-9, 1e-9, 1e-9))), (1.0, 0.0, 0.0))
    drv.H_amp, drv.freq, drv.t0 = 1.0e4, 30e6, 2e-9
    T = 1.0 / drv.freq
    v = np.array([drv.value(drv.t0 + (k + 1) * T / 256) for k in range(256)])
    assert abs(v.mean()) < 1e-9 * drv.H_amp
    assert abs(drv.value(drv.t0 + T / 4) - drv.H_amp) < 1e-9 * drv.H_amp
    assert abs(drv.value(drv.t0 + 3 * T / 4) + drv.H_amp) < 1e-9 * drv.H_amp


def test_dissipation_formulas_agree():
    """Physics invariant (PROTOKOLL 7.40): the dissipation of run_loops, p = alpha gamma mu0 Ms |m x H|^2 / (1 + alpha^2)
    (magnum.np gamma = 2.21e5 m/(A s), includes mu0), equals alpha mu0 Ms / gamma |dm/dt|^2 of the LLG equation, which
    phasemap_batch / phase_extra / acc_reduce use (GAMMA = 2.21276157e5), and equals mu0 Ms H . dm/dt (= -dE/dt)."""
    from magnumnp import constants
    import acc_reduce
    rng = np.random.default_rng(3)
    alpha, Ms = 0.1, 1.5 / constants.mu_0
    m = rng.standard_normal((50, 3)); m /= np.linalg.norm(m, axis=1, keepdims=True)
    H = 1e4 * rng.standard_normal((50, 3))
    g = constants.gamma
    mxH = np.cross(m, H)
    dmdt = -g / (1 + alpha ** 2) * (mxH + alpha * np.cross(m, mxH))
    p_run = alpha * g * constants.mu_0 / (1 + alpha ** 2) * Ms * (mxH ** 2).sum(1)
    p_llg = alpha * acc_reduce.MU0 * Ms / acc_reduce.GAMMA * (dmdt ** 2).sum(1)
    p_work = constants.mu_0 * Ms * (H * dmdt).sum(1)          # -dE/dt with dE/dt = -mu0 Ms H . dm/dt
    assert abs(acc_reduce.GAMMA / g - 1) < 1e-9
    assert np.allclose(p_run, p_llg, rtol=1e-9) and np.allclose(p_run, p_work, rtol=1e-9)


def test_demag_pbc_slab_and_zero_mean():
    """Physics invariant (PROTOKOLL 7.40): true periodic demag (k = 0 removed). A periodic stack of slabs (half of
    the box along z, m along z) has H = -Ms/2 inside and +Ms/2 outside (zero mean: an infinite medium has no
    sample-shape demag); a uniformly magnetised full box has H = 0; for any m the box mean of H is 0."""
    import torch
    from magnumnp import Mesh, State
    from magnumnp.field_terms.demagPBC import DemagFieldPBC
    n = (8, 8, 8)
    mesh = Mesh(n, (1e-9, 1e-9, 1e-9))
    state = State(mesh)
    Ms = 1.0e6
    ms = torch.zeros(n + (1,), dtype=torch.float64)
    ms[:, :, :4] = Ms
    state.material = {"Ms": ms}
    m = torch.zeros(n + (3,), dtype=torch.float64)
    m[..., 2] = 1.0
    state.m = m
    h = DemagFieldPBC().h(state)
    assert torch.allclose(h[:, :, :4, 2], torch.full_like(h[:, :, :4, 2], -Ms / 2), rtol=1e-9, atol=1e-6)
    assert torch.allclose(h[:, :, 4:, 2], torch.full_like(h[:, :, 4:, 2], Ms / 2), rtol=1e-9, atol=1e-6)
    state.material = {"Ms": torch.full(n + (1,), Ms, dtype=torch.float64)}
    assert float(DemagFieldPBC().h(state).abs().max()) < 1e-6 * Ms
    g = torch.Generator().manual_seed(0)
    r = torch.randn(n + (3,), generator=g, dtype=torch.float64)
    state.m = r / torch.linalg.norm(r, dim=-1, keepdim=True)
    assert float(DemagFieldPBC().h(state).mean(dim=(0, 1, 2)).abs().max()) < 1e-9 * Ms


def test_stress_term_uses_its_own_parameters():
    """Physics invariant (PROTOKOLL 7.40): the second uniaxial term of run_loops, UniaxialAnisotropyField(Ku="Kp",
    Ku_axis="Kp_axis"), acts with K_p along the particle axis and does not fall back to Ku / Ku_axis (else r_p would
    be a copy of the Herzer term); h = 2 K (m.u) u / (mu0 Ms); easy axis (energy lower for m || u than m _|_ u)."""
    import torch
    from magnumnp import Mesh, State, constants
    from magnumnp.field_terms import UniaxialAnisotropyField
    n = (2, 2, 2)
    state = State(Mesh(n, (1e-9, 1e-9, 1e-9)))
    Ms, K = 1.0e6, 5.0e3
    u = torch.tensor([0.6, 0.0, 0.8], dtype=torch.float64)
    v = torch.tensor([0.0, 1.0, 0.0], dtype=torch.float64)            # the Herzer axis, different from u
    state.material = {"Ms": torch.full(n + (1,), Ms, dtype=torch.float64),
                      "Ku": torch.zeros(n + (1,), dtype=torch.float64), "Ku_axis": v.expand(n + (3,)).clone(),
                      "Kp": torch.full(n + (1,), K, dtype=torch.float64), "Kp_axis": u.expand(n + (3,)).clone()}
    m = torch.tensor([1.0, 0.0, 0.0], dtype=torch.float64)
    state.m = m.expand(n + (3,)).clone()
    herzer, stress = UniaxialAnisotropyField(), UniaxialAnisotropyField(Ku="Kp", Ku_axis="Kp_axis")
    assert float(herzer.h(state).abs().max()) == 0.0                  # Ku = 0: the default term is silent
    h = stress.h(state)[0, 0, 0]
    expect = 2 * K / (constants.mu_0 * Ms) * float(m @ u) * u
    assert torch.allclose(h, expect, rtol=1e-12)
    state.m = u.expand(n + (3,)).clone()
    E_par = float(stress.E(state))
    state.m = v.expand(n + (3,)).clone()
    E_perp = float(stress.E(state))
    assert E_par < E_perp                                               # easy axis


def _harm_signal(theta, amps):
    """m(theta) of shape (2, 2, 2, 3): a constant plus harmonics (k, component, amplitude, phase); the same in all
    cells. Returns m and dm/dtheta."""
    m = np.zeros((2, 2, 2, 3)) + np.array([0.1, -0.2, 0.3])
    dm = np.zeros_like(m)
    for k, comp, a, ph in amps:
        m[..., comp] += a * np.cos(k * theta + ph)
        dm[..., comp] += -a * k * np.sin(k * theta + ph)
    return m, dm


def test_harmacc_recovers_harmonics_and_skips_event_cycles():
    """run_loops.HarmAcc (--acc_harm, PROTOKOLL 7.72): |m_k|^2 of known harmonics, no leak into other k, an event
    cycle does not enter the stage sums, the samples after each step (theta = 2 pi (j + 1) / n) as in run_loops."""
    torch = pytest.importorskip("torch")
    import run_loops
    K, n = 7, 256
    amps = [(1, 0, 0.30, 0.2), (3, 1, 0.05, -1.0), (5, 2, 0.02, 0.5), (7, 0, 0.01, 2.0)]
    h = run_loops.HarmAcc((2, 2, 2, 3), K, "cpu", torch.ones((2, 2, 2), dtype=torch.bool))
    h.new_stage()
    for cyc, calm in enumerate((True, False, True)):
        for j in range(n):
            th = 2 * math.pi * (cyc + (j + 1) / n)
            amp = amps if calm else [(1, 0, 0.9, 0.0)]           # the event cycle has another motion
            m, _ = _harm_signal(th, amp)
            h.add(torch.as_tensor(m), torch.ones((2, 2, 2), dtype=torch.float64), th, 1.0, 1.0)
        box = h.end_cycle(calm, cyc)
        assert box["calm"] == calm and box["n"] == n
    assert h.kept == [0, 2] and h.n_cycles == 3 and h.s_n == 2 * n
    mk2 = (h.s_mk2 / len(h.kept)).numpy()[:, 0, 0, 0]
    want = np.zeros(K)
    for k, comp, a, ph in amps:
        want[k - 1] += a * a
    assert np.allclose(mk2, want, atol=1e-6)
    assert float(h.s["p"][0, 0, 0]) == 2 * n


def test_harmonic_loss_parseval():
    """acc_reduce.harmonic_loss on |m_k|^2 from HarmAcc equals the mean Gilbert loss alpha mu0 Ms / gamma |dm/dt|^2
    of a motion with harmonics up to K only (Parseval); the rest is 0."""
    torch = pytest.importorskip("torch")
    import run_loops
    import acc_reduce
    K, n, f, alpha, Ms = 7, 256, 3.0e7, 0.1, 1.5 / acc_reduce.MU0
    w = 2 * math.pi * f
    amps = [(1, 0, 0.30, 0.2), (2, 1, 0.04, 0.3), (3, 1, 0.05, -1.0), (7, 2, 0.01, 2.0)]
    h = run_loops.HarmAcc((2, 2, 2, 3), K, "cpu", torch.ones((2, 2, 2), dtype=torch.bool))
    h.new_stage()
    p_mean = 0.0
    for j in range(n):
        th = 2 * math.pi * (j + 1) / n
        m, dm = _harm_signal(th, amps)
        p = alpha * acc_reduce.MU0 * Ms / acc_reduce.GAMMA * w * w * (dm[0, 0, 0] ** 2).sum()
        p_mean += p / n
        h.add(torch.as_tensor(m), torch.full((2, 2, 2), p, dtype=torch.float64), th, 1.0, p)
    h.end_cycle(True, 0)
    mk2 = (h.s_mk2 / len(h.kept)).numpy().reshape(K, -1).T
    pk = acc_reduce.harmonic_loss(mk2, alpha, Ms, w)
    assert abs(pk[0].sum() / p_mean - 1) < 1e-5
    p_cell = float(h.s["p"][0, 0, 0]) / h.s_n
    assert abs(p_cell / p_mean - 1) < 1e-12


def test_jump_measures_any_samples_per_cycle(tmp_path):
    """eval_jumps_time reads the samples per cycle from the file (C′ runs with 256 … 1024 samples, 7.75): a smooth
    periodic dissipation gives RF = RN = 0 for 512 samples per cycle."""
    import eval_jumps_time as EJ
    n, cyc = 512, 4
    rows = []
    for c in range(cyc):
        for j in range(n + 1):
            th = 2 * math.pi * j / n
            rows.append([c, 0, 0, np.cos(th), 0, 0, 0, 0, 1.0 + 0.5 * np.cos(2 * th), 0, 0, 0, 0, 10])
    f = tmp_path / "samples_00_amp00.csv"
    np.savetxt(f, np.array(rows), delimiter=",", header="cycle,...")
    P, N60, dM = EJ.series(f)
    assert P.shape == (EJ.KEEP_LAST, n)
    m = EJ.measures(P, N60, dM)
    assert m["RF"] < 1e-9 and m["RN"] < 1e-9


def test_transfer_m_block_average_mask_and_norm():
    """run_loops.transfer_m (--init_transfer, PLAN 4.10): a constant field stays constant; a fine field that is constant
    on each coarse block is recovered exactly, also when the void sub-cells hold other vectors (mask); every output
    vector has norm 1; a coarse cell without a sub-cell of its own index takes the plain block mean and is counted; a
    non-integer cell ratio is rejected."""
    pytest.importorskip("torch")
    from run_loops import transfer_m
    rng = np.random.default_rng(3)
    a, d = box_edge(96.0, 0.65)                                  # a = 144 l_ex, the B0 geometry
    ic = fcc_box(48, d, 3.0)["ids"]
    i_f = fcc_box(96, d, 1.5)["ids"]
    m_c, st = transfer_m(np.broadcast_to([0.0, 0.6, 0.8], (96, 96, 96, 3)), i_f, ic)
    assert np.allclose(m_c, [0.0, 0.6, 0.8]) and st["factor"] == 2
    v = rng.normal(size=(48, 48, 48, 3))
    v /= np.linalg.norm(v, axis=-1, keepdims=True)
    m_f = np.repeat(np.repeat(np.repeat(v, 2, 0), 2, 1), 2, 2)
    m_f[i_f < 0] = rng.normal(size=((i_f < 0).sum(), 3))        # garbage in the void
    m_c, st = transfer_m(m_f, i_f, ic)
    mag = ic >= 0
    assert np.allclose(m_c[mag], v[mag], atol=1e-12)
    assert np.allclose(np.linalg.norm(m_c, axis=-1), 1.0)
    assert st["no_same_index_subcell"] == 0 and st["zero_vector"] == 0 and st["partial_blocks"] > 0
    ids_c = np.zeros((2, 2, 2), np.int8)                         # forced empty block: coarse cell (0, 0, 0) is
    ids_f = np.zeros((4, 4, 4), np.int8)                         # particle 0, its fine sub-cells are all void
    ids_f[:2, :2, :2] = -1
    m_f = np.zeros((4, 4, 4, 3)); m_f[..., 2] = 1.0
    m_f[:2, :2, :2] = [1.0, 0.0, 0.0]
    m_c, st = transfer_m(m_f, ids_f, ids_c)
    assert st["no_same_index_subcell"] == 1 and np.allclose(m_c[0, 0, 0], [1.0, 0.0, 0.0])
    with pytest.raises(ValueError):
        transfer_m(np.zeros((6, 6, 6, 3)), np.zeros((6, 6, 6), np.int8), np.zeros((4, 4, 4), np.int8))


@pytest.mark.parametrize("nc, nf", [(48, 96), (96, 192), (48, 192)])
def test_transfer_m_surfaces_nest_on_study_grids(nc, nf):
    """On the B0 grids (d/l_ex 96, dx 3 / 1.5 / 0.75 l_ex) every coarse particle cell has at least one fine sub-cell of
    the same particle: the staircase surfaces nest, no fallback and no zero vector for a random state."""
    pytest.importorskip("torch")
    from run_loops import transfer_m
    a, d = box_edge(96.0, 0.65)
    ic = fcc_box(nc, d, a / nc)["ids"]
    i_f = fcc_box(nf, d, a / nf)["ids"]
    m_f = np.random.default_rng(nf).normal(size=(nf, nf, nf, 3))
    _, st = transfer_m(m_f, i_f, ic)
    assert st["no_same_index_subcell"] == 0 and st["zero_vector"] == 0


def _conv_recs(F, T15, T3, I15, I3):
    import conv_jobs as CJ
    vals = {"F": F, "T15": T15, "T3": T3, "I15": I15, "I3": I3}
    recs = {}
    for arm, v in vals.items():
        for s, x in zip(CJ.SEEDS, v):
            ok = x is not None
            recs[(arm, s)] = {"lnP": x if ok else float("nan"), "status": "settled" if ok else "no calm tail (< 2 calm cycles)",
                              "s_cyc": 0.03, "dof": 3, "iters": 14000 if arm == "F" else 2000,
                              "transfer": ({"no_same_index_subcell": 0, "zero_vector": 0} if arm[0] == "T" else None),
                              "n_cycles": 6}
    return recs


def test_eval_conv_reading_on_synthetic_data():
    """eval_conv.py (B0, PLAN 4.10): the pre-registered reading on synthetic ln P with a known answer:
    (a) the transfer removes a large state scatter -> go_transfer, Δ(1.5 − 0.75) recovered; (b) the transferred runs
    scatter like independent ones -> stop; (c) the scatter is small with either start -> go_small; (d) two settled
    pairs only -> undetermined; (e) go with fewer than 3 independent pairs -> go_undetermined. The paired statistics equal scipy's one-sample t interval."""
    import eval_conv as EC
    import conv_jobs as CJ
    from scipy import stats as st
    F = np.array([0.0, 0.4, -0.3, 0.2])
    e = np.array([0.02, -0.01, 0.03, -0.02])
    I15 = [0.3, -0.25, 0.25, -0.1]
    I3 = [0.5, 0.1, -0.2, 0.4]
    recs = _conv_recs(list(F), list(F + 0.05 + e), list(F + 0.25 - e), I15, I3)
    text, key = EC.report(recs)
    assert key == "go_transfer"
    lnP = {k: r["lnP"] for k, r in recs.items()}
    p = EC.paired(lnP, "T15", "F", CJ.SEEDS)
    d = 0.05 + e
    lo, hi = st.t.interval(0.90, len(d) - 1, loc=d.mean(), scale=st.sem(d))
    assert abs(p["mean"] - d.mean()) < 1e-12 and abs(p["ci"][0] - lo) < 1e-12 and abs(p["ci"][1] - hi) < 1e-12
    assert p["sd_ci"][0] < p["sd"] < p["sd_ci"][1]
    assert "go:" in text and "inside ± 0.10" in text
    recs = _conv_recs(list(F), [0.35, -0.2, 0.1, -0.3], list(F + 0.25), I15, I3)
    assert EC.report(recs)[1] == "stop"
    recs = _conv_recs(list(F), list(F + 0.05 + e), list(F + 0.25 - e), list(F + 0.04 - e), list(F + 0.3 + e))
    assert EC.report(recs)[1] == "go_small"
    recs = _conv_recs(list(F), [0.05, 0.45, None, None], list(F), I15, I3)
    assert EC.report(recs)[1] == "undetermined"
    recs = _conv_recs(list(F), list(F + 0.05 + e), list(F + 0.25 - e), [0.3, None, None, -0.1], I3)
    assert EC.report(recs)[1] == "go_undetermined"


def test_conv_jobs_change_only_the_planned_arguments():
    """conv_jobs.py (B0): 20 jobs, the F runs first; every line equals the mesh212 template M212_dx15_s931_b9 except
    d_lex, dx_lex, seed, init_transfer, init_wait_h; the T runs start from the F run of the same seed."""
    import shlex
    import conv_jobs as CJ
    from rerun_jobs import get_arg, job_lines
    here = pathlib.Path(__file__).resolve().parent
    tmpl = shlex.split(dict(job_lines(here / "jobs" / "mesh212.txt"))["M212_dx15_s931_b9"])
    keys = lambda tok: {t[2:] for t in tok if t.startswith("--")}
    lines = CJ.lines(str(here / "jobs" / "mesh212.txt"))
    assert len(lines) == 20 and all(n.startswith(CJ.PREFIX + "_F_") for n, _ in lines[:4])
    for n, a in lines:
        tok = shlex.split(a)
        for k in keys(tmpl) | keys(tok):
            if k not in CJ.CHANGED:
                assert get_arg(tok, k) == get_arg(tmpl, k), (n, k)
        arm, s = n.split("_")[1], n.split("_s")[1]
        assert get_arg(tok, "seed") == [s] and get_arg(tok, "d_lex") == ["96.0"]
        assert get_arg(tok, "dx_lex") == [{"F": "0.75", "T15": "1.5", "I15": "1.5", "T3": "3.0", "I3": "3.0"}[arm]]
        assert (get_arg(tok, "init_transfer") == ["%s_F_s%s/init.pt" % (CJ.PREFIX, s)]) == arm.startswith("T")
        assert get_arg(tok, "atol") == [CJ.ATOL]
        assert get_arg(tok, "init_from") is None


def test_eval_conv_run_record_on_real_run_folders():
    """eval_conv.run_record on real run folders (local raw data; skipped without them): the steady rule ≥ 2 calm
    cycles is applied inside run_record and the module state of analyze is restored; P in W/cm³ (M212_dx15_s931_b9:
    4.64 W/cm³, results/mesh212 data); a run with 1 calm cycle at the end (M300_dx15_s922_b9) has no value."""
    import analyze
    import eval_conv as EC
    if not (EC.R / "M212_dx15_s931_b9" / "DONE").exists() or not (EC.R / "M300_dx15_s922_b9" / "DONE").exists():
        pytest.skip("local run folders not present")
    old = analyze.STEADY_TAIL, analyze.STEADY_MIN
    r = EC.run_record(EC.R / "M212_dx15_s931_b9")
    assert r["status"] == "settled" and abs(math.exp(r["lnP"]) / 4.64 - 1) < 0.005
    assert abs(r["b_ratio"] - 1) < 0.01 and r["iters"] == 13500
    assert EC.run_record(EC.R / "M300_dx15_s922_b9")["status"].startswith("no calm tail")
    assert (analyze.STEADY_TAIL, analyze.STEADY_MIN) == old
    assert EC.run_record(EC.R / "no_such_run")["status"] == "not run"


def test_cleanup_protects_init_transfer_sources(tmp_path):
    """cleanup_instance.dependents knows the sources of --init_transfer (B0 F states) as well as of --init_from."""
    import cleanup_instance as CI
    (tmp_path / "a.txt").write_text("B0_T15_s941 --d_lex 96 --init_transfer B0_F_s941/init.pt\n"
                                    "X_b50 --d_lex 96 --init_from X_b9/init.pt\n")
    dep = CI.dependents(tmp_path)
    assert dep["B0_F_s941/init.pt"] == {"B0_T15_s941"} and dep["X_b9/init.pt"] == {"X_b50"}


def test_eval_sweep_reading_on_synthetic_resonances():
    """eval_sweep.py (PROTOKOLL 7.89): the pre-registered reading on synthetic loss-per-cycle curves w(f) of a driven
    damped oscillator: a lower eigen-frequency at A 10 -> supports; the same eigen-frequency -> against; both far
    above the range with a slow rise -> against; the vertex of the parabola recovers the peak within one grid step."""
    import eval_sweep as ES
    f = np.array(ES.F_MHZ)

    def w(f0, zeta=0.15):                          # loss per cycle ∝ χ'' of a damped oscillator
        r = f / f0
        return 2 * zeta * r / ((1 - r ** 2) ** 2 + (2 * zeta * r) ** 2)
    key, f10, f20 = ES.reading(f, w(35.0), w(150.0))
    assert key == "supports" and abs(math.log(f10 / 35.0)) < 0.35 and abs(math.log(f20 / 150.0)) < 0.35
    assert ES.reading(f, w(60.0), w(60.0))[0] == "against"
    slow = lambda: f ** 0.8                        # n = 1.8 everywhere, no maximum
    assert ES.reading(f, slow(), slow())[0] == "against"
    assert ES.reading(f, w(35.0), w(1000.0))[0] == "supports"          # A 20 still rising at 240 MHz


def test_eval_atol60_reading_on_synthetic_series():
    """PROTOKOLL 7.100: the classes and the reading of the long tolerance test (fixed before the runs)."""
    import eval_atol60 as EA
    n = 38                                                        # windows ending at cycles 2 … 39
    rn_ok, clo_ok, w_ok = [0.01] * n, [0.0002] * n, [1.0] * n
    s, c = EA.windows([0.2] * 5 + rn_ok[5:], [0.005] * 5 + clo_ok[5:], [1.3] * 5 + w_ok[5:])
    assert c == "periodic" and s == 7                             # first good window ends at cycle 2 + 5
    s, c = EA.windows([0.2] * n, [0.004] * n, w_ok)
    assert c == "stays irregular" and s is None
    s, c = EA.windows([0.07] * n, [0.0015] * n, w_ok)
    assert c == "marginal"
    early, late, irr = (5, "periodic"), (20, "periodic"), (None, "stays irregular")
    assert EA.reading([early, early, irr], [irr, irr, early]) == "tolerance artefact"
    assert EA.reading([late, late, early], [late, early, late]) == "settles slowly"
    assert EA.reading([irr, irr, early], [irr, irr, irr]) == "stays irregular"
    assert EA.reading([early, irr, (None, "marginal")], [early, early, early]) == "not determined"
    assert EA.value_reading([10.0, 10.5, 9.5], [10.2, 9.8, 10.0], [10.0, 10.0, 10.0]) == "the C′ value holds"
    assert EA.value_reading([13.0, 13.5, 9.5], [13.2, 13.8, 10.0], [10.0, 10.0, 10.0]) == "the C′ value is off"


def test_init_from_jobs_have_the_init_keys_of_their_source():
    """run_loops stops a job whose --init_from state was relaxed with other INIT_KEYS (atol since 7.88). Every job line
    with --init_from <run>/init.pt must have the INIT_KEYS of the line of <run> (PROTOKOLL 7.101: the L6 lines of the
    first jobs/atol60.txt had atol 1e-6 against 1e-5 and would have stopped at once)."""
    import shlex
    import run_loops
    from rerun_jobs import job_lines
    jobs = pathlib.Path(__file__).resolve().parent / "jobs"
    src = {}
    for f in sorted(jobs.glob("*.txt")):
        for name, a in job_lines(f):
            assert src.setdefault(name, a) == a, "two different lines for %s" % name
    parse = lambda a: vars(run_loops.parse_args(shlex.split(a) + ["--out", "x"]))
    n = 0
    for f in sorted(jobs.glob("*.txt")):
        for name, a in job_lines(f):
            t = parse(a)
            init = t.get("init_from") or ""
            run = init.split("/")[0]
            if not init.endswith("/init.pt") or run not in src:
                continue                                          # no init_from, or a state file without a job line
            s = parse(src[run])
            diff = [k for k in run_loops.INIT_KEYS if k != "relax_maxiter" and t.get(k) != s.get(k)]
            assert not diff, "%s (%s): %s differ from %s" % (name, f.name, diff, run)
            n += 1
    assert n >= 190


def test_atol_job_lines():
    """PROTOKOLL 7.101: atol60 L5 starts from K1_k/init.pt, L6 relaxes itself; the atol30 twins are the production line
    of K1_k with only atol 1e-6 added and the snapshots removed. The job files equal the generator output."""
    import shlex
    import pair_jobs
    from rerun_jobs import job_lines, strip_snap
    here = pathlib.Path(__file__).resolve().parent
    for fn, lines in (("atol60.txt", pair_jobs.atol60_lines()), ("atol30.txt", pair_jobs.atol30_lines())):
        assert [(n, " ".join(t)) for n, t in lines] == job_lines(here / "jobs" / fn)
    for name, tok in pair_jobs.atol60_lines():
        assert ("--init_from" in tok) == name.endswith("L5")
    src = dict(job_lines(here / "jobs" / "sobol_corners.txt"))
    for k, (name, tok) in enumerate(pair_jobs.atol30_lines()):
        prod = shlex.split(src["K1_%d" % k])
        assert name == "K1_%d_L6" % k and tok[-2:] == ["--atol", "1e-06"]
        assert tok[:-2] == strip_snap(prod) and "--init_from" not in tok
        assert {x for x in prod if x.startswith("--")} - set(tok) == {"--snap_phases", "--snap_cycles", "--snap_mT"}


def test_eval_atol30_reading_and_closure_median(tmp_path, monkeypatch):
    """PROTOKOLL 7.101: the reading of the tolerance twins (fixed before the runs) and the measured cycles."""
    import json
    import eval_atol30 as E3
    assert E3.reading([0.9, 0.7, 0.2]) == "physical"
    assert E3.reading([0.1, 0.3, 0.8]) == "numerical"
    assert E3.reading([0.6, 0.55, None]) == "physical"
    assert E3.reading([0.9, None, 0.1]) == "not determined"
    assert E3.reading([None, None, None]) == "not determined"
    cyc = lambda clo: {"closure": clo, "w_loop": 1e-6, "w_dis": 1.02e-6}
    corrected = {"h_amp_initial": 1.0, "cycles": [cyc(0.5), cyc(0.4)] + [cyc(0.001)] * 5}
    plain = {"cycles": [cyc(0.5)] + [cyc(0.003)] * 6}
    assert len(E3.measured(corrected)) == 5 and len(E3.measured(plain)) == 6
    assert E3.closure_median([corrected, plain]) == 0.003
    st = [dict(corrected, measure=True, complete=True) for _ in E3.AMPS]
    for name, clo in (("K1_0", 0.001), ("K1_0_L6", 0.0008)):
        d = tmp_path / name
        d.mkdir()
        for s in st:
            s["cycles"] = [cyc(0.5), cyc(0.4)] + [cyc(clo)] * 5
        (d / "summary.json").write_text(json.dumps({"stages": st, "init": {"E_over_Kd": 5e-4, "M_over_Ms": [0, 0, 0]}}))
        (d / "DONE").touch()
    monkeypatch.setattr(E3, "R", tmp_path)
    a, b = E3.load(tmp_path / "K1_0"), E3.load(tmp_path / "K1_0_L6")
    assert abs(E3.closure_median(b[0]) / E3.closure_median(a[0]) - 0.8) < 1e-12
    assert E3.load(tmp_path / "K1_1") is None


def test_eval_closure_state_same_sign_references():
    """PROTOKOLL 7.101: consecutive ΔM of independent values of M have the same sign in 1/3 of the pairs, a wandering
    state (independent steps) in 1/2, a monotonic change in all pairs; the stage bootstrap gives a finite SE."""
    import eval_closure_state as EC
    rng = np.random.default_rng(3)
    def share(m_series):
        st = [EC.same_sign_pairs(list(np.sign(np.diff(m)))) for m in m_series]
        return EC.boot_share(st, np.random.default_rng(1), n_boot=200)
    p, se, sb, n = share([rng.normal(size=6) for _ in range(3000)])
    assert abs(p - 1 / 3) < 0.015 and 0 < se < 0.02 and n == 3000 * 4
    p, _, _, _ = share([np.cumsum(rng.normal(size=6)) for _ in range(3000)])
    assert abs(p - 0.5) < 0.015
    p, _, _, _ = share([np.arange(6) * 10.0 + rng.normal(size=6) for _ in range(300)])
    assert p > 0.99
    assert EC.same_sign_pairs([1, 1, -1, -1, 1]) == (2, 4) and EC.same_sign_pairs([1]) == (0, 0)
    assert EC.reading_drift(-0.005, 0.002, 0.33, 0.024) == "recurrent jumps"
    assert EC.reading_drift(-0.05, 0.002, 0.33, 0.024) == "still relaxing"
    assert EC.reading_drift(0.0, 0.002, 0.62, 0.024) == "still relaxing"


def test_state_files_are_on_the_disk_before_the_rename(tmp_path, monkeypatch):
    """PROTOKOLL 7.102: a hard restart of the instance left a renamed init.pt with 0 bytes. save_atomic, write_atomic
    and run_queue.save_queue_state flush the data (os.fsync of the file) before the rename and the directory after it."""
    import json
    torch = pytest.importorskip("torch")
    import run_loops
    import run_queue
    calls = []
    real = os.fsync
    monkeypatch.setattr(os, "fsync", lambda fd: (calls.append(fd), real(fd))[1])
    p = tmp_path / "init.pt"
    run_loops.save_atomic({"m": torch.ones(2)}, p)
    assert float(torch.load(p)["m"].sum()) == 2.0 and len(calls) == 2 and not (tmp_path / "init.pt.tmp").exists()
    run_loops.write_atomic(tmp_path / "summary.json", "{}")
    assert (tmp_path / "summary.json").read_text() == "{}" and len(calls) == 4
    run_queue.save_queue_state(tmp_path / "run", {"state": "done"})
    assert json.loads((tmp_path / "run" / "QUEUE_STATE.json").read_text()) == {"state": "done"} and len(calls) == 6


def test_eval_b0b_check_attempts_of_a_restarted_run(tmp_path):
    """PROTOKOLL 7.102: eval_b0b_check.attempts splits the cycle lines of stdout.log at each start of run_loops."""
    import eval_b0b_check as EB
    head = "2026-10-14 08:37:13  magnum.np:INFO magnum.np 2.2.0 (run_loops.py ...)"
    cyc = "[amp00 c%d] h=0.00075  b_peak=0.006  w_loop=%s  w_dis=%s  closure=%s  dW=-  max_angle=57 deg"
    lines = [head, "[relax] it 500", cyc % (0, "4.289e-07", "8.636e-07", "0.032"), cyc % (1, "5.967e-07", "1.13e-06", "0.011"),
             head, cyc % (0, "4.289e-07", "8.7e-07", "0.032")]
    (tmp_path / "stdout.log").write_text("\n".join(lines) + "\n")
    at = EB.attempts(tmp_path)
    assert len(at) == 2 and sorted(at[0]) == [0, 1] and sorted(at[1]) == [0]
    assert at[0][0] == (4.289e-07, 8.636e-07, 0.032) and at[1][0][1] == 8.7e-07


def test_long60_jobs_and_readings():
    """PROTOKOLL 7.103: rerun_jobs.cond finds a tested amplitude by its place in the nominal list (amps may start above
    9 mT), the long60 lines equal the C′ lines up to b_list / n_amp / cycles_last_amp, and the fixed readings."""
    import shlex
    import analyze_design as AD
    import eval_long60 as EL
    import pair_jobs
    from rerun_jobs import NOMINAL_MT, get_arg, job_lines, set_arg
    assert tuple(NOMINAL_MT) == tuple(int(x) for x in AD.AMPS_MT)
    here = pathlib.Path(__file__).resolve().parent
    src = dict(job_lines(here / "jobs" / "freq_K1pair.txt"))
    lines = pair_jobs.long60_lines()
    assert [(n, a) for n, a in lines] == job_lines(here / "jobs" / "long60.txt")
    assert sorted(n for n, _ in lines) == sorted(["cond_K1A10_%d_f2_b9" % k for k in range(3)] +
                                                 ["cond_K1_%d_f2_b%d" % (k, m) for k in range(3) for m in (20, 35)])
    for name, args in lines:
        base, mT = name[5:].rsplit("_b", 1)
        tok, ref = shlex.split(args), shlex.split(src[base])
        k = NOMINAL_MT.index(int(mT))
        assert get_arg(tok, "b_list") == get_arg(ref, "b_list")[:k + 1] and get_arg(tok, "n_amp") == [str(k + 1)]
        assert get_arg(tok, "cycles_last_amp") == ["40"] and get_arg(tok, "init_from") is None
        strip = lambda t: set_arg(set_arg(set_arg(list(t), "b_list", None), "n_amp", None), "cycles_last_amp", None)
        assert sorted(strip(tok)) == sorted(strip(ref))
    assert EL.motion_reading(["periodic", "periodic", "stays irregular"]) == "periodic"
    assert EL.motion_reading(["stays irregular", "marginal", "stays irregular"]) == "stays irregular"
    assert EL.motion_reading(["periodic", "marginal", None]) == "not determined"
    assert EL.value_reading([10.0, 10.5, 12.0], [10.0, 10.0, 10.0]) == "the C′ value holds"
    assert EL.value_reading([13.0, 13.5, 10.0], [10.0, 10.0, 10.0]) == "the C′ value is off"
    assert EL.value_reading([10.0, None, 13.0], [10.0, 10.0, 10.0]) == "not determined"
    assert EL.sign_reading([-0.35, -0.36, -0.31])[0] == "the sign change holds"
    assert EL.sign_reading([0.10, -0.05, 0.02])[0] == "the sign change does not hold"
    assert EL.sign_reading([-0.6, 0.1, -0.2])[0] == "not determined"
    assert EL.sign_reading([-0.35, -0.36, None])[0] == "not determined (runs missing)"
    assert EL.history_class(0.0, True, False) == "bit-identical"
    assert EL.history_class(3e-6, True, False) == "same within the solver tolerance"
    assert EL.history_class(3e-6, False, False) == "differs"
    assert EL.history_class(2e-3, True, True) == "differs, restart"


def test_newell_pbc_compensation_exact_and_fd_off_axis_error():
    """Demag solver check (PROTOKOLL 7.105): uniform m in an aligned cube of edge a/2 in a periodic box has the exact
    energy f (1 - f) / 3 K_d per box volume (f = 1/8, every direction of m: cubic symmetry, zero mean field). The Newell
    demag with one image shell and compensated outer surface charges (check_mesh_energy.newell_pbc) gives it for every
    direction; DemagFieldPBC gives it for m along an axis and deviates for m along (111) (+13.7 % for 4 cells per cube
    edge; +4.6 / 1.4 / 0.43 % for 8 / 16 / 32 cells)."""
    import torch
    import check_mesh_energy as C
    from magnumnp import Mesh, State
    from magnumnp.field_terms.demagPBC import DemagFieldPBC
    from units import units
    old = torch.get_default_dtype()
    torch.set_default_dtype(torch.float64)
    try:
        U = units(1.5, 1e-11)
        Ms, Kd, dx, n = U["Ms"], U["Kd"], 1e-9, 8
        ms = torch.zeros((n, n, n, 1))
        ms[:n // 2, :n // 2, :n // 2] = Ms
        exact = (1 / 8) * (7 / 8) / 3
        st_n, term_n = C.newell_pbc(n, dx, ms, 1)
        st_f = State(Mesh((n, n, n), (dx,) * 3))
        st_f.material = {"Ms": ms}
        fd = DemagFieldPBC()
        for u, fd_exact in (((1.0, 0.0, 0.0), True), ((1.0, 1.0, 1.0), False)):
            m = torch.tensor(u).expand(n, n, n, 3) / math.sqrt(sum(x * x for x in u))
            assert abs(C.energy_newell_pbc(st_n, term_n, m.clone(), Ms, Kd) / exact - 1) < 1e-3
            st_f.m = m.clone()
            e_fd = float(fd.E(st_f)) / (Kd * dx ** 3 * n ** 3)
            assert (abs(e_fd / exact - 1) < 1e-9) if fd_exact else (e_fd / exact - 1 > 0.02)
    finally:
        torch.set_default_dtype(old)


def test_strayfield_setup_and_init_keys():
    """--strayfield (PROTOKOLL 7.105): fd gives DemagFieldPBC on a mesh with pbc 1; newell gives DemagField with
    true_pbc and --strayfield_shells image shells. Both keys set the relaxed state (INIT_KEYS); an old config without
    them counts as fd with 4 shells (NEW_KEYS), so old runs resume and old init states stay usable for fd runs."""
    import types
    import run_loops
    from magnumnp import DemagField, DemagFieldPBC
    pbc, term = run_loops.strayfield_setup(types.SimpleNamespace(strayfield="fd", strayfield_shells=4))
    assert pbc == (1, 1, 1) and isinstance(term, DemagFieldPBC)
    pbc, term = run_loops.strayfield_setup(types.SimpleNamespace(strayfield="newell", strayfield_shells=2))
    assert pbc == (2, 2, 2) and isinstance(term, DemagField) and term._true_pbc
    with pytest.raises(ValueError):
        run_loops.strayfield_setup(types.SimpleNamespace(strayfield="newell", strayfield_shells=0))
    assert {"strayfield", "strayfield_shells"} <= set(run_loops.INIT_KEYS)
    assert run_loops.NEW_KEYS["strayfield"] == "fd" and run_loops.NEW_KEYS["strayfield_shells"] == 4
    a = run_loops.parse_args(["--out", "x", "--d_lex", "96"])
    assert a.strayfield == "fd" and a.strayfield_shells == 4
    a = run_loops.parse_args(["--out", "x", "--d_lex", "96", "--strayfield", "newell", "--strayfield_shells", "2"])
    assert a.strayfield == "newell" and a.strayfield_shells == 2


def test_true_pbc_equals_the_solver_check_helper():
    """DemagField(true_pbc=True) of magnum.np gives the same energy as check_mesh_energy.energy_newell_pbc, which made
    the solver check of PROTOKOLL 7.105 (random m, partly filled box, one image shell)."""
    import torch
    import check_mesh_energy as C
    from magnumnp import Mesh, State, DemagField, constants
    old = torch.get_default_dtype()
    torch.set_default_dtype(torch.float64)
    try:
        Ms = 1. / constants.mu_0
        Kd = 0.5 * constants.mu_0 * Ms ** 2
        ms = torch.full((8, 8, 8, 1), 1e-8 * Ms)
        ms[1:6, 2:7, 0:5] = Ms
        g = torch.Generator().manual_seed(1)
        r = torch.randn((8, 8, 8, 3), generator=g)
        m = r / r.norm(dim=-1, keepdim=True)
        st = State(Mesh((8, 8, 8), (1e-9,) * 3, pbc=(1, 1, 1)))
        st.material = {"Ms": ms}
        st.m = m
        e1 = float(DemagField(true_pbc=True).E(st)) / (Kd * 512e-27)
        e2 = C.energy_newell_pbc(*C.newell_pbc(8, 1e-9, ms, 1), m, Ms, Kd)
        assert e1 == pytest.approx(e2, rel=1e-12)
    finally:
        torch.set_default_dtype(old)
