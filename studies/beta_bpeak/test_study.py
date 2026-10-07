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
