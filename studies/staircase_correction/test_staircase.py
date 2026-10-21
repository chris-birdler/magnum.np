"""Unit tests T0 of PLAN.md section 5 (code tests on small boxes, not physics checks)."""
import os
import sys
import pathlib

import numpy as np
import pytest
import torch

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "beta_bpeak"))
os.environ.setdefault("CUDA_DEVICE", "-1")

from magnumnp import set_precision  # noqa: E402
set_precision("double")

from geometry import fcc_box  # noqa: E402
from nearfield import NearCorrection, FFTField, v1_cells, fine_owner_m, coarse_mean, centroid_offsets, centroid_field  # noqa: E402
from magnumnp.field_terms.demag_staircase import centroid_field_fourier, pair_matrix, sym_matvec  # noqa: E402


def one_sphere(N, q, r, centre):
    """Fine body of one sphere (radius r in coarse cells, centre in coarse units) in an open box of N^3 cells."""
    c = (np.arange(q * N) + 0.5) / q
    X, Y, Z = np.meshgrid(c, c, c, indexing="ij")
    inside = (X - centre[0]) ** 2 + (Y - centre[1]) ** 2 + (Z - centre[2]) ** 2 <= r ** 2
    return np.where(inside, 0, -1).astype(np.int8)


def random_m(shape, seed):
    g = np.random.default_rng(seed)
    m = g.normal(size=shape + (3,))
    return torch.as_tensor(m / np.linalg.norm(m, axis=-1, keepdims=True))


def test_no_modified_cells_gives_zero():
    N, q = 6, 2
    ids = -np.ones((q * N,) * 3, dtype=np.int8)
    ids[2:8, 4:10, 2:6] = 0                                   # whole coarse cells only
    n, pid, occ = v1_cells(ids, q)
    nc = NearCorrection(n, occ, q, p=2, periodic=False)
    m = random_m((N, N, N), 1)
    assert len(nc.J) == 0
    assert torch.count_nonzero(nc.h_pair(m)) == 0 and torch.count_nonzero(nc.h_local(m)) == 0


@pytest.mark.parametrize("periodic", [False, True])
def test_correction_is_symmetric(periodic):
    N, q = 8, 2
    ids = one_sphere(N, q, 2.7, (4.0, 4.0, 4.0))
    n, pid, occ = v1_cells(ids, q)
    nc = NearCorrection(n, occ, q, p=2, periodic=periodic)
    m1, m2 = random_m((N, N, N), 2), random_m((N, N, N), 3)
    v = nc.v.reshape(N, N, N, 1)
    for h in (nc.h_pair, nc.h_local):
        a = float((v * m1 * h(m2)).sum())
        b = float((v * m2 * h(m1)).sum())
        assert abs(a - b) <= 1e-12 * (abs(a) + abs(b))


def test_full_window_open_boundaries_equals_fine_reference():
    N, q = 8, 2
    ids = one_sphere(N, q, 2.9, (4.0, 4.3, 3.8))
    n, pid, occ = v1_cells(ids, q)
    nc = NearCorrection(n, occ, q, p=N - 1, periodic=False)
    m = random_m((N, N, N), 4)
    m = torch.where(torch.as_tensor(n > 0)[..., None], m, torch.zeros_like(m))
    h_c = FFTField(N, "open").h(nc.v.reshape(N, N, N), m) + nc.h_pair(m)
    body = ids >= 0
    h_f = FFTField(q * N, "open").h(body.astype(float), fine_owner_m(m, q))
    h_ref = coarse_mean(h_f, body, q)
    mag = torch.as_tensor(n > 0)
    err = (h_c - h_ref)[mag].abs().max()
    assert err < 1e-10, float(err)


def test_shift_invariance_under_pbc():
    N, q = 16, 2
    geo = fcc_box(q * N, 6.0, 0.5)                            # 4 spheres, d 6 cells, fine cells of 0.5
    ids = geo["ids"]
    m = random_m((N, N, N), 5)
    s = (5, 3, 7)
    E = []
    for shift in ((0, 0, 0), s):
        idsr = np.roll(ids, tuple(q * k for k in shift), axis=(0, 1, 2))
        mr = torch.roll(m, shift, dims=(0, 1, 2))
        n, pid, occ = v1_cells(idsr, q)
        nc = NearCorrection(n, occ, q, p=2, periodic=True)
        h = FFTField(N, ("newell", 1)).h(nc.v.reshape(N, N, N), mr) + nc.h_pair(mr)
        E.append(nc.energy(mr, h))
    assert abs(E[0] - E[1]) <= 1e-10 * abs(E[0])


def test_uniform_m_cubic_body_exact_energy():
    N, q = 16, 2
    ids = fcc_box(q * N, 6.0, 0.5)["ids"]
    body = ids >= 0
    phi = body.mean()
    m = torch.zeros((q * N,) * 3 + (3,), dtype=torch.float64)
    m[..., :] = torch.tensor([1.0, 1.0, 1.0]) / np.sqrt(3.0)
    h = FFTField(q * N, ("newell_ext", 4)).h(body.astype(float), m)
    E = float(-(torch.as_tensor(body)[..., None] * m * h).sum() / (q * N) ** 3)
    exact = phi * (1 - phi) / 3
    assert abs(E / exact - 1) < 5e-4, (E, exact)


@pytest.mark.parametrize("form,periodic", [("pair", True), ("pair", False), ("local", True)])
def test_field_term_equals_correction(form, periodic):
    """DemagStaircaseCorrection (the field term of magnum.np; the pair form as one sparse matrix, PROTOKOLL 3.36) gives
    M_s times the correction of NearCorrection (gather and index_add per offset)."""
    from magnumnp import Mesh, State, DemagStaircaseCorrection
    N, q, Ms = 16, 2, 8e5
    ids = fcc_box(q * N, 6.0, 0.5)["ids"]
    n, pid, occ = v1_cells(ids, q)
    nc = NearCorrection(n, occ, q, p=2, periodic=periodic)
    m = random_m((N, N, N), 6)
    state = State(Mesh((N, N, N), (2e-9,) * 3, pbc=(1, 1, 1)))
    state.m = m
    state.material["Ms"] = torch.as_tensor(n / q ** 3 * Ms).reshape(N, N, N, 1)
    term = DemagStaircaseCorrection(n, occ, q, 2, Ms, form=form, periodic=periodic)
    h = term.h(state)
    ref = Ms * (nc.h_pair(m) if form == "pair" else nc.h_local(m))
    assert float((h - ref).abs().max()) <= 1e-12 * Ms
    E = float(term.E(state))
    E_ref = -0.5 * 4e-7 * np.pi * float((state.material["Ms"] * m * ref).sum()) * 8e-27
    assert abs(E - E_ref) <= 1e-9 * abs(E_ref)


@pytest.mark.parametrize("q", [2, 4])
def test_centroid_full_window_open_boundaries_equals_fine_reference(q):
    """With the first-moment far field (centroid) and a window over the whole box the correction is still exact."""
    N = 8
    ids = one_sphere(N, q, 2.9, (4.0, 4.3, 3.8))
    n, pid, occ = v1_cells(ids, q)
    nc = NearCorrection(n, occ, q, p=N - 1, periodic=False, centroid=True)
    m = random_m((N, N, N), 7)
    m = torch.where(torch.as_tensor(n > 0)[..., None], m, torch.zeros_like(m))
    co = FFTField(N, "open")
    v = nc.v.reshape(N, N, N)
    h0 = co.h(v, m)
    hc = centroid_field(co.term._N, v, centroid_offsets(occ, q), m, shape=co.fft_shape)
    h_c = h0 + hc + nc.h_pair(m)
    body = ids >= 0
    h_ref = coarse_mean(FFTField(q * N, "open").h(body.astype(float), fine_owner_m(m, q)), body, q)
    mag = torch.as_tensor(n > 0)
    assert float((h_c - h_ref)[mag].abs().max()) < 1e-10


def test_centroid_operator_is_symmetric():
    N, q = 16, 2
    ids = fcc_box(q * N, 6.0, 0.5)["ids"]
    n, pid, occ = v1_cells(ids, q)
    nc = NearCorrection(n, occ, q, p=2, periodic=True, centroid=True)
    co = FFTField(N, ("newell", 1))
    v = nc.v.reshape(N, N, N)
    d = centroid_offsets(occ, q)
    m1, m2 = random_m((N, N, N), 8), random_m((N, N, N), 9)
    co.h(v, m1)                                                    # sets up the kernel

    def h(m):
        return co.h(v, m) + centroid_field(co.term._N, v, d, m) + nc.h_pair(m)

    a = float((v[..., None] * m1 * h(m2)).sum())
    b = float((v[..., None] * m2 * h(m1)).sum())
    assert abs(a - b) <= 1e-11 * (abs(a) + abs(b))


@pytest.mark.parametrize("form", ["pair", "local"])
def test_field_term_with_centroid_equals_static_model(form):
    """DemagStaircaseCorrection(centroid=True) gives M_s times (window term + centroid far-field term)."""
    from magnumnp import Mesh, State, DemagField, DemagStaircaseCorrection
    N, q, Ms = 16, 2, 8e5
    ids = fcc_box(q * N, 6.0, 0.5)["ids"]
    n, pid, occ = v1_cells(ids, q)
    nc = NearCorrection(n, occ, q, p=2, periodic=True, centroid=True)
    m = random_m((N, N, N), 10)
    state = State(Mesh((N, N, N), (2e-9,) * 3, pbc=(1, 1, 1)))
    state.m = m
    state.material["Ms"] = torch.as_tensor(n / q ** 3 * Ms).reshape(N, N, N, 1)
    dem = DemagField(true_pbc=True)
    dem.h(state)
    term = DemagStaircaseCorrection(n, occ, q, 2, Ms, form=form, centroid=True, demag=dem)
    v = torch.as_tensor(n / q ** 3)
    ref = Ms * ((nc.h_pair(m) if form == "pair" else nc.h_local(m))
                + centroid_field(dem._N, v, centroid_offsets(occ, q), m))
    assert float((term.h(state) - ref).abs().max()) <= 1e-10 * Ms


def test_fine_term_open_equals_full_window_correction():
    """DemagStaircaseFine (exact field of the fine staircase, PROTOKOLL 3.18) equals DemagField with the moments v M_s
    plus the pairwise correction with the full window, for open boundaries."""
    from magnumnp import Mesh, State, DemagStaircaseFine
    N, q, Ms = 8, 2, 8e5
    ids = one_sphere(N, q, 2.9, (4.0, 4.3, 3.8))
    n, pid, occ = v1_cells(ids, q)
    nc = NearCorrection(n, occ, q, p=N - 1, periodic=False)
    mag = torch.as_tensor(n > 0)
    m = torch.where(mag[..., None], random_m((N, N, N), 11), torch.zeros(N, N, N, 3, dtype=torch.float64))
    state = State(Mesh((N, N, N), (2e-9,) * 3))
    state.m = m
    state.material["Ms"] = torch.as_tensor(n / q ** 3 * Ms).reshape(N, N, N, 1)
    term = DemagStaircaseFine(ids >= 0, q, Ms)
    ref = Ms * (FFTField(N, "open").h(nc.v.reshape(N, N, N), m) + nc.h_pair(m))
    assert float((term.h(state) - ref)[mag].abs().max()) <= 1e-10 * Ms


def test_fine_term_energy_and_reciprocity_under_pbc():
    """With PBC: the field is the fine field averaged over the body sub-cells (static reference of check_static), the
    energy is the stray-field energy of the fine body, and the operator is symmetric in the weights v."""
    from magnumnp import Mesh, State, DemagStaircaseFine
    N, q, Ms = 16, 2, 8e5
    ids = fcc_box(q * N, 6.0, 0.5)["ids"]
    body = ids >= 0
    n, pid, occ = v1_cells(ids, q)
    v = torch.as_tensor(n / q ** 3).reshape(N, N, N, 1)
    state = State(Mesh((N, N, N), (2e-9,) * 3, pbc=(1, 1, 1)))
    state.material["Ms"] = v * Ms
    term = DemagStaircaseFine(body, q, Ms, true_pbc=True)
    m1, m2 = random_m((N, N, N), 12), random_m((N, N, N), 13)
    state.m = m1
    h1 = term.h(state)
    fine = FFTField(q * N, ("newell", 1))
    hf = fine.h(body.astype(np.float64), fine_owner_m(m1, q))
    assert float((h1 - Ms * coarse_mean(hf, body, q))[torch.as_tensor(n > 0)].abs().max()) <= 1e-10 * Ms
    E_f = -0.5 * 4e-7 * np.pi * Ms ** 2 * float((torch.as_tensor(body)[..., None] * fine_owner_m(m1, q) * hf).sum()) * 1e-27
    assert abs(float(term.E(state)) - E_f) <= 1e-9 * abs(E_f)
    state.m = m2
    h2 = term.h(state)
    a, b = float((v * m1 * h2).sum()), float((v * m2 * h1).sum())
    assert abs(a - b) <= 1e-11 * (abs(a) + abs(b))


def test_centroid_pbc_differences_equal_fourier():
    """Under PBC the centroid term with central differences in real space (12 FFTs) equals the form with i sin w in
    Fourier space (24 FFTs) (PROTOKOLL 3.36, criterion C2 of PLAN 10)."""
    N, q = 16, 4
    ids = fcc_box(q * N, 6.0, 0.25)["ids"]
    n, pid, occ = v1_cells(ids, q)
    co = FFTField(N, ("newell", 1))
    v = torch.as_tensor(n / q ** 3)
    d = centroid_offsets(occ, q)
    m = random_m((N, N, N), 11)
    co.h(v, m)                                                     # sets up the kernel
    new = centroid_field(co.term._N, v, d, m)
    old = centroid_field_fourier(co.term._N, v, d, m)
    assert float((new - old).abs().max()) <= 1e-12 * float(old.abs().max())


def test_pair_matrix_equals_index_add_and_sums_duplicates():
    """pair_matrix gives the pair sum of gather, sym_matvec and index_add, also with repeated (target, source) pairs."""
    g = torch.Generator().manual_seed(12)
    n_cells, P = 50, 400
    tgt = torch.randint(0, n_cells, (P,), generator=g)
    src = torch.randint(0, n_cells, (P,), generator=g)
    tgt, src = torch.cat([tgt, tgt[:30]]), torch.cat([src, src[:30]])          # 30 repeated pairs
    cf = torch.randn(len(tgt), 6, generator=g, dtype=torch.float64)
    m = torch.randn(n_cells, 3, generator=g, dtype=torch.float64)
    ref = torch.zeros_like(m).index_add_(0, tgt, sym_matvec(cf, m[src]))
    A = pair_matrix(tgt, src, cf, n_cells)
    new = (A @ m.reshape(-1)).reshape(-1, 3)
    assert float((new - ref).abs().max()) <= 1e-12 * float(ref.abs().max())


def test_pair_matrix_builds_on_the_cpu_with_another_default_device():
    """magnum.np sets the default device to the GPU; pair_matrix must still build on the CPU (audit 2026-10-21).
    The meta device stands for the GPU: any tensor made on the default device would fail to mix with the CPU data."""
    g = torch.Generator().manual_seed(13)
    n_cells, P = 40, 300
    tgt = torch.randint(0, n_cells, (P,), generator=g)
    src = torch.randint(0, n_cells, (P,), generator=g)
    tgt, src = torch.cat([tgt, tgt[:20]]), torch.cat([src, src[:20]])
    cf = torch.randn(len(tgt), 6, generator=g, dtype=torch.float64)
    m = torch.randn(n_cells, 3, generator=g, dtype=torch.float64)
    ref = (pair_matrix(tgt, src, cf, n_cells) @ m.reshape(-1)).reshape(-1, 3)
    torch.set_default_device("meta")
    try:
        A = pair_matrix(tgt, src, cf, n_cells)
    finally:
        torch.set_default_device("cpu")
    assert A.device.type == "cpu"
    assert float(((A @ m.reshape(-1)).reshape(-1, 3) - ref).abs().max()) == 0.0
