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


@pytest.mark.parametrize("form", ["pair", "local"])
def test_field_term_equals_correction(form):
    """DemagStaircaseCorrection (the field term of magnum.np) gives M_s times the correction of NearCorrection."""
    from magnumnp import Mesh, State, DemagStaircaseCorrection
    N, q, Ms = 16, 2, 8e5
    ids = fcc_box(q * N, 6.0, 0.5)["ids"]
    n, pid, occ = v1_cells(ids, q)
    nc = NearCorrection(n, occ, q, p=2, periodic=True)
    m = random_m((N, N, N), 6)
    state = State(Mesh((N, N, N), (2e-9,) * 3, pbc=(1, 1, 1)))
    state.m = m
    state.material["Ms"] = torch.as_tensor(n / q ** 3 * Ms).reshape(N, N, N, 1)
    term = DemagStaircaseCorrection(n, occ, q, 2, Ms, form=form)
    h = term.h(state)
    ref = Ms * (nc.h_pair(m) if form == "pair" else nc.h_local(m))
    assert float((h - ref).abs().max()) <= 1e-12 * Ms
    E = float(term.E(state))
    E_ref = -0.5 * 4e-7 * np.pi * float((state.material["Ms"] * m * ref).sum()) * 8e-27
    assert abs(E - E_ref) <= 1e-9 * abs(E_ref)


def test_centroid_full_window_open_boundaries_equals_fine_reference():
    """With the first-moment far field (centroid) and a window over the whole box the correction is still exact."""
    N, q = 8, 2
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
