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

from geometry import fcc_box, phi_from, a_from, random_rotations
from stress import random_stress_field


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
    ids[2, 1, 1] = -1; ids[2, 2, 1] = 1  # edge contact only: no exchange coupling
    _check_no_contact(ids)


def test_rotations():
    R = random_rotations(20, seed=3)
    for r in R:
        assert np.allclose(r @ r.T, np.eye(3), atol=1e-12)
        assert abs(np.linalg.det(r) - 1.0) < 1e-12


def test_stress_statistics_and_correlation():
    N, dx, xi, s = 64, 5e-9, 20e-9, 100e6
    sig = random_stress_field(N, dx, xi, s, seed=7)
    assert sig.shape == (N, N, N, 6)
    assert np.allclose(sig.reshape(-1, 6).std(axis=0), s, rtol=1e-6)
    assert np.allclose(sig.reshape(-1, 6).mean(axis=0), 0.0, atol=1e-6 * s)
    # autocorrelation along x at lag r = xi (4 cells): expect exp(-1/2) = 0.61
    f = sig[..., 0]
    lag = int(round(xi / dx))
    c = np.mean(f * np.roll(f, lag, axis=0)) / np.mean(f * f)
    assert abs(c - math.exp(-0.5)) < 0.1


def test_stress_masked_normalisation():
    N, dx = 32, 10e-9
    mask = np.zeros((N, N, N), bool); mask[:16] = True
    sig = random_stress_field(N, dx, 40e-9, 50e6, seed=1, mask=mask)
    assert np.allclose(sig[mask].std(axis=0), 50e6, rtol=1e-6)


def test_stress_field_term_is_energy_gradient():
    import torch
    from magnumnp import Mesh, State, constants
    from field_terms_extra import StressAnisotropyField
    torch.set_default_dtype(torch.float64)

    n = (4, 3, 2)
    mesh = Mesh(n, (2e-9, 2e-9, 2e-9))
    state = State(mesh)
    Ms, lam = 1.4e6, 5e-6
    state.material = {"Ms": Ms}
    rng = np.random.default_rng(0)
    sig = torch.tensor(rng.normal(0, 1e8, n + (6,)))
    m = torch.tensor(rng.normal(size=n + (3,)))
    m = m / torch.linalg.norm(m, dim=-1, keepdim=True)
    state.m = m.clone()

    term = StressAnisotropyField(sig, lam)
    h = term.h(state)

    # energy density e = -(3/2) lam sum_ij sigma_ij m_i m_j (Voigt: off-diagonals twice)
    def e_density(mm):
        sx, sy, sz, syz, sxz, sxy = sig.unbind(-1)
        mx, my, mz = mm.unbind(-1)
        q = sx*mx*mx + sy*my*my + sz*mz*mz + 2*(syz*my*mz + sxz*mx*mz + sxy*mx*my)
        return -1.5 * lam * q

    mm = m.clone().requires_grad_(True)
    e_density(mm).sum().backward()
    h_ref = -mm.grad / (constants.mu_0 * Ms)
    assert torch.allclose(h, h_ref, rtol=1e-10, atol=1e-6)

    E_ref = float(e_density(m).sum()) * mesh.cell_volumes
    assert abs(float(term.E(state)) - E_ref) < 1e-10 * abs(E_ref)


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
