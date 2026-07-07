import pytest
import torch
from magnumnp import *
from magnumnp.common import accumulate_h


def make_state():
    mesh = Mesh((8, 6, 4), (2e-9, 2e-9, 2e-9))
    state = State(mesh)
    state.material = {"Ms": 8e5, "alpha": 0.1}
    state.m = state.Constant([0., 0., 1.])
    return state


def test_thermal_prefactor_invalidation():
    state = make_state()
    state.T = 300.
    state._dt = 1e-13
    state._step = 1

    thermal = ThermalField()
    h0 = thermal.h(state)

    # in-place domain write keeps the tensor identity -> version counter must invalidate
    x, y, z = state.mesh.SpatialCoordinate()
    domain = x < 4e-9
    state.material["alpha"][domain] = 0.4
    state._step = 2
    h1 = thermal.h(state)

    ratio = (h1[0, 0, 0] / (thermal._sigma[0, 0, 0] + 1e-300)).abs().mean()
    pref = thermal._pref
    assert pref[0, 0, 0, 0] != pref[-1, 0, 0, 0] # prefactor reflects the domain write

    # full reassignment (new object) also invalidates
    state.material["alpha"] = 0.2
    state._step = 3
    thermal.h(state)
    assert torch.allclose(thermal._pref[0, 0, 0], thermal._pref[-1, 0, 0])


def test_rux_coefficient_invalidation():
    state = make_state()
    mat = state.Constant([1], dtype=int)
    J = torch.zeros((3, 2), dtype=torch.get_default_dtype())
    J[2, :] = torch.tensor([7e-21, 3e-21])
    state.material["RuxDistribution"] = mat

    rux = AtomisticRuXExchangeField(J)
    h0 = rux.h(state)

    state.material["Ms"] = 4e5 # halved Ms -> doubled field
    h1 = rux.h(state)
    assert torch.allclose(h1, 2. * h0)


def test_external_cache_write_isolation():
    state = make_state()
    external = ExternalField([1e5, 0., 0.])
    h0 = external.h(state)
    h1 = external.h(state)
    assert h0.data_ptr() == h1.data_ptr() # cached storage is shared between calls

    # an indexed write by the caller materializes a private copy...
    h0[0, 0, 0, 0] = 0.
    assert h0[0, 0, 0, 0] == 0.
    # ...and must not corrupt the cached field
    torch.testing.assert_close(state.avg(external.h(state)), torch.tensor([1e5, 0., 0.]))


def test_material_decoupled_from_source_buffer():
    state = make_state()
    axis = torch.tensor([0., 0., 1.])
    state.material["Ku_axis"] = axis
    axis[2] = -1. # user mutates their own buffer afterwards
    torch.testing.assert_close(state.avg(state.material["Ku_axis"]), torch.tensor([0., 0., 1.]))


def test_accumulate_h_single_term():
    state = make_state()
    external = ExternalField([1e5, 0., 0.])
    total = accumulate_h(term.h(state) for term in [external])
    assert total is not external.h(state) # single term: caller gets a private tensor
    total[0, 0, 0, 0] = 0. # and may modify it freely
    torch.testing.assert_close(state.avg(external.h(state)), torch.tensor([1e5, 0., 0.]))
