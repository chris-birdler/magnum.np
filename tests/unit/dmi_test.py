import pytest
import pathlib
import torch
from magnumnp import *

# move to fieldterms test
def test_call():
    n = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.m = state.Constant([1,0,0])
    state.material = {"Ms": 1./constants.mu_0,
                      "xDi": 1.,
                      "xDb": 1.,
                      "xDD2d": 1.}

    dmi = InterfaceDMIField(Di = "xDi")
    dmi.h(state)
    dmi = BulkDMIField(Db = "xDb")
    dmi.h(state)
    dmi = D2dDMIField(DD2d = "xDD2d")
    dmi.h(state)


# move to demo as integration test
def test_interface_1D():
    """ 1D Interface DMI Testcase
        taken from 'Proposal for a micromagnetic standard problem for materials with Dzyaloshinskii–Moriya interaction'
    """
    n = (20, 1, 1) # (100, 1, 1) for discretization as in Paper
    dx = (5e-9, 1e-9, 1e-9)
    origin = (-n[0]*dx[0]/2., -n[1]*dx[1]/2., -n[2]*dx[2]/2.,)
    mesh = Mesh(n, dx, origin)
    state = State(mesh)
    state.material = {"Ms": 0.86e6,
                      "A": 13e-12,
                      "Ku": 0.4e6,
                      "Ku_axis": [0,0,1],
                      "Di": -3e-3,
                      "alpha": 1.}
    state.m = state.Constant([0.1,0.0,1.0])
    normalize(state.m)

    exchange = ExchangeField()
    aniso = UniaxialAnisotropyField()
    dmi = InterfaceDMIField()

    llg = LLGSolver([exchange, aniso, dmi])
    logger = ScalarLogger("data/m_relax.dat", ['t', 'm'])
    while state.t <= 2e-10:
        logger << state
        llg.step(state, 1e-11)
    assert state.m[...,0].max().cpu() == pytest.approx(0.443, rel=1e-2)
    assert state.m[...,2].min().cpu() == pytest.approx(0.896, rel=1e-2)

@pytest.mark.parametrize("dmi_term", [InterfaceDMIField, BulkDMIField, D2dDMIField])
def test_nonequi_vs_equi(dmi_term):
    n  = (10, 10, 10)
    dx1 = (1e-9, 2e-9, 5e-9)
    mesh1 = Mesh(n, dx1)
    state1 = State(mesh1)
    state1.material = {"Ms": 1./constants.mu_0,
                       "Di": 1.,
                       "Db": 1.,
                       "DD2d": 1.}

    dx2 = (torch.ones(n[0]) * 1e-9, 2e-9, 5e-9)
    mesh2 = Mesh(n, dx2)
    state2 = State(mesh2)
    state2.material = {"Ms": 1./constants.mu_0,
                       "Di": 1.,
                       "Db": 1.,
                       "DD2d": 1.}
    x, y, z = state1.SpatialCoordinate()
    state1.m = torch.stack([x*y, y*z, z*x], dim=-1)
    state2.m = torch.stack([x*y, y*z, z*x], dim=-1)

    dmi = dmi_term()
    h1 = dmi.h(state1).cpu()
    h2 = dmi.h(state2).cpu()
    torch.testing.assert_close(h1/h1.max(), h2/h1.max(), atol=1e-10, rtol=0)
