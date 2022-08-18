import pytest
import torch
from math import pi, cos, sin
from magnumnp import *

def test_cubic_anisotropy():
    n  = (1,1,1)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx)
    state = State(mesh, {})
    state._material = {"Kc_alpha": state.Constant([0.]),
                       "Kc_beta": state.Constant([0.]),
                       "Kc_gamma": state.Constant([0.]),
                       "Kc1": state.tensor([1e3]),
                       "Kc2": state.tensor([0.]),
                       "Ms": 800e3}
    aniso = CubicAnisotropyField(state)
    for phi in torch.linspace(0., pi, steps=10):
        mx = cos(phi)
        my = sin(phi)
        mz = 0.
        state.m = state.Constant((mx, my, mz))

        h_sim = aniso.h(0, state.m)
        E_sim = aniso.E(0, state.m)

        hx_analytic = -2. / constants.mu_0 / state._material["Ms"] * state._material["Kc1"] * mx * (my**2 + mz**2)
        E_analytic = state._material["Kc1"] * (mx**2 * my**2 + mx**2 * mz**2 + my**2 * mz**2) * mesh.cell_volume
        assert h_sim[0,0,0,0] == pytest.approx(hx_analytic[0].cpu().numpy())
        assert E_sim == pytest.approx(E_analytic.cpu().numpy())


    state._material["Kc_alpha"] = state.Constant([pi/4.])
    aniso = CubicAnisotropyField(state)
    for phi in torch.linspace(0., pi, steps=10):
        mx = cos(phi - pi/4)
        my = sin(phi - pi/4)
        mz = 0.
        state.m = state.Constant((mx, my, mz))
        h_sim = aniso.h(0, state.m)
        E_sim = aniso.E(0, state.m)

        mx = cos(phi)
        my = sin(phi)
        mz = 0.
        hx_analytic = -2. / constants.mu_0 / state._material["Ms"] * state._material["Kc1"] * mx * (my**2 + mz**2)
        E_analytic = state._material["Kc1"] * (mx**2 * my**2 + mx**2 * mz**2 + my**2 * mz**2) * mesh.cell_volume

#        assert h_sim[0,0,0,0] == pytest.approx(hx_analytic[0].cpu().numpy())
        assert E_sim == pytest.approx(E_analytic.cpu().numpy())
