import pytest
import torch
from math import pi, cos, sin
from magnumnp import *

def test_cubic_energy():
    n  = (1, 1, 1)
    dx = (1, 2, 5)
    cell_volume = dx[0] * dx[1] * dx[2]
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Kc_alpha": state.Constant([0.]),
                      "Kc_beta": state.Constant([0.]),
                      "Kc_gamma": state.Constant([0.]),
                      "Kc1": state.Tensor([1e3]),
                      "Kc2": state.Tensor([0.]),
                      "Ms": 800e3}
    aniso = CubicAnisotropyField()
    for phi in torch.linspace(0., pi, steps=10):
        mx = cos(phi)
        my = sin(phi)
        mz = 0.
        state.m = state.Constant((mx, my, mz))

        h_sim = aniso.h(state)
        E_sim = aniso.E(state)

        hx_analytic = -2. / constants.mu_0 / state.material["Ms"] * state.material["Kc1"] * mx * (my**2 + mz**2)
        E_analytic = (state.material["Kc1"] * (mx**2 * my**2 + mx**2 * mz**2 + my**2 * mz**2) * cell_volume).sum()

        torch.testing.assert_close(h_sim[0,0,0,0], hx_analytic[0,0,0,0], atol=1e-4, rtol=1e-4)
        torch.testing.assert_close(E_sim, E_analytic, atol=1e-4, rtol=1e-4)

def test_cubic_energy_rotated():
    n  = (1, 1, 1)
    dx = (1, 2, 5)
    cell_volume = dx[0] * dx[1] * dx[2]
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Kc_alpha": state.Constant([pi/4.]),
                      "Kc_beta": state.Constant([0.]),
                      "Kc_gamma": state.Constant([0.]),
                      "Kc1": state.Tensor([1e3]),
                      "Kc2": state.Tensor([0.]),
                      "Ms": 800e3}
    aniso = CubicAnisotropyField()
    for phi in torch.linspace(0., pi, steps=10):
        mx = cos(phi - pi/4)
        my = sin(phi - pi/4)
        mz = 0.
        state.m = state.Constant((mx, my, mz))
        h_sim = aniso.h(state)
        E_sim = aniso.E(state)

        mx = cos(phi)
        my = sin(phi)
        mz = 0.
        E_analytic = (state.material["Kc1"] * (mx**2 * my**2 + mx**2 * mz**2 + my**2 * mz**2) * cell_volume).sum()

        torch.testing.assert_close(E_sim, E_analytic, atol=1e-4, rtol=1e-4)

# TODO: generalize in Fieldterms test
def test_cubic_material_tensor():
    n  = (2, 3, 4)
    dx = (1, 2, 5)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Kc_alpha": state.Tensor([pi/4.]),
                      "Kc_beta": state.Tensor([0.]),
                      "Kc_gamma": state.Tensor([0.]),
                      "Kc1": state.Tensor([1e3]),
                      "Kc2": state.Tensor([0.]),
                      "Ms": 800e3}
    aniso = CubicAnisotropyField()
    phi = 0.123
    mx = cos(phi - pi/4)
    my = sin(phi - pi/4)
    mz = 0.
    state.m = state.Constant((mx, my, mz))
    torch.testing.assert_close(avg(aniso.h(state)), state.Tensor([-191.01109252, -148.98001006, 0.]), atol=1e-4, rtol=1e-4)

def test_renamed_parameters():
    n  = (2, 3, 4)
    dx = (1, 2, 5)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"K_alpha": state.Tensor([pi/4.]),
                      "Kc_beta": state.Tensor([0.]),
                      "Kc_gamma": state.Tensor([0.]),
                      "Kc1": state.Tensor([1e3]),
                      "Kc2": state.Tensor([0.]),
                      "Ms": 800e3}
    aniso = CubicAnisotropyField(Kc_alpha = "K_alpha")
    phi = 0.123
    mx = cos(phi - pi/4)
    my = sin(phi - pi/4)
    mz = 0.
    state.m = state.Constant((mx, my, mz))
    torch.testing.assert_close(avg(aniso.h(state)), state.Tensor([-191.01109252, -148.98001006, 0.]), atol=1e-4, rtol=1e-4)
