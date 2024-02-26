import pytest
import torch
from magnumnp import *
from helpers import *

def test_zhangli(simple_state):
    state = simple_state
    state.m = state.Constant([0,1,0])
    state.m[:1,:,:,1] = -1.

    torque = SpinTorqueZhangLi()

    state.j = state.Constant([1e12,0,0])
    h = torque.h(state)

def test_sot(simple_state):
    state = simple_state
    state.m = state.Constant([0,1,0])
    state.m[:1,:,:,1] = -1.

    torque = SpinOrbitTorque()

    state.material["p"] = state.Constant([0, -1, 0])
    h = torque.h(state)

def test_slonczewski(simple_state):
    state = simple_state
    state.m = state.Constant([0,1,0])
    state.m[:1,:,:,1] = -1.

    state.material = {
        "Ms": 8e5,
        "A": 1.3e-11,
        "alpha": 0.01,
        "P": 0.5669,
        "Lambda": 2,
        "epsilon_prime": 1,
        "mp": state.Tensor((0.93, 0.34, 0)),
        "d": 5e-9,
        "J": -4e11,
        }

    torque = SpinTorqueSlonczewski()

    state.material["mp"] = state.Tensor([0, -1, 0])
    h1 = torque.h(state)
    state.material["mp"] = state.Constant([0, -1, 0])
    h2 = torque.h(state)

    assert torch.allclose(h1, h2)
