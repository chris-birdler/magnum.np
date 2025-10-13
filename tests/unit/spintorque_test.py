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
        "Ms": state.Constant(8e5),
        "A": state.Constant(1.3e-11),
        "alpha": state.Constant(0.01),
        "P": state.Constant(0.5669),
        "Lambda": state.Constant(2),
        "epsilon_prime": state.Constant(1),
        "mp": state.Constant([0.93, 0.34, 0]),
        "d": state.Constant(5e-9),
        "J": state.Constant(-4e11),
        }

    torque = SpinTorqueSlonczewski()

    state.material["mp"] = state.Constant([0, -1, 0])
    h = torque.h(state)
