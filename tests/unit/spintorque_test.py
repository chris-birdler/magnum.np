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

