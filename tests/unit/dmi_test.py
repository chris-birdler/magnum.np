import pytest
import pathlib
import torch
from magnumnp import *

def test_call():
    n  = (8,10,12)
    dx = (1e-9, 2e-9, 5e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.m = state.Constant([1,0,0])
    state.material = {"Ms": Ms,
                      "Di": 1.,
                      "Db": 1.,
                      "DD2d": 1.}

    dmi = InterfaceDMIField()
    dmi.h(state)
    dmi = BulkDMIField()
    dmi.h(state)
    dmi = D2dDMIField()
    dmi.h(state)

