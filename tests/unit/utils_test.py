import pytest
import torch
from magnumnp import *

def test_MFM():
    n = (100, 20, 10)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2., -n[1]*dx[1]/2., -n[2]*dx[2]/2.))
    state = State(mesh)
    state.material = {"Ms": 8e5}
    x, y, z = state.SpatialCoordinate()

    state.m = torch.stack([y, -x, 0*z], dim=-1)
    state.m.normalize()

    mfm = MFM(height=10e-9, mm_tip = 10e-9, dm_tip=20e-9)
    phi = mfm.PhaseShift(state)

    print ("phi:", phi.shape)
