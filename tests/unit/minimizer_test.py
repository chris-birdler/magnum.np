import pytest
import torch
from magnumnp import *

def test_sp4():
    # initialize mesh
    eps = 1e-15
    n  = (100, 25, 1)
    dx = (5e-9, 5e-9, 3e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)

    state.material = {
        "Ms": state.Constant([8e5]),
        "A": state.Constant([1.3e-11]),
        "alpha": state.Constant([0.02])
        }

    # initialize field terms
    demag    = DemagField()
    exchange = ExchangeField()
    external = ExternalField([-24.6e-3/constants.mu_0,
                              +4.3e-3/constants.mu_0,
                              0.0])

    # initialize magnetization that relaxes into s-state
    state.m = state.Constant([0,0,0])
    state.m[1:-1,:,:,0]   = 1.0
    state.m[(-1,0),:,:,1] = 1.0

    # relax without external field
    minimizer = MinimizerBB([demag, exchange])
    minimizer.minimize(state)

    torch.testing.assert_close(avg(state.m), state.Tensor([0.96720725, 0.12482232, 0.]), atol=1e-3, rtol=1e-3)
