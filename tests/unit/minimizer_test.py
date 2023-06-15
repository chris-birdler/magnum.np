import pytest
import torch
from magnumnp import *

def test_descent():
    n  = (100, 25, 1)
    dx = (5e-9, 5e-9, 3e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {
        "Ms": 8e5,
        "A": 1.3e-11,
        }
    demag    = DemagField()
    exchange = ExchangeField()
    
    m0 = state.Constant([0,0,0])
    m0[1:-1,:,:,0]   = 1.0
    m0[(-1,0),:,:,1] = 1.0
    
    minimizer = Minimizer_BB([demag, exchange])

    def finite_grad(state):
        eps = 1e-5
        y0 = minimizer._E(state)
        m0 = state.m.clone().reshape(-1)
        grad = torch.zeros_like(m0)
        for i in range(int(m0.shape[0])):
            m = m0.clone()
            m[i] += eps
            state.m = m.reshape(state.m.shape)
            E2 = minimizer._E(state)
            m[i] -= 2*eps
            state.m = m.reshape(state.m.shape)
            E1 = minimizer._E(state)
            grad[i] = (E2-E1)/(2*eps)
        return grad.reshape(state.m.shape)

    def dm(state):
        h = sum([term.h(state) for term in minimizer._terms])
        return constants.mu_0*state.material["Ms"]*state.cell_volumes*torch.cross(state.m, torch.cross(state.m, h))

    state.m = m0.clone()
    grad1 = minimizer._dm(state)

    #state.m = m0.clone()
    #grad2 = finite_grad(state)

    state.m = m0.clone()
    grad3 = dm(state)

    print("grad1:", grad1.max().numpy())
#    print("grad2:", grad2.max().numpy())
    print("grad3:", grad3.max().numpy())
    write_vti(grad1, "data/grad1.vti", state)
#    write_vti(grad2, "data/grad2.vti", state)
    write_vti(grad3, "data/grad3.vti", state)
    #assert E.cpu() == pytest.approx(1./6.*mesh_volume*constants.mu_0*Ms**2, abs=0, rel=1e-3)

