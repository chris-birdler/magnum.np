import pytest
import torch
from magnumnp import *
from helpers import *

def finite_grad(op, x0):
    eps = 1. #e-5
    y0 = op(x0)
    grad = torch.zeros_like(x0)
    for i in range(int(x0.shape[0])):
        x = x0.detach().clone()
        x[i] += eps
        y2 = op(x)
        x[i] -= 2*eps
        y1 = op(x)
        grad[i] = (y2-y1)/(2*eps)
    return grad

def test_hext(simple_state):
    h_ext = torch.tensor([1.,0.,0.], requires_grad = True)
    external = ExternalField(h_ext)
    h_target = torch.tensor([0.,0.,1.])

    def forward(h_ext):
        external.h = h_ext
        L = ((external.h(simple_state)-h_target)**2).sum()
        return L

    L = forward(h_ext)
    L.backward()
    grad = h_ext.grad

    h_ext2 = torch.tensor([1.,0.,0.])
    grad2 = finite_grad(forward, h_ext2)

    torch.testing.assert_close(grad, grad2)

@pytest.mark.parametrize("fieldterm", [UniaxialAnisotropyField(),
                                       ExchangeField(),
                                       DemagField(),
                                       InterfaceDMIField(),
                                       BulkDMIField(),
                                       D2dDMIField()])
def test_linear_fieldterms(simple_state, fieldterm):
    cell_volume = simple_state.mesh.dx[0] * simple_state.mesh.dx[1] * simple_state.mesh.dx[2]
    x,y,z = simple_state.SpatialCoordinate()
    m0 = torch.stack([x,2*x,y], dim=-1)
    m0.requires_grad = True

    simple_state.m = m0
    h0 = fieldterm.h(simple_state)
    def forward(m):
        simple_state.m = m
        return fieldterm.E(simple_state)

    L = forward(m0)
    L.backward()
    h1 = -m0.grad / (constants.mu_0 * simple_state.material["Ms"] * cell_volume)

    write_vti({"h1":h1, "h0":h0}, "data/results.vti", simple_state)
    torch.testing.assert_close(h0/h0.abs().max(), h1/h0.abs().max(), atol=1e-4, rtol=1e-4)
