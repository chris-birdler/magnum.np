import torch
import pytest
from magnumnp import *


def test_energy_sums_in_float64():
    """The energies of the field terms sum the cell energies in float64, also for a float32 state (PROTOKOLL 7.110:
    the relaxation test compares energy changes of about 1e-6 relative, the size of the rounding of a float32 sum over
    millions of cells)."""
    old = torch.get_default_dtype()
    torch.set_default_dtype(torch.float32)
    try:
        n, dx = (16, 16, 16), (2e-9, 2e-9, 2e-9)
        state = State(Mesh(n, dx, pbc=(1, 1, 1)))
        state.material = {"Ms": state.Constant(8e5), "A": state.Constant(1e-11), "Ku": state.Constant(1e4),
                          "Ku_axis": state.Constant([0., 0., 1.])}
        g = torch.Generator().manual_seed(1)
        r = torch.randn(n + (3,), generator=g).to(state.device)
        state.m = r / r.norm(dim=-1, keepdim=True)
        assert state.m.dtype == torch.float32
        for term in (ExchangeField(), UniaxialAnisotropyField(), DemagFieldPBC()):
            E = term.E(state)
            assert E.dtype == torch.float64
            h = term.h(state)
            ref = (-0.5 * constants.mu_0 * state.material["Ms"] * state.m * h * state.mesh.cell_volumes).double().sum()
            assert float(E) == pytest.approx(float(ref), rel=1e-6)
    finally:
        torch.set_default_dtype(old)
