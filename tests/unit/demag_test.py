import pytest
import pathlib
import torch
from magnumnp import *
from magnumnp.field_terms.demag import newell, f, g
import numpy as np


@pytest.mark.parametrize("dtype", [torch.float32, torch.float64])
def test_energy_cube(dtype):
    torch.set_default_dtype(dtype)
    n  = (20, 20, 20)
    dx = (1e-9, 1e-9, 1e-9)
    mesh_volume = n[0] * n[1] * n[2] * dx[0] * dx[1] * dx[2]
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": state.Constant(Ms)}
    demag = DemagField()

    state.m = state.Constant([1,0,0])

    h = demag.h(state)
    torch.testing.assert_close(constants.mu_0*h[10,10,10], torch.tensor([-1./3., 0.0, 0.0]), atol=1e-3, rtol=1e-3)

    E = demag.E(state)
    assert E.cpu() == pytest.approx(1./6.*mesh_volume*constants.mu_0*Ms**2, abs=0, rel=1e-3)

def test_single():
    n  = (1, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {'Ms': state.Constant(1.)}
    demag = DemagField()

    state.m = state.Constant([1,0,0])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([-1./3.,0.,0.]), atol=1e-10, rtol=1e-10)

def test_truePBC():
    n  = (4, 4, 4)
    dx = (1e-9, 1e-9, 1e-9)
    Ms = 1./constants.mu_0
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {'Ms': state.Constant(Ms)}
    state.m = state.Constant([1,0,0])
    demag = DemagFieldPBC()

    h = demag.h(state)
    torch.testing.assert_close(h, state.Constant([0.0, 0.0, 0.0]), atol=1e-15, rtol=0)
    E = demag.E(state)
    assert E.cpu() == pytest.approx(0., abs=1e-15, rel=0.)

def test_pseudoPBC():
    n  = (1, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)

    # cube
    mesh = Mesh(n, dx, pbc = (0,0,0))
    state = State(mesh)
    state.material = {'Ms': state.Constant(1.)}
    demag = DemagField()
    state.m = state.Constant([1,0,0])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([-1./3.,0.,0.]), atol=1e-10, rtol=1e-10)
    state.m = state.Constant([0,1,0])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([0.,-1./3.,0.]), atol=1e-10, rtol=1e-10)
    state.m = state.Constant([0,0,1])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([0.,0.,-1./3.]), atol=1e-10, rtol=1e-10)

    # long cylinder
    mesh = Mesh(n, dx, pbc = (10,0,0))
    state = State(mesh)
    state.material = {'Ms': state.Constant(1.)}
    demag = DemagField()
    state.m = state.Constant([1,0,0])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([0.,0.,0.]), atol=1e-2, rtol=1e-2)
    state.m = state.Constant([0,1,0])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([0.,-1./2.,0.]), atol=1e-2, rtol=1e-2)
    state.m = state.Constant([0,0,1])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([0.,0.,-1./2.]), atol=1e-2, rtol=1e-2)

    # thin film
    mesh = Mesh(n, dx, pbc = (10,10,0))
    state = State(mesh)
    state.material = {'Ms': state.Constant(1.)}
    demag = DemagField()
    state.m = state.Constant([1,0,0])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([0.,0.,0.]), atol=5e-2, rtol=5e-2)
    state.m = state.Constant([0,1,0])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([0.,0.,0.]), atol=5e-2, rtol=5e-2)
    state.m = state.Constant([0,0,1])
    h = demag.h(state)
    torch.testing.assert_close(h[0,0,0,:], torch.tensor([0.,0.,-1.]), atol=5e-2, rtol=5e-2)

@pytest.mark.parametrize("nx", [5,6])
def test_pseudoPBC_single(nx):
    n  = (nx, 4, 3)
    dx = (1e-9, 2e-9, 5e-9)

    mesh = Mesh(n, dx, pbc = (1,0,0))
    state = State(mesh)
    state.material = {'Ms': state.Constant(1.)}
    demag = DemagField()
    state.m = state.Constant([0,0,0])
    state.m[0,0,0,0] = 1.

    h = demag.h(state)

    h2xa = newell(f, torch.tensor(2*dx[0]), torch.tensor(2*dx[1]), torch.tensor(0), *dx, *dx)
    h2xb = newell(f, torch.tensor((2+n[0])*dx[0]), torch.tensor(2*dx[1]), torch.tensor(0), *dx, *dx)
    h2xc = newell(f, torch.tensor((2-n[0])*dx[0]), torch.tensor(2*dx[1]), torch.tensor(0), *dx, *dx)
    h2x = h2xa + h2xb + h2xc
    torch.testing.assert_close(h[2,2,0,0]/h.max(), h2x/h.max(), atol=1e-10, rtol=1e-10)

    h2ya = newell(g, torch.tensor(2*dx[0]), torch.tensor(2*dx[1]), torch.tensor(0), *dx, *dx)
    h2yb = newell(g, torch.tensor((2+n[0])*dx[0]), torch.tensor(2*dx[1]), torch.tensor(0), *dx, *dx)
    h2yc = newell(g, torch.tensor((2-n[0])*dx[0]), torch.tensor(2*dx[1]), torch.tensor(0), *dx, *dx)
    h2y = h2ya + h2yb + h2yc

    torch.testing.assert_close(h[2,2,0,1]/h.max(), h2y/h.max(), atol=1e-10, rtol=1e-10)

def test_pseudoPBC_vs_repeat():
    n  = (108, 1, 50)
    dx = (1e-9, 16e-9, 1e-9)

    # pbc
    mesh = Mesh(n, dx, pbc = (1,0,0))
    state = State(mesh)
    state.material = {'Ms': state.Constant(2.1/constants.mu_0)}

    demag = DemagField()
    state.m = state.Constant([0.0, 1.0, 0.0])
    state.m[:20,:,:,:] = 0
    state.m[-20:,:,:,:] = 0
    h = demag.h(state)

    # nopbc (repeat)
    mesh2 = Mesh((3*n[0], n[1], n[2]), dx, origin=(-n[0]*dx[0],0.,0.))
    state2 = State(mesh2)
    state2.material = {'Ms': state2.Constant(2.1/constants.mu_0)}
    demag2 = DemagField()
    state2.m = torch.vstack([state.m,state.m,state.m])
    write_vti(state.m, "data/m.vti", state)
    write_vti(state2.m, "data/m2.vti", state2)
    h2 = demag2.h(state2)

    torch.testing.assert_close(h/h.abs().max(), h2[n[0]:2*n[0],:,:,:]/h.abs().max(), atol=1e-3, rtol=1e-3)

def test_nonequi_vs_equi():
    n  = (4, 4, 4)
    dx = (1., 2., 5.)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {'Ms': state.Constant(1.)}

    x, y, z = mesh.SpatialCoordinate()
    state.m = Expression([x*y, y*z, z*x])

    demag = DemagField()
    h1 = demag.h(state)

    demag2 = DemagFieldNonEquidistant()
    h2 = demag2.h(state)

    torch.testing.assert_close(h1, h2, atol=1e-10, rtol=1e-10)

def test_nonequi_vs_equi_far_layers():
    # layer pairs beyond the Newell radius (p = 20 cell diagonals): the
    # far-field on-axis coupling must be kept. Previous versions zeroed the
    # local [0,0,0] entry of the dipole formula positionally, which dropped
    # the dominant interlayer coupling for distant layers; the field of a
    # single magnetized layer then vanished in layers further away than the
    # Newell radius.
    n  = (4, 4, 50)
    dx = (1., 1., 1.)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {'Ms': state.Constant(1.)}
    state.m = state.Constant([0., 0., 0.])
    state.m[:,:,0,2] = 1. # magnetize only the bottom layer

    h1 = DemagField().h(state)
    h2 = DemagFieldNonEquidistant().h(state)

    # the far layers carry a small but essential dipole field
    assert h1[:,:,-1,:].abs().max() > 0.
    torch.testing.assert_close(h1, h2, atol=1e-10, rtol=1e-6)

def test_nonequidistant():
    n  = (20, 4, 6)
    dx1 = (5., 3., 1.)
    mesh1 = Mesh(n, dx1)
    state1 = State(mesh1)
    state1.material = {'Ms': state1.Constant(1.)}

    state1.m = state1.Constant([0.,0.,0.])
    state1.m[0,0,0:2,0] = 1.

    demag1 = DemagField()
    h1 = demag1.h(state1)

    n  = (20, 4, 4)
    dx2 = np.ones(n[2]) * 2.
    dx2[2:] = 1.
    dx = (5., 3., dx2)
    mesh2 = Mesh(n, dx)
    state2 = State(mesh2)
    state2.material = {'Ms': state2.Constant(1.)}

    state2.m = state2.Constant([0.,0.,0.])
    state2.m[0,0,0,0] = 1.

    demag2 = DemagFieldNonEquidistant()
    h2 = demag2.h(state2)

    torch.testing.assert_close(h1[:,:,-2:,:], h2[:,:,-2:,:], atol=1e-10, rtol=1e-10)

def test_regression():
    n  = (101, 101, 101)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"Ms": state.Constant(1.)}
    state.m = state.Constant([0,0,0])
    state.m[50,50,50] = 1.

    demag = DemagField()
    h_demag = demag.h(state)

    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "ref" / "h_demag_regression.vti"
    ## Uncomment to updated reference data
    #write_vti({"h_demag":h_demag}, filename)
    mesh, ref = read_vti(filename)

    torch.testing.assert_close(h_demag, ref["h_demag"], atol=1e-10, rtol=1e-6)

def test_precision():
    ''' if Newell formula is evaluated with single precision large fluctuations are expected '''
    n = (100,1,1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    dtype = torch.get_default_dtype()
    torch.set_default_dtype(torch.float32)
    state = State(mesh)
    state.material = {"Ms": state.Constant(1.)}
    state.m = state.Constant([0,0,0])
    state.m[0,0,0,2] = 1.

    demag = DemagField()
    h_demag = demag.h(state).abs().log()[:,0,0,2]

    ref = torch.arange(n[0])
    ref = torch.log(1/4./torch.pi/ref**3)

    torch.set_default_dtype(dtype)
    torch.testing.assert_close(h_demag[20:], ref[20:], atol=1e-2, rtol=1e-2)


def test_block_demag_factors():
    from magnumnp.field_terms.demag import block_demag_factors
    torch.set_default_dtype(torch.float64)
    N = block_demag_factors(State(Mesh((4, 4, 4), (1e-9, 1e-9, 1e-9), pbc=(2, 2, 2))))
    assert N == pytest.approx([1./3.]*3, abs=1e-14)
    N = block_demag_factors(State(Mesh((4, 2, 8), (1e-9, 2e-9, 3e-9), pbc=(1, 3, 2))))
    assert sum(N) == pytest.approx(1., abs=1e-14)
    N = block_demag_factors(State(Mesh((40, 1, 1), (1e-9, 1e-9, 1e-9), pbc=(1, 1, 1))))   # long along x
    assert N[0] < 0.001 and N[1] == pytest.approx(N[2])
    with pytest.raises(ValueError):
        block_demag_factors(State(Mesh((4, 4, 4), (1e-9, 1e-9, 1e-9), pbc=(1, 1, 0))))


@pytest.mark.parametrize("u", [[1., 0., 0.], [1., 1., 0.], [1., 1., 1.]])
def test_true_pbc_newell_aligned_cube(u):
    """Uniform m in an aligned cube of edge L/2 in a periodic cube: the exact energy is f (1 - f) / 3 * mu0 Ms^2 / 2 * V
    with the filling f = 1/8, for every direction of m (cubic symmetry, zero mean field). DemagField(true_pbc=True)
    gives it; DemagFieldPBC (finite differences) gives it only for m along an axis (+6.9 % along (110), +13.7 % along
    (111) for 4 cells per cube edge)."""
    torch.set_default_dtype(torch.float64)
    n, dx, Ms = (8, 8, 8), (1e-9, 1e-9, 1e-9), 1./constants.mu_0
    Kd_V = 0.5 * constants.mu_0 * Ms**2 * 512 * 1e-27
    ms = torch.zeros(n + (1,))
    ms[:4, :4, :4] = Ms
    state = State(Mesh(n, dx, pbc=(1, 1, 1)))
    state.material = {"Ms": ms}
    state.m = (torch.tensor(u) / torch.tensor(u).norm()).expand(n + (3,)).clone()
    exact = (1/8) * (7/8) / 3
    assert float(DemagField(true_pbc=True).E(state)) / Kd_V == pytest.approx(exact, rel=1e-3)
    e_fd = float(DemagFieldPBC().E(state)) / Kd_V
    if u.count(0.) == 2:
        assert e_fd == pytest.approx(exact, rel=1e-9)
    else:
        assert e_fd / exact - 1 > 0.05


def test_true_pbc_newell_checkerboard_and_zero_mean():
    """A cell-scale pattern m = (1, -1, 0)/sqrt2 (-1)^(i+j) in a full box: DemagFieldPBC sees no charge (energy 0);
    the exact energy (Fourier sum over the aliases of the cell form factor) is 0.50 mu0 Ms^2 / 2 per volume. For a
    random m the box mean of h is 0 (zero mean field of the infinite medium) within the image truncation; without
    the compensation it is the field of the outer surface of the image block."""
    torch.set_default_dtype(torch.float64)
    n, dx, Ms = (8, 8, 8), (1e-9, 1e-9, 1e-9), 1./constants.mu_0
    Kd_V = 0.5 * constants.mu_0 * Ms**2 * 512 * 1e-27
    i, j, k = torch.meshgrid(*[torch.arange(8)] * 3, indexing="ij")
    s = (-1.)**(i + j)
    state = State(Mesh(n, dx, pbc=(1, 1, 1)))
    state.material = {"Ms": state.Constant(Ms)}
    state.m = torch.stack([s, -s, 0. * s], -1) / 2**0.5
    assert float(DemagField(true_pbc=True).E(state)) / Kd_V == pytest.approx(0.50, rel=0.01)
    assert abs(float(DemagFieldPBC().E(state)) / Kd_V) < 1e-12
    g = torch.Generator().manual_seed(0)
    r = torch.randn(n + (3,), generator=g)
    state.m = r / r.norm(dim=-1, keepdim=True)
    assert float(DemagField(true_pbc=True).h(state).mean(dim=(0, 1, 2)).abs().max()) < 1e-3 * Ms
    assert float(DemagField().h(state).mean(dim=(0, 1, 2)).abs().max()) > 3e-3 * Ms


def test_true_pbc_extrapolated_kernel():
    """DemagField(true_pbc=True, extrapolate=True) with k shells: the kernel N_k + (N_k - N_i) x_k / (x_i - x_k),
    x_j = 1/(2j+1)^2, k_i = k // 2. The energy is linear in the kernel, so E_ext = E_k + (E_k - E_i) x_k / (x_i - x_k).
    For two domains (zero mean, large first moment) the image sum converges as 1/(2k+1)^2: the extrapolated
    energies with 2 and 4 shells agree much better than the plain ones (PROTOKOLL 7.108). p = 2 keeps the test fast."""
    torch.set_default_dtype(torch.float64)
    n, dx, Ms = (4, 4, 4), (1e-9, 1e-9, 1e-9), 1./constants.mu_0
    m = torch.zeros(n + (3,))
    m[:2, ..., 2] = 1.
    m[2:, ..., 2] = -1.

    def energy(k, ext=False):
        state = State(Mesh(n, dx, pbc=(k, k, k)))
        state.material = {"Ms": state.Constant(Ms)}
        state.m = m.clone()
        return float(DemagField(p=2, true_pbc=True, extrapolate=ext).E(state))

    E = {k: energy(k) for k in (1, 2, 4)}
    X = {k: energy(k, True) for k in (2, 4)}
    assert X[4] == pytest.approx(E[4] + (E[4] - E[2]) * (1/81) / (1/25 - 1/81), rel=1e-10)
    assert X[2] == pytest.approx(E[2] + (E[2] - E[1]) * (1/25) / (1/9 - 1/25), rel=1e-10)
    assert abs(E[4] - E[2]) > 1e-3 * abs(E[4])
    assert abs(X[4] - X[2]) < 0.1 * abs(E[4] - E[2])
    with pytest.raises(ValueError):
        DemagField(extrapolate=True)
    state = State(Mesh(n, dx, pbc=(1, 1, 1)))
    state.material = {"Ms": state.Constant(Ms)}
    state.m = m.clone()
    with pytest.raises(ValueError):
        DemagField(true_pbc=True, extrapolate=True).h(state)


def _h_reference(dem, state):
    """h() of DemagField before the fused kernel product (audit item A2): loop over the components, the block term
    of true_pbc as a separate term."""
    dim = [i for i in range(3) if state.mesh.n[i] > 1]
    s = [dem._shape(state)[i] for i in dim]
    M = state.material["Ms"] * state.m
    hk = [0, 0, 0]
    for ax in range(3):
        mk = torch.fft.rfftn(M[..., ax], dim=dim, s=s)
        for i in range(3):
            hk[i] = hk[i] + dem._N[i][ax] * mk
    h = torch.stack([torch.fft.irfftn(x, s=s, dim=dim)[:state.mesh.n[0], :state.mesh.n[1], :state.mesh.n[2]]
                     for x in hk], dim=3)
    if dem._true_pbc:
        from magnumnp.field_terms.demag import block_demag_factors
        h = h + torch.tensor(block_demag_factors(state), dtype=h.dtype) * M.mean(dim=(0, 1, 2))
    return h


@pytest.mark.parametrize("dtype, tol", [(torch.float64, 1e-12), (torch.float32, 1e-5)])
@pytest.mark.parametrize("pbc, kw", [((2, 2, 2), {"true_pbc": True, "extrapolate": True}),
                                     ((1, 1, 1), {"true_pbc": True}),
                                     ((1, 1, 1), {}),
                                     ((0, 0, 0), {}),
                                     ((2, 0, 1), {})])
def test_fused_h_against_reference(dtype, tol, pbc, kw):
    """The fused kernel product (1/n and the block term of true_pbc in the kernel) gives the field of the loop over
    the components within round-off (PLAN 11, S1): random m, a varying Ms with a void, periodic and open."""
    default = torch.get_default_dtype()
    torch.set_default_dtype(dtype)
    try:
        _check_fused_h(dtype, tol, pbc, kw)
    finally:
        torch.set_default_dtype(default)


def _check_fused_h(dtype, tol, pbc, kw):
    n, dx = (8, 6, 10), (1e-9, 1.5e-9, 1e-9)
    g = torch.Generator().manual_seed(1)
    state = State(Mesh(n, dx, pbc=pbc))
    ms = 1e6 * (0.5 + torch.rand(n + (1,), generator=g, dtype=torch.float64)).to(dtype)
    ms[:2, :2, :2] = 0.
    state.material = {"Ms": ms}
    r = torch.randn(n + (3,), generator=g, dtype=torch.float64).to(dtype)
    state.m = r / r.norm(dim=-1, keepdim=True)
    dem = DemagField(p=2, **kw)
    h = dem.h(state)
    h_ref = _h_reference(dem, state)
    assert float((h - h_ref).abs().max() / h_ref.abs().max()) < tol
    E, E_ref = float(dem.E(state)), float((-0.5 * constants.mu_0 * ms * state.m * h_ref * state.mesh.cell_volumes).sum())
    assert abs(E / E_ref - 1) < tol
