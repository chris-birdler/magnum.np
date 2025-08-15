import pytest
import torch
from magnumnp import *
import numpy as np

def test_C_mask():
    mask_iso = np.zeros((6,6), dtype=int)
    mask_iso[:3,:3] = 1
    for i in range(3):
        mask_iso[-i-1,-i-1] = 1
    
    llg_plain = LLGWithLESolver([])
    assert(np.all(llg_plain._C_mask == np.ones((6,6))))

    llg_str = LLGWithLESolver([], C_sym="isotropic")
    llg_np = LLGWithLESolver([], C_sym=mask_iso)
    llg_torch = LLGWithLESolver([], C_sym=torch.tensor(mask_iso))

    assert(np.all(llg_str._C_mask == mask_iso))
    assert(np.all(llg_np._C_mask == mask_iso))
    assert(np.all(llg_torch._C_mask == mask_iso))

    with pytest.raises(Exception):
        llg = LLGWithLESolver([], C_sym="invalid")
    with pytest.raises(Exception):
        llg = LLGWithLESolver([], C_sym=np.ones((3,4)))

def test_bcs():
    n = (10,11,12)
    dx0 = 3e-9
    dx = (dx0, dx0, dx0)

    mesh = Mesh(n, dx, pbc=(0,0,1))
    state = State(mesh)

    llg = LLGWithLESolver([])

    # Test 1:
    # no user set boundary values
    #############################

    llg._update_neumann_bcs(state)

    # Expected:
    # * Neumann BCs at all non-pbc boundaries
    # * Overlapp at the edges
    proof = torch.zeros(n, dtype=torch.int)
    proof[0] += 1
    proof[-1] += 1
    proof[:,0] += 1
    proof[:,-1] += 1

    for bc in llg._neumann_bcs:
        proof[bc.plane.slices] -= 1

    torch.testing.assert_close(proof.sum(), torch.tensor(0))
    torch.testing.assert_close(torch.max(proof), torch.tensor(0, dtype=torch.int))
    torch.testing.assert_close(torch.min(proof), torch.tensor(0, dtype=torch.int))

    # Test 2:
    # User set boundary values
    ##########################

    # Expected:
    # * Mostly as before
    # * Override natural BCs with user set values

    # with user set boundarie values
    # full x at 0
    t_xm = torch.zeros((n[1],n[2],3))
    t_xm[...,0] = torch.ones((n[1],n[2]))

    plane_xm = Plane(0,0,-1)
    bc_full_xm = PlaneBC(plane_xm, t_xm)

    # no conflic x interiour
    t_xi = torch.zeros((5,5,3))
    t_xi[...,1] = torch.ones((5,5))

    plane_xi = Plane(0,2,1,[1,6],[1,6])
    bc_xi = PlaneBC(plane_xi, t_xi)

    # conflicts on y at -1
    t_yp1 = torch.zeros((3,3,3))
    t_yp1[...,1] = torch.ones((3,3))

    plane_yp1 = Plane(1,-1,1,[1,4],[1,4])
    bc_yp1 = PlaneBC(plane_yp1, t_yp1)

    t_yp2 = torch.zeros((3,3,3))
    t_yp2[...,0] = -torch.ones((3,3))

    plane_yp2 = Plane(1,n[1]-1,1,[3,6],[3,6])
    bc_yp2 = PlaneBC(plane_yp2, t_yp2)

    t_bcs = [bc_full_xm, bc_xi, bc_yp1, bc_yp2]
    state.bcs = {"t" : t_bcs}

    llg._update_neumann_bcs(state)

    # test set values
    proof = torch.zeros(n+(3,))
    proof[0,:,:,0] = 1. # full x at 0
    proof[2,1:6,1:6,1] = 1. # x interiour
    proof[1:4,-1,1:4,1] = 1. # y conflict at -1
    proof[3:6,-1,3:6,0] = -1 # y conflict at -1

    for bc in llg._neumann_bcs:
        proof[bc.plane.slices] -= bc.condition(state)

    torch.testing.assert_close(proof.sum(), torch.tensor(0.))
    torch.testing.assert_close(torch.max(proof), torch.tensor(0.))
    torch.testing.assert_close(torch.min(proof), torch.tensor(0.))

    # test assignment
    proof = torch.zeros(n, dtype=torch.int)
    proof[0] += 1
    proof[-1] += 1
    proof[:,0] += 1
    proof[:,-1] += 1
    proof[2,1:6,1:6] += 1 # x interiour

    for bc in llg._neumann_bcs:
        proof[bc.plane.slices] -= 1

    torch.testing.assert_close(proof.sum(), torch.tensor(0))
    torch.testing.assert_close(torch.max(proof), torch.tensor(0, dtype=torch.int))
    torch.testing.assert_close(torch.min(proof), torch.tensor(0, dtype=torch.int))

def test_jump_condition():
    #######################
    # Temperature flow test
    #######################
    # Test jump conditions in first derivatives on the problem of heat flux
    # Problem setup: 
    #   1D heat flux. Two materials with coefficents CI and CII. 
    #   T(0) and T(L) are given. 
    #   In equilibrium, second derivative of T is zero
    #   At the interface, CI*dT/dx(x-) = CII*dT/dx(x+)

    # CASE I: No heat flux source
    # ---------------------------

    T0 = -10. # Temperature left boundary
    TL = 9. # Temperature right boundary

    CI = 30. # Heat flux coefficent left material
    CII = 150. # Heat flux coefficent right material

    dxL = 1e-3
    dxR = 4e-3

    N = 30
    NL = 20
    NR = N-NL
    a = NL*dxL # position of the interface

    cells_L = dxL*np.ones(NL)
    cells_R = dxR*np.ones(NR)

    dx_vals = np.concatenate((cells_L,cells_R))
    L = np.sum(dx_vals)

    n = (N,1,1)
    dx = (dx_vals, dxL, dxL)

    mesh = Mesh(n, dx)
    state = State(mesh)

    # Using stiffness matrix C11 component for heat flux coefficient
    C = state.Constant(torch.zeros((6,6)))
    C[...,0,0] = torch.ones(n)
    C[:NL] *= CI
    C[NL:] *= CII

    # Analytical solution T1(x) = A1*x+D1, T2(x) = A2*x+D2
    A1 = CII*(TL-T0)/(CI*(L-a)+CII*a)
    A2 = CI*A1/CII

    D1 = T0
    D2 = A1*a*(1.-CI/CII) + T0

    x = mesh.SpatialCoordinate()[0]
    state.ud = state.Constant((0.,0.,0.))
    state.ud[...,0] = x*A1+D1
    state.ud[NL:,...,0] = x[NL:]*A2+D2

    # no heat flux source
    state.m = state.Constant((1.,0.,0.))

    state.material = {"C" : C,
                      "lambda_100" : state.Constant(0.), # no heat flux source
                      "lambda_111" : state.Constant(0.), # no heat flux source
                      "A" : state.Constant(1.)} # no heat flux source

    ###################################

    llg = LLGWithLESolver([])
    llg._update_diff_data(state)

    # first derivative
    grad = llg.diff_data.gradient_ud[0][0]
    grad_eps = epsilon(state)[...,0]
    # second derivative
    laplace = llg._2nd_derivative_homogeneous_Neumann(state, 0, 0, [0,0])
    laplace_pbc = llg._2nd_derivative_with_pbc(state, 0, 0, [0,0])

    ###################################

    rel = 1e-6
    abs = rel*np.max([A1, A2])
    # epsilon and diff data contain same data
    torch.testing.assert_close(grad, grad_eps, rtol=rel, atol=abs)
    # compare first derivative
    g = torch.ones(n)
    g[:NL] *= A1 
    g[NL:] *= A2
    torch.testing.assert_close(grad, g, rtol=rel, atol=abs)
    # compare second derivative
    l = torch.zeros(n)[1:-1]
    torch.testing.assert_close(laplace[1:-1], l, rtol=0, atol=abs)
    torch.testing.assert_close(laplace_pbc[1:-1], l, rtol=0, atol=abs)

    # CASE II: Heat flux source
    # -------------------------

    q = 0.8*CII*TL/L

    A1 = (CII*(TL-T0)+q*(L-a)) / (CI*(L-a) + CII*a)
    A2 = (CI*A1-q)/CII

    D1 = T0
    D2 = TL-A2*L

    state.ud[...,0] = x*A1+D1
    state.ud[NL:,...,0] = x[NL:]*A2+D2

    # heat flux source
    state.material["lambda_100"] = state.Constant(0.)
    state.material["lambda_100"][NL:] = state.Constant(-q/CII)[NL:]

    llg._update_diff_data(state)

    # first derivative
    grad = llg.diff_data.gradient_ud[0][0]
    grad_eps = epsilon(state)[...,0]
    # second derivative
    laplace = llg._2nd_derivative_homogeneous_Neumann(state, 0, 0, [0,0])
    laplace_pbc = llg._2nd_derivative_with_pbc(state, 0, 0, [0,0])

    ###################################

    rel = 1e-6
    abs = rel*np.max([A1, A2])
    # epsilon and diff data contain same data
    torch.testing.assert_close(grad, grad_eps, rtol=rel, atol=abs)
    # compare first derivative
    g = torch.ones(n)
    g[:NL] *= A1 
    g[NL:] *= A2
    torch.testing.assert_close(grad, g, rtol=rel, atol=abs)
    # compare second derivative
    l = torch.zeros(n)[1:-1]
    torch.testing.assert_close(laplace[1:-1], l, rtol=0, atol=abs)
    torch.testing.assert_close(laplace_pbc[1:-1], l, rtol=0, atol=abs)

def test_bulk_modes():
    rho = 7500.
    C = C_isotropic(2e11, 0.3)
    C11 = C[0,0]
    C44 = C[3,3]

    cl = np.sqrt(C11/rho)
    ct = np.sqrt(C44/rho)

    # Setup mesh
    n = (40,1,1)
    dx = (4e-9, 4e-9, 4e-9)
    mesh = Mesh(n, dx, pbc=(1,1,1))

    # compute eigenvalues
    L = n[0]*dx[0]
    k = 2.*np.pi / L

    omega_l = k*cl
    omega_t = k*ct

    ev_l = -rho*omega_l**2.
    ev_t = -rho*omega_t**2.

    # setup state
    state = State(mesh)
    state.m = state.Constant((1.,0.,0.))

    state.material = {
        "Ms" : state.Constant(1./constants.mu_0),
        "A" : state.Constant(1.),
        "C" : state.Constant(C),
        "rho" : state.Constant(rho),
        "lambda_100" : state.Constant(0.),
        "lambda_111" : state.Constant(0.),
        "eta" : state.Constant(0),
    }

    x = mesh.SpatialCoordinate()[0]

    u0 = 5e-12
    state.ud = state.Constant((0.,0.,0.))
    state.ud[...,0] = u0*torch.sin(k*x)
    state.ud[...,2] = u0*torch.cos(k*x)

    state.pd = state.Constant((0.,0.,0.))
    state.pd[...,0] = k*u0*torch.cos(k*x)
    state.pd[...,2] = -k*u0*torch.sin(k*x)

    llg = LLGWithLESolver([])
    v = llg.dpd(state)

    rel = 1e-2
    abstol = rel*torch.max(v)

    torch.testing.assert_close(v[...,0], ev_l*state.ud[...,0], rtol=0, atol=abstol)
    torch.testing.assert_close(v[...,2], ev_t*state.ud[...,2], rtol=0, atol=abstol)
