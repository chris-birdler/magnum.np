#
# This file is part of the magnum.np distribution
# (https://gitlab.com/magnum.np/magnum.np).
# Copyright (c) 2023 magnum.np team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

import torch 
import numpy as np 
from magnumnp.linear_elasticity.utils import gradient_with_pbc
from magnumnp.linear_elasticity.stress import sigma

__all__ = ["epsilon", "epsilon_m", "epsilon_el", "_get_C_jump_conditions", "_get_sigM_jump_conditions", "_get_B_jump_conditions"]

#@torch.compile
def epsilon(state, ud=None, iteration_depht = 1, second_order_boundary = False):
    if ud is None:
        ud = state.ud
    else:
        assert isinstance(ud, torch.Tensor) and ud.shape == state.mesh.n + (3,)

    mxl = torch.clone(state.m)
    mxr = torch.clone(state.m)
    myl = torch.clone(state.m)
    myr = torch.clone(state.m)
    mzl = torch.clone(state.m)
    mzr = torch.clone(state.m)

    dx_exp = state.mesh.dx_tensor[0].reshape(-1,1,1,1)
    dy_exp = state.mesh.dx_tensor[1].reshape(1,-1,1,1)
    dz_exp = state.mesh.dx_tensor[2].reshape(1,1,-1,1)

    A = state.material["A"][...,0]
    grad_mx = gradient_with_pbc(state.m[...,0], state.mesh, [0,1,2], [A, A, A], second_order_boundary=second_order_boundary)
    grad_my = gradient_with_pbc(state.m[...,1], state.mesh, [0,1,2], [A, A, A], second_order_boundary=second_order_boundary)
    grad_mz = gradient_with_pbc(state.m[...,2], state.mesh, [0,1,2], [A, A, A], second_order_boundary=second_order_boundary)

    n = state.mesh.n
    dmdx = torch.zeros(n+(3,))
    dmdx[...,0] = grad_mx[0]
    dmdx[...,1] = grad_my[0]
    dmdx[...,2] = grad_mz[0]

    dmdy = torch.zeros(n+(3,))
    dmdy[...,0] = grad_mx[1]
    dmdy[...,1] = grad_my[1]
    dmdy[...,2] = grad_mz[1]

    dmdz = torch.zeros(n+(3,))
    dmdz[...,0] = grad_mx[2]
    dmdz[...,1] = grad_my[2]
    dmdz[...,2] = grad_mz[2]

    mxl -= 0.5*dmdx*dx_exp
    mxr += 0.5*dmdx*dx_exp
    myl -= 0.5*dmdy*dy_exp
    myr += 0.5*dmdy*dy_exp
    mzl -= 0.5*dmdz*dz_exp
    mzr += 0.5*dmdz*dz_exp

    m_data = [mxl, mxr], [myl, myr], [mzl, mzr]
    C = _get_C_jump_conditions(state)
    Bl_sigM, Br_sigM = _get_sigM_jump_conditions(state, m_data)

    grad_ud_x = gradient_with_pbc(ud[...,0], state.mesh, [0,1,2], C[0], Bl_sigM[0], Br_sigM[0], second_order_boundary=second_order_boundary)
    grad_ud_y = gradient_with_pbc(ud[...,1], state.mesh, [0,1,2], C[1], Bl_sigM[1], Br_sigM[1], second_order_boundary=second_order_boundary)
    grad_ud_z = gradient_with_pbc(ud[...,2], state.mesh, [0,1,2], C[2], Bl_sigM[2], Br_sigM[2], second_order_boundary=second_order_boundary)

    if (iteration_depht > 0):
        for iter in range(iteration_depht):
            gradient_data = grad_ud_x, grad_ud_y, grad_ud_z

            Bl_eps, Br_eps = _get_B_jump_conditions(state, gradient_data)

            Bxl = [Bl_eps[0][0] + Bl_sigM[0][0], Bl_eps[0][1] + Bl_sigM[0][1], Bl_eps[0][2] + Bl_sigM[0][2]]
            Byl = [Bl_eps[1][0] + Bl_sigM[1][0], Bl_eps[1][1] + Bl_sigM[1][1], Bl_eps[1][2] + Bl_sigM[1][2]]
            Bzl = [Bl_eps[2][0] + Bl_sigM[2][0], Bl_eps[2][1] + Bl_sigM[2][1], Bl_eps[2][2] + Bl_sigM[2][2]]

            Bxr = [Br_eps[0][0] + Br_sigM[0][0], Br_eps[0][1] + Br_sigM[0][1], Br_eps[0][2] + Br_sigM[0][2]]
            Byr = [Br_eps[1][0] + Br_sigM[1][0], Br_eps[1][1] + Br_sigM[1][1], Br_eps[1][2] + Br_sigM[1][2]]
            Bzr = [Br_eps[2][0] + Br_sigM[2][0], Br_eps[2][1] + Br_sigM[2][1], Br_eps[2][2] + Br_sigM[2][2]]

            grad_ud_x[:] = gradient_with_pbc(ud[...,0], state.mesh, [0,1,2], C[0], Bxl, Bxr, second_order_boundary=second_order_boundary)[:]
            grad_ud_y[:] = gradient_with_pbc(ud[...,1], state.mesh, [0,1,2], C[1], Byl, Byr, second_order_boundary=second_order_boundary)[:]
            grad_ud_z[:] = gradient_with_pbc(ud[...,2], state.mesh, [0,1,2], C[2], Bzl, Bzr, second_order_boundary=second_order_boundary)[:]

    n = state.mesh.n
    eps = torch.zeros(n+(6,))
    eps[:,:,:,0] = grad_ud_x[0] # eps_0 = eps_xx
    eps[:,:,:,1] = grad_ud_y[1] # eps_1 = eps_yy
    eps[:,:,:,2] = grad_ud_z[2] # eps_2 = eps_zz
    eps[:,:,:,3] = grad_ud_y[2] + grad_ud_z[1] # eps_3 = 2*eps_yz
    eps[:,:,:,4] = grad_ud_x[2] + grad_ud_z[0] # eps_4 = 2*eps_xz
    eps[:,:,:,5] = grad_ud_x[1] + grad_ud_y[0] # eps_5 = 2*eps_xy

    return eps

@torch.compile
def epsilon_m(state, m=None):
    if m==None:
        m = state.m

    lambda_100 = state.material["lambda_100"][:,:,:,0]
    lambda_111 = state.material["lambda_111"][:,:,:,0]

    n = state.mesh.n
    eps_m = torch.zeros(n+(6,))

    eps_m[:,:,:,0] = (3./2.)*lambda_100*(m[:,:,:,0]**2. - 1./3.)
    eps_m[:,:,:,1] = (3./2.)*lambda_100*(m[:,:,:,1]**2. - 1./3.)
    eps_m[:,:,:,2] = (3./2.)*lambda_100*(m[:,:,:,2]**2. - 1./3.)
    eps_m[:,:,:,3] = 3.*lambda_111*m[:,:,:,1]*m[:,:,:,2]
    eps_m[:,:,:,4] = 3.*lambda_111*m[:,:,:,0]*m[:,:,:,2]
    eps_m[:,:,:,5] = 3.*lambda_111*m[:,:,:,0]*m[:,:,:,1]

    return eps_m

def epsilon_el(state, ud=None, m=None):
    return epsilon(state, ud) - epsilon_m(state, m)

def _get_C_jump_conditions(state):
    C = state.material["C"]
    C11 = C[...,0,0]
    C22 = C[...,1,1]
    C33 = C[...,2,2]
    C44 = C[...,3,3]
    C55 = C[...,4,4]
    C66 = C[...,5,5]

    # for x derivatives
    Cxx = C11
    Cyx = C66
    Czx = C55

    # for y derivatives
    Cxy = C66
    Cyy = C22
    Czy = C44

    # for z derivatives
    Cxz = C55
    Cyz = C44
    Czz = C33

    return [Cxx, Cxy, Cxz], [Cyx, Cyy, Cyz], [Czx, Czy, Czz]

def _get_sigM_jump_conditions(state, m_data):
    mxl, mxr = m_data[0]
    myl, myr = m_data[1]
    mzl, mzr = m_data[2]
    
    # for x derivatives
    eps_m = epsilon_m(state, mxl)
    sig_m = sigma(state,eps_m)
    Bx_xl = -sig_m[...,0]
    By_xl = -sig_m[...,5]
    Bz_xl = -sig_m[...,4]

    eps_m[:] = epsilon_m(state, mxr)
    sig_m[:] = sigma(state,eps_m)
    Bx_xr = -sig_m[...,0]
    By_xr = -sig_m[...,5]
    Bz_xr = -sig_m[...,4]

    # for y derivatives
    eps_m[:] = epsilon_m(state, myl)
    sig_m[:] = sigma(state,eps_m)
    Bx_yl = -sig_m[...,5]
    By_yl = -sig_m[...,1]
    Bz_yl = -sig_m[...,4]

    eps_m[:] = epsilon_m(state, myr)
    sig_m[:] = sigma(state,eps_m)
    Bx_yr = -sig_m[...,5]
    By_yr = -sig_m[...,1]
    Bz_yr = -sig_m[...,4]

    # for z derivatives
    eps_m[:] = epsilon_m(state, mzl)
    sig_m[:] = sigma(state,eps_m)
    Bx_zl = -sig_m[...,4]
    By_zl = -sig_m[...,3]
    Bz_zl = -sig_m[...,2]

    eps_m[:] = epsilon_m(state, mzr)
    sig_m[:] = sigma(state,eps_m)
    Bx_zr = -sig_m[...,4]
    By_zr = -sig_m[...,3]
    Bz_zr = -sig_m[...,2]

    Bl = [[Bx_xl, Bx_yl, Bx_zl], [By_xl, By_yl, By_zl], [Bz_xl, Bz_yl, Bz_zl]]
    Br = [[Bx_xr, Bx_yr, Bx_zr], [By_xr, By_yr, By_zr], [Bz_xr, Bz_yr, Bz_zr]]

    return Bl, Br

def _get_B_jump_conditions(state, gradient_data):
    gux, guy, guz = gradient_data

    C = state.material["C"]
    C12 = C[...,0,1]
    C13 = C[...,0,2]
    C23 = C[...,1,2]
    C44 = C[...,3,3]
    C55 = C[...,4,4]
    C66 = C[...,5,5]

    def harmonic_mean(g, C, shift, dim):
        a = C*g
        mean = (a+torch.roll(a, shift, dim)) / (C + torch.roll(C, shift, dim))
        return mean.nan_to_num(posinf=0, neginf=0)
    
    dyux_xl = harmonic_mean(gux[1], C66, 1, 0)
    dzux_xl = harmonic_mean(gux[2], C55, 1, 0)
    dxux_yl = harmonic_mean(gux[0], C12, 1, 1)
    dxux_zl = harmonic_mean(gux[0], C13, 1, 2)

    dxuy_yl = harmonic_mean(guy[0], C66, 1, 1)
    dzuy_yl = harmonic_mean(guy[2], C44, 1, 1)
    dyuy_xl = harmonic_mean(guy[1], C12, 1, 0)
    dyuy_zl = harmonic_mean(guy[1], C23, 1, 2)

    dxuz_zl = harmonic_mean(guz[0], C55, 1, 2)
    dyuz_zl = harmonic_mean(guz[1], C44, 1, 2)
    dzuz_xl = harmonic_mean(guz[2], C13, 1, 0)
    dzuz_yl = harmonic_mean(guz[2], C23, 1, 1)

    dyux_xr = torch.roll(dyux_xl, -1, 0)
    dzux_xr = torch.roll(dzux_xl, -1, 0)
    dxux_yr = torch.roll(dxux_yl, -1, 1)
    dxux_zr = torch.roll(dxux_zl, -1, 2)

    dxuy_yr = torch.roll(dxuy_yl, -1, 1)
    dzuy_yr = torch.roll(dzuy_yl, -1, 1)
    dyuy_xr = torch.roll(dyuy_xl, -1, 0)
    dyuy_zr = torch.roll(dyuy_zl, -1, 2)

    dxuz_zr = torch.roll(dxuz_zl, -1, 2)
    dyuz_zr = torch.roll(dyuz_zl, -1, 2)
    dzuz_xr = torch.roll(dzuz_xl, -1, 0)
    dzuz_yr = torch.roll(dzuz_yl, -1, 1)

    # for x derivatives
    Bx_xl = C12*dyuy_xl + C13*dzuz_xl
    By_xl = C66*dyux_xl
    Bz_xl = C55*dzux_xl

    Bx_xr = C12*dyuy_xr + C13*dzuz_xr
    By_xr = C66*dyux_xr
    Bz_xr = C55*dzux_xr

    # for y derivatives
    Bx_yl = C66*dxuy_yl
    By_yl = C12*dxux_yl + C23*dzuz_yl
    Bz_yl = C44*dzuy_yl

    Bx_yr = C66*dxuy_yr
    By_yr = C12*dxux_yr + C23*dzuz_yr
    Bz_yr = C44*dzuy_yr

    # for z derivatives
    Bx_zl = C55*dxuz_zl
    By_zl = C44*dyuz_zl
    Bz_zl = C13*dxux_zl + C23*dyuy_zl

    Bx_zr = C55*dxuz_zr
    By_zr = C44*dyuz_zr
    Bz_zr = C13*dxux_zr + C23*dyuy_zr

    Bl = [[Bx_xl, Bx_yl, Bx_zl], [By_xl, By_yl, By_zl], [Bz_xl, Bz_yl, Bz_zl]]
    Br = [[Bx_xr, Bx_yr, Bx_zr], [By_xr, By_yr, By_zr], [Bz_xr, Bz_yr, Bz_zr]]
    
    return Bl, Br