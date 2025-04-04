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

__all__ = ["epsilon", "epsilon_m", "epsilon_el", "_get_jump_conditions"]

#@torch.compile
def epsilon(state, ud=None):
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
    grad_mx = gradient_with_pbc(state.m[...,0], state.mesh, [0,1,2], [A, A, A])
    grad_my = gradient_with_pbc(state.m[...,1], state.mesh, [0,1,2], [A, A, A])
    grad_mz = gradient_with_pbc(state.m[...,2], state.mesh, [0,1,2], [A, A, A])

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

    C, Bl, Br = _get_jump_conditions(state, mxl, mxr, myl, myr, mzl, mzr)

    grad_ud_x = gradient_with_pbc(ud[...,0], state.mesh, [0,1,2], C[0], Bl[0], Br[0])
    grad_ud_y = gradient_with_pbc(ud[...,1], state.mesh, [0,1,2], C[1], Bl[1], Br[1])
    grad_ud_z = gradient_with_pbc(ud[...,2], state.mesh, [0,1,2], C[2], Bl[2], Br[2])

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
    eps_m = torch.zeros((n[0], n[1], n[2], 6), device=m.device)

    eps_m[:,:,:,0] = (3./2.)*lambda_100*(m[:,:,:,0]**2. - 1./3.)
    eps_m[:,:,:,1] = (3./2.)*lambda_100*(m[:,:,:,1]**2. - 1./3.)
    eps_m[:,:,:,2] = (3./2.)*lambda_100*(m[:,:,:,2]**2. - 1./3.)
    eps_m[:,:,:,3] = 3.*lambda_111*m[:,:,:,1]*m[:,:,:,2]
    eps_m[:,:,:,4] = 3.*lambda_111*m[:,:,:,0]*m[:,:,:,2]
    eps_m[:,:,:,5] = 3.*lambda_111*m[:,:,:,0]*m[:,:,:,1]

    return eps_m

def epsilon_el(state, ud=None, m=None):
    return epsilon(state, ud) - epsilon_m(state, m)


def _get_jump_conditions(state, mxl, mxr, myl, myr, mzl, mzr):
    C = state.material["C"]
    C11 = C[...,0,0]
    C22 = C[...,1,1]
    C33 = C[...,2,2]
    C44 = C[...,3,3]
    C55 = C[...,4,4]
    C66 = C[...,5,5]

    C15 = C[...,0,4]
    C16 = C[...,0,5]
    C24 = C[...,1,3]
    C26 = C[...,1,5]
    C35 = C[...,2,4]
    C34 = C[...,2,3]
    C45 = C[...,3,4]
    C46 = C[...,3,5]
    C56 = C[...,4,5]

    # for x derivatives
    Cxx = C11+C15+C16
    Cyx = C16+C66+C56
    Czx = C15+C56+C55

    eps_m = epsilon_m(state, mxl)
    sig_m = sigma(state,eps_m)
    # B_xl = sig_m[...,0] + sig_m[...,5] + sig_m[...,4]
    Bx_xl = sig_m[...,0]
    By_xl = sig_m[...,5]
    Bz_xl = sig_m[...,4]
    """
    epsMXX = eps_m[...,0]
    epsMXY = eps_m[...,5]
    epsMXZ = eps_m[...,4]
    
    Bx_xl = C11*epsMXX + C15*epsMXZ + C16*epsMXY
    By_xl = C16*epsMXX + C66*epsMXY + C56*epsMXZ
    Bz_xl = C15*epsMXX + C56*epsMXY + C55*epsMXZ
    """

    eps_m = epsilon_m(state, mxr)
    sig_m = sigma(state,eps_m)
    #B_xr = sig_m[...,0] + sig_m[...,5] + sig_m[...,4]
    Bx_xr = sig_m[...,0]
    By_xr = sig_m[...,5]
    Bz_xr = sig_m[...,4]
    """
    epsMXX = eps_m[...,0]
    epsMXY = eps_m[...,5]
    epsMXZ = eps_m[...,4]

    Bx_xr = C11*epsMXX + C15*epsMXZ + C16*epsMXY
    By_xr = C16*epsMXX + C66*epsMXY + C56*epsMXZ
    Bz_xr = C15*epsMXX + C56*epsMXY + C55*epsMXZ
    """

    # for y derivatives
    Cxy = C66+C26+C46
    Cyy = C26+C22+C24
    Czy = C46+C24+C44

    eps_m = epsilon_m(state, myl)
    sig_m = sigma(state,eps_m)
    #B_yl = sig_m[...,1] + sig_m[...,5] + sig_m[...,3]
    Bx_yl = sig_m[...,5]
    By_yl = sig_m[...,1]
    Bz_yl = sig_m[...,4]
    """
    epsMYY = eps_m[...,1]
    epsMXY = eps_m[...,5]
    epsMYZ = eps_m[...,3] 

    Bx_yl = C66*epsMXY + C26*epsMYY + C46*epsMYZ
    By_yl = C26*epsMXY + C22*epsMYY + C24*epsMYZ
    Bz_yl = C46*epsMXY + C24*epsMYY + C44*epsMYZ
    """

    eps_m = epsilon_m(state, myr)
    sig_m = sigma(state,eps_m)
    #B_yr = sig_m[...,1] + sig_m[...,5] + sig_m[...,3]
    Bx_yr = sig_m[...,5]
    By_yr = sig_m[...,1]
    Bz_yr = sig_m[...,4]
    """
    epsMYY = eps_m[...,1]
    epsMXY = eps_m[...,5]
    epsMYZ = eps_m[...,3] 

    Bx_yr = C66*epsMXY + C26*epsMYY + C46*epsMYZ
    By_yr = C26*epsMXY + C22*epsMYY + C24*epsMYZ
    Bz_yr = C46*epsMXY + C24*epsMYY + C44*epsMYZ
    """

    # for z derivatives
    Cxz = C55+C45+C35
    Cyz = C45+C44+C34
    Czz = C35+C34+C33

    eps_m = epsilon_m(state, mzl)
    sig_m = sigma(state,eps_m)
    #B_zl = sig_m[...,2] + sig_m[...,4] + sig_m[...,3]
    Bx_zl = sig_m[...,4]
    By_zl = sig_m[...,3]
    Bz_zl = sig_m[...,2]
    """
    epsMZZ = eps_m[...,2]
    epsMXZ = eps_m[...,4]
    epsMYZ = eps_m[...,3] 

    Bx_zl = C55*epsMXZ + C45*epsMYZ + C35*epsMZZ
    By_zl = C45*epsMXZ + C44*epsMYZ + C34*epsMZZ
    Bz_zl = C35*epsMXZ + C34*epsMYZ + C33*epsMZZ
    """

    eps_m = epsilon_m(state, mzr)
    sig_m = sigma(state,eps_m)
    #B_zr = sig_m[...,2] + sig_m[...,4] + sig_m[...,3]
    Bx_zr = sig_m[...,4]
    By_zr = sig_m[...,3]
    Bz_zr = sig_m[...,2]
    """
    epsMZZ = eps_m[...,2]
    epsMXZ = eps_m[...,4]
    epsMYZ = eps_m[...,3] 

    Bx_zr = C55*epsMXZ + C45*epsMYZ + C35*epsMZZ
    By_zr = C45*epsMXZ + C44*epsMYZ + C34*epsMZZ
    Bz_zr = C35*epsMXZ + C34*epsMYZ + C33*epsMZZ
    """

    C = [[Cxx, Cxy, Cxz], [Cyx, Cyy, Cyz], [Czx, Czy, Czz]]
    Bl = [[-Bx_xl, -Bx_yl, -Bx_zl], [-By_xl, -By_yl, -By_zl], [-Bz_xl, -Bz_yl, -Bz_zl]]
    Br = [[-Bx_xr, -Bx_yr, -Bx_zr], [-By_xr, -By_yr, -By_zr], [-Bz_xr, -Bz_yr, -Bz_zr]]
    
    return C, Bl, Br