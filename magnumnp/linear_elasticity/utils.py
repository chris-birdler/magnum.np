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

__all__ = ["gradient_with_pbc",
           "first_derivative_with_pbc"]

def _get_diff_slices(ndim, xdim):
    # ndim ... number of dimensions
    # xdim ... dimension that is used for differentiation
    slice_0 = [slice(None)] * ndim # for a nx,ny,nz,nd tensor, this is eqivalent to :,:,:,:
    slice_1 = [slice(None)] * ndim
    slice_0[xdim] = slice(0, -1)
    slice_1[xdim] = slice(1, None)
    slice_0 = tuple(slice_0)
    slice_1 = tuple(slice_1)

    return slice_0, slice_1

def _get_derivative_on_boundary(f0,f1,f2,h0,h1):
    g_bdr = -f2*h0**2. + f1*(h0+h1)**2. - f0*h1*(2*h0+h1)
    g_bdr /= h0*h1*(h0 + h1)
    return g_bdr

def gradient_with_pbc(f, mesh, dim=[0,1,2], C=None, Bl=None, Br=None, slices=None, second_order_boundary=False):
    """
    f      ... the gradient of this tensor is calculated 
    C,B    ... assumed to be constants on cells, gives jump conditions of the type
            lim x->a C_0 * df(x)/dx + B0 = lim a<-x C_1 * df(x)/dx + B1 
    mesh   ... the mesh on which f lives
    dim    ... the dimensions in which the first derivatives are calculated
    slices ... slicing of f, default to :,:,:
    """
    # validate input for dim
    dims = [dim] if isinstance(dim, int) else dim
    assert isinstance(dims, list) and all(isinstance(d, int) for d in dims)

    # slices 
    if slices == None:
        slices = tuple([slice(None)]*3)
    assert isinstance(slices, tuple) and all(isinstance(s, slice) for s in slices)

    # get jump condition data
    if (Bl is None) == (Br is None):
        if Bl is None:
            Bl = [torch.zeros(mesh.n)]*len(dims)
            Br = [torch.zeros(mesh.n)]*len(dims)
        Bl = Bl if isinstance(Bl, list) else [Bl]
        Br = Br if isinstance(Br, list) else [Br]
        assert len(Bl)==len(dims) and len(Br)==len(dims)

    else:
        raise Exception("If Bl is set, Br needs to be set too and vice-versa.")
    
    if C is None:
        C = [torch.ones(mesh.n)]*len(dims)
    C = C if isinstance(C, list) else [C]
    assert len(C)==len(dims)

    # slice data 
    fs = f[slices]
    n_f = fs.shape # used nx, ny, nz + data point dimensions
    
    # get expanded dx and slice it
    dx = mesh.dx_tensor[0][slices[0]].reshape(-1,1,1)
    dy = mesh.dx_tensor[1][slices[1]].reshape(1,-1,1)
    dz = mesh.dx_tensor[2][slices[2]].reshape(1,1,-1)
    dx_exp = dx,dy,dz

    # get gradient
    output = []
    for i in range(len(dims)):
        d = dims[i]
        # if there is only one data point in this dimension, the first derivative is 0
        if (n_f[d]) < 2:
            output.append(torch.zeros(n_f))
        # case no pbc: boundary node order is 1, obtained from forward and backward difference
        # this method is also used when the slice limits the input to within the mesh
        elif ((mesh.pbc[d] == 0) or (fs.shape[d] != f.shape[d])):
            g = first_derivative_with_jump_conditions(fs, C[i][slices], Bl[i][slices], Br[i][slices], dx_exp[d], d)
            
            if (n_f[d] > 2) and second_order_boundary:
                dx_s = dx_exp[d]
                if (d == 0):
                    f0 = fs[0]
                    f1 = fs[1]
                    f2 = fs[2]
                    h0 = 0.5*(dx_s[0]+dx_s[1])
                    h1 = 0.5*(dx_s[1]+dx_s[2])
                    g[0] = _get_derivative_on_boundary(f0, f1, f2, h0, h1)

                    f0 = fs[-1]
                    f1 = fs[-2]
                    f2 = fs[-3]
                    h0 = 0.5*(dx_s[-1]+dx_s[-2])
                    h1 = 0.5*(dx_s[-2]+dx_s[-3])
                    g[-1] = -_get_derivative_on_boundary(f0, f1, f2, h0, h1)

                if (d == 1):
                    f0 = fs[:,0]
                    f1 = fs[:,1]
                    f2 = fs[:,2]
                    h0 = 0.5*(dx_s[:,0]+dx_s[:,1])
                    h1 = 0.5*(dx_s[:,1]+dx_s[:,2])
                    g[:,0] = _get_derivative_on_boundary(f0, f1, f2, h0, h1)

                    f0 = fs[:,-1]
                    f1 = fs[:,-2]
                    f2 = fs[:,-3]
                    h0 = 0.5*(dx_s[:,-1]+dx_s[:,-2])
                    h1 = 0.5*(dx_s[:,-2]+dx_s[:,-3])
                    g[:,-1] = -_get_derivative_on_boundary(f0, f1, f2, h0, h1)

                if (d == 2):
                    f0 = fs[:,:,0]
                    f1 = fs[:,:,1]
                    f2 = fs[:,:,2]
                    h0 = 0.5*(dx_s[:,:,0]+dx_s[:,:,1])
                    h1 = 0.5*(dx_s[:,:,1]+dx_s[:,:,2])
                    g[:,:,0] = _get_derivative_on_boundary(f0, f1, f2, h0, h1)

                    f0 = fs[:,:,-1]
                    f1 = fs[:,:,-2]
                    f2 = fs[:,:,-3]
                    h0 = 0.5*(dx_s[:,:,-1]+dx_s[:,:,-2])
                    h1 = 0.5*(dx_s[:,:,-2]+dx_s[:,:,-3])
                    g[:,:,-1] = -_get_derivative_on_boundary(f0, f1, f2, h0, h1)

            output.append(g)
        # case pbc
        else:
            g = first_derivative_with_jump_conditions_and_pbc(fs, C[i][slices], Bl[i][slices], Br[i][slices], dx_exp[d], d)
            output.append(g)

    return output

def first_derivative_with_jump_conditions(f, C, Bl, Br, dx, dim):   
    slice_0, slice_1 = _get_diff_slices(f.ndim, dim)

    # for non-periodic meshes: 
    # the virtual node outside the boundary is assumed to have the same dx as the boundary node

    h = 0.5*(dx[slice_0] + dx[slice_1])
    h_sum = torch.zeros(dx.shape)
    h_sum[slice_0] = h 
    h_sum[slice_1] += h

    # other aliases
    C0 = C[slice_0]
    C1 = C[slice_1]

    B0 = Br[slice_0]
    B1 = Bl[slice_1]

    f0 = f[slice_0]
    f1 = f[slice_1]

    # output gradient component should be equivalent to torch.gradient if C and B are constant
    g = torch.zeros(f.shape)
    
    denom = (C0*dx[slice_1]+C1*dx[slice_0])

    # backward difference
    g[slice_1] += h*(2.*C0*(f1-f0) + dx[slice_0]*(B0-B1)) / denom

    # forward difference
    g[slice_0] += h*(2.*C1*(f1-f0) + dx[slice_1]*(B1-B0)) / denom

    # weighted mean between forward and backward difference 
    return g.nan_to_num(posinf=0, neginf=0)/h_sum

def first_derivative_with_jump_conditions_and_pbc(f, C, Bl, Br, dx, dim):
    fl = torch.roll(f, 1, dims=dim)
    fr = torch.roll(f, -1, dims=dim)

    Cl = torch.roll(C, 1, dims=dim)
    Cr = torch.roll(C, -1, dims=dim)

    dxl = torch.roll(dx, 1, dims=dim)
    dxr = torch.roll(dx, -1, dims=dim)

    # spacing between the data points
    hl = 0.5*(dx+dxl)
    hr = 0.5*(dx+dxr)

    # first-order accurate approximations of first derivative
    g_bkw = (2.*Cl*(f-fl) + dxl*(torch.roll(Br,  1, dims=dim)-Bl)) / (Cl*dx+C*dxl)
    g_fwd = (2.*Cr*(fr-f) + dxr*(torch.roll(Bl, -1, dims=dim)-Br)) / (C*dxr+Cr*dx)
    
    g_bkw = g_bkw.nan_to_num(posinf=0, neginf=0)
    g_fwd = g_fwd.nan_to_num(posinf=0, neginf=0)

    # weighted mean between forward and backward difference 
    # this is second-order accurate, when C and B are constant
    return (hr*g_bkw + hl*g_fwd)/(hl+hr)

@torch.compile
def first_derivative_with_pbc(f, dx_expanded, dim):
    fl = torch.roll(f, 1, dims=dim)
    fr = torch.roll(f, -1, dims=dim)

    hl = 0.5*(dx_expanded + torch.roll(dx_expanded, 1, dims=dim))
    hr = torch.roll(hl, -1, dims=dim)

    return  (hl**2.*fr - hr**2.*fl + (hr**2.-hl**2.)*f) / (hr*hl**2. + hl*hr**2)