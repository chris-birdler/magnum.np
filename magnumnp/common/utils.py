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
from magnumnp.common import logging, Material
from magnumnp.common.io import write_vti, write_vtr

__all__ = ["complex_dtype", "avg", "normalize", "Expression"]


complex_dtype = {
    torch.float: torch.complex,
    torch.float32: torch.complex64,
    torch.float64: torch.complex128
    }


def avg(data, cell_volumes = None, dim=(0,1,2)):
    r"""
    Average over spatial dimensions of tensor fields.

    :param data: tensor field to average
    :type A: :class:`Tensor`
    :param dim: dimensions to average over
    :type dim: tuple, optional
    :param cell_volumes: volume of each cell (required only in case of non-equidistant meshes)
    :type cell_volumes: :class:`Tensor`, optional
    
    :Examples:

    .. code::
        Ms_avg = avg(state.material["Ms"])
        m_avg = avg(state.m)
    """
    if cell_volumes == None:
        if data.dim() <= 1: # e.g. [0,0,1]
            return data
        elif data.dim() == 2: # state.m[domain]
            return data.mean(dim=0)
        else:                 # [nx,ny,nz,...]
            return data.mean(dim=dim)
    else: # non-equidistant
        if data.dim() <= 1: # e.g. [0,0,1]
            return data
        elif data.dim() == 2: # state.m[domain]
            return (data * cell_volumes).sum(dim=0) / cell_volumes.sum(dim=0)
        elif data.dim() == 3: # [nx,ny,nz]
            return (data * cell_volumes.squeeze(-1)).sum(dim=dim) / cell_volumes.sum()
        else:                 # [nx,ny,nz,...]
            return (data * cell_volumes).sum(dim=dim) / cell_volumes.sum()


def normalize(data):
    data /= torch.linalg.norm(data, dim = -1, keepdim = True)
    data[...] = torch.nan_to_num(data, posinf=0, neginf=0)
    return data


def Expression(comps):
    r"""
    Helper function to create scalar or vector fields by stacking the corresponding components.

    :param comps: list of components or single torch.Tensor
    :type comps: list, :class:`Tensor`

    :Examples:

    .. code::
        n  = (1, 1, 10)
        dx = (2e-9, 2e-9, 2e-9)
        mesh = Mesh(n, dx)
        x, y, z = mesh.SpatialCoordinate()
        
        Ms = Expression(x)            # create scalar field
        Ku_axis = Expression((x,y,z)) # create vector field
    """
    if isinstance(comps, torch.Tensor):
        comps = [comps]
    return torch.stack(comps, dim=-1)
