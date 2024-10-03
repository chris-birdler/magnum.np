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

__all__ = ["complex_dtype", "normalize", "randM", "Expression"]


complex_dtype = {
    torch.float: torch.complex,
    torch.float32: torch.complex64,
    torch.float64: torch.complex128
    }


def normalize(data):
    r"""
    Helper function to normalize vectorial data inplace
    """
    data /= torch.linalg.norm(data, dim = -1, keepdim = True)
    data[...] = torch.nan_to_num(data, posinf=0, neginf=0)
    return data

def randM(data):
    r"""
    Helper function to generate uniform distibution on the unit-sphere
    """
    if data.shape[-1] != 3:
        raise ValueError("Input tensor's last dimension needs to be 3 (shape='%s')" % str(data.shape))

    theta = 2.*torch.pi*torch.rand(data.shape[:-1])
    phi = torch.acos(2.*torch.rand(data.shape[:-1])-1.)

    data[...,0] = torch.sin(phi) * torch.cos(theta)
    data[...,1] = torch.sin(phi) * torch.sin(theta)
    data[...,2] = torch.cos(phi)
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


def get_gpu_with_least_memory():
    if not torch.cuda.is_available():
        return -1

    import pynvml
    pynvml.nvmlInit()
    num_gpus = pynvml.nvmlDeviceGetCount()

    gpu_memory = []
    for i in range(num_gpus):
        handle = pynvml.nvmlDeviceGetHandleByIndex(i)
        mem_info = pynvml.nvmlDeviceGetMemoryInfo(handle)
        gpu_memory.append(mem_info.used)

    pynvml.nvmlShutdown()
    return gpu_memory.index(min(gpu_memory))
