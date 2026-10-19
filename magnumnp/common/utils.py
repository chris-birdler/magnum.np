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
import sys
import subprocess
from magnumnp.common import logging, Material
from magnumnp.common.io import write_vti, write_vtr

__all__ = ["complex_dtype", "normalize", "randM", "Expression", "accumulate_h", "set_precision"]


complex_dtype = {
    torch.float: torch.complex,
    torch.float32: torch.complex64,
    torch.float64: torch.complex128
    }


def set_precision(precision, force = False):
    r"""
    Set the global floating point precision.

    magnum.np defaults to double precision. Single precision halves the
    memory consumption and is considerably faster (especially on GPUs,
    where the FP64 throughput is a fraction of FP32), while the demag /
    Oersted kernels are still evaluated in double precision internally.
    Single precision is usually sufficient for LLG time integration;
    tightly converged energy minimization or comparing energies of nearly
    degenerate states may still require double precision.

    Needs to be called before the first :class:`Mesh` or :class:`State`
    is created, since the dtype is baked into the created tensors (which
    would silently mix precisions otherwise).

    :param precision: "single" / "float32" or "double" / "float64"
    :type precision: str
    :param force: skip the mesh-already-created check (use with care)
    :type force: bool, optional

    :Example:
        .. code::
            from magnumnp import *
            set_precision("single")

            mesh = Mesh(n, dx)
            state = State(mesh)
    """
    dtypes = {"single": torch.float32, "float32": torch.float32,
              "double": torch.float64, "float64": torch.float64}
    if precision not in dtypes:
        raise ValueError("Unknown precision '%s' (use one of %s)" % (precision, sorted(dtypes.keys())))

    from magnumnp.common import mesh as _mesh
    if _mesh._MESH_CREATED and not force:
        raise RuntimeError("set_precision() has to be called before the first Mesh/State is created "
                           "(existing tensors keep their dtype, which would silently mix precisions). "
                           "Use force=True to override.")

    torch.set_default_dtype(dtypes[precision])
    logging.info_green("[magnum.np] floating point precision: %s" % precision)


def accumulate_h(fields):
    r"""
    Helper function to sum up field contributions with minimal allocations.

    Consumes the iterable lazily and accumulates in place, so at most two
    full fields are alive at the same time (in contrast to sum([...]), which
    materializes all contributions and allocates a new tensor for each add).
    The first field is never modified, since field terms may return cached
    or shared tensors.
    """
    it = iter(fields)
    try:
        total = next(it)
    except StopIteration:
        raise ValueError("accumulate_h() requires at least one field contribution")
    accumulated = False
    for f in it:
        if not accumulated:
            total = total + f # allocate the accumulator; never mutate the first field
            accumulated = True
        else:
            total.add_(f)
    if not accumulated:
        total = total.clone() # single contribution: never hand out a term's own (possibly cached) tensor
    return total


def normalize(data):
    r"""
    Helper function to normalize vectorial data inplace
    """
    data /= torch.linalg.norm(data, dim = -1, keepdim = True)
    data.nan_to_num_(posinf=0, neginf=0)
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
    code = """
           import torch
           def get_gpu_with_least_memory():
               if not torch.cuda.is_available():
                   print(-1)
                   return
               num_gpus = torch.cuda.device_count()
               if num_gpus == 1:
                   print(0)
                   return
               gpu_memory = [torch.cuda.mem_get_info(i)[0] for i in range(num_gpus)]
               print(gpu_memory.index(max(gpu_memory)))
           get_gpu_with_least_memory()
           """

    # unindent code
    lines = code.splitlines()
    non_empty = [line for line in lines if line.strip()]
    min_indent = min((len(line) - len(line.lstrip()) for line in non_empty), default=0)
    code = "\n".join(line[min_indent:] for line in lines)

    # run in seperate process in order to prevent nvidia-smi show the process on all GPUs
    result = subprocess.check_output([sys.executable, "-c", code])
    return int(result.decode().strip())
