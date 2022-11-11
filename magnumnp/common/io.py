import torch
import numpy as np
import pyvista as pv
import os
from . import Mesh, DecoratedTensor

__all__ = ["write_vti", "read_vti"]

r"""
write vti files (compressed) using pyvista

:Examples:

  .. code::
    # write single scalar or vector
    write_vti(state.material.Ms, "scalar.vti")
    write_vti(state.m, "vector.vti")

    # use dictinary or list
    write_vti([state.m, h], "list.vti")
    write_vti({'m':state.m, 'h':h}, "dict.vti")
"""
def write_vti(fields, filename, state = None):
    dirname = os.path.dirname(filename)
    if dirname and not os.path.isdir(dirname):
        os.makedirs(dirname)

    if not (isinstance(fields, list) or isinstance(fields, dict)):
        fields = [fields]
    if not isinstance(fields, dict):
        fields = {"f%03d"%i:f for (i,f) in enumerate(fields)}

    if state is None:
        n = list(fields.values())[0].shape[:3]
        dx = (1., 1., 1.)
        origin = (0., 0., 0.)
    else:
        n = state.mesh.n
        dx = state.mesh.dx
        origin = state.mesh.origin

    grid = pv.UniformGrid(dims = np.array(n) + 1,
                          spacing = dx,
                          origin = origin)

    for name, f in fields.items():
        if len(f.shape) == 0 or len(f.shape) == 1: # expand constant tensor to tensorfield
            f = f.expand(n + f.shape)
        if len(f.shape) == 4 and f.shape[-1] == 1: # remove dim for scalar field (nx,ny,nz,1) => (nx,ny,nz)
            f = f[: ,:, :, 0]
        if len(f.shape) == 3: # scalar data
            grid.cell_data.set_array(f.detach().cpu().numpy().flatten('F'), name)
        elif len(f.shape) == 4: # vector data
            grid.cell_data.set_array(f.detach().cpu().numpy().reshape(-1,3,order='F'), name)
        else:
            raise ValueError("write_vti: unsupported data format (", name, f.shape, ")")

    grid.save(filename)


def read_vti(filename):
    fields = {}
    data = pv.read(filename)

    mesh = Mesh(np.array(data.dimensions)-1, data.spacing, data.origin)

    for name in data.array_names:
        f = data.get_array(name)
        vals = data.get_array(name)
        if len(vals.shape) == 1:
            dim = mesh.n
        else:
            dim = mesh.n + (vals.shape[-1],)
        f = torch.from_numpy(vals.reshape(dim, order="F")).as_subclass(DecoratedTensor)
        fields[name] = f
    return mesh, fields
