import torch
import numpy as np
from pyevtk.hl import gridToVTK
import pyvista as pv
import os
from . import Mesh

__all__ = ["write_vtr", "write_vti", "read_vti"]

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
#TODO: maybe move to state, or introduce decorated function
def write_vti(fields, filename, state = None):
    dirname = os.path.dirname(filename)
    if not os.path.isdir(dirname):
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
        n = state._mesh.n
        dx = state._mesh.dx
        origin = state._mesh.origin

    grid = pv.UniformGrid(dims = np.array(n) + 1,
                          spacing = dx,
                          origin = origin)

    for name, f in fields.items():
        if len(f.shape) == 4 and f.shape[-1] == 1:
            f = f[: ,:, :, 0]
        if len(f.shape) == 3:
            grid.cell_data.set_scalars(f.detach().cpu().numpy().flatten('F'), name)
        else:
            grid.cell_data.set_vectors(f.detach().cpu().numpy().reshape(-1,3,order='F'), name)

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
        f = torch.from_numpy(vals.reshape(dim, order="F"))
        fields[name] = f
    return mesh, fields


### deprecated functions (will be removed)
def write_vtr(field, filename, mesh = None, name="f"):
    dirname = os.path.dirname(filename)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    if mesh is None:
        n = field.shape[:3]
        dx = (1., 1., 1.)
    else:
        n = mesh.n
        dx = mesh.dx

    x = np.arange(0, (0.1 + n[0]) * dx[0], dx[0], dtype='float64')
    y = np.arange(0, (0.1 + n[1]) * dx[1], dx[1], dtype='float64')
    z = np.arange(0, (0.1 + n[2]) * dx[2], dx[2], dtype='float64')

    if field.shape[-1] == 1:
        field = field[...,0]

    if len(field.shape) == 3:
        gridToVTK(filename, x, y, z, cellData = {name: field.detach().cpu().numpy().copy()})
    else:
        gridToVTK(filename, x, y, z, cellData = {
                name: (field[:,:,:,0].detach().cpu().numpy().copy(),
                       field[:,:,:,1].detach().cpu().numpy().copy(),
                       field[:,:,:,2].detach().cpu().numpy().copy())})
