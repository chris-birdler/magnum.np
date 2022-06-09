import numpy as np
from pyevtk.hl import gridToVTK
import os

__all__ = ["write_vtr"]

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
        gridToVTK(filename, x, y, z, cellData = {"f" : field.detach().cpu().numpy().copy()})
    else:
        gridToVTK(filename, x, y, z, cellData = {
                "f": (field[:,:,:,0].detach().cpu().numpy().copy(),
                      field[:,:,:,1].detach().cpu().numpy().copy(),
                      field[:,:,:,2].detach().cpu().numpy().copy())})
