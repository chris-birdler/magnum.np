import numpy as np
from pyevtk.hl import gridToVTK

def write_vtr(field, filename, mesh = None):
    if mesh is None:
        n = field.shape[:3]
        dx = (1., 1., 1.)
    else:
        n = mesh.n
        dx = mesh.dx

    # save result
    x = np.arange(0, (0.1 + n[0]) * dx[0], dx[0], dtype='float64')
    y = np.arange(0, (0.1 + n[1]) * dx[1], dx[1], dtype='float64')
    z = np.arange(0, (0.1 + n[2]) * dx[2], dx[2], dtype='float64')

    # scalar data
    if len(field.shape) == 3:
        gridToVTK(filename, x, y, z, cellData = {"f" : field.detach().cpu().numpy().copy()})
    else:
        gridToVTK(filename, x, y, z, cellData = {
                "x": field[:,:,:,0].detach().cpu().numpy().copy(),
                "y": field[:,:,:,1].detach().cpu().numpy().copy(),
                "z": field[:,:,:,2].detach().cpu().numpy().copy()})
