import torch
import numpy as np
import scipy
import pyvista as pv
import os
from . import Mesh, DecoratedTensor

__all__ = ["write_vti", "read_vti", "read_image"]

def write_vti(fields, filename, state = None):
    r"""
    write vti files (compressed) using pyvista

    :param list/dict/:class:`torch.Tensor`: List of Fields to be written

    :Examples:
      .. code::
        # write single scalar or vector
        write_vti(state.material.Ms, "scalar.vti")
        write_vti(state.m, "vector.vti")

        # use dictinary or list
        write_vti([state.m, h], "list.vti")
        write_vti({'m':state.m, 'h':h}, "dict.vti")
    """
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
    r"""
    Read vti files using pyvista

    :param str filename: Filename to be read
    :return :class:`Mesh` & dict: Mesh object and dictionary containing all data tensors

    :Examples:
      .. code::
        mesh, fields = read_vti("m0.vti")
    """
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


#TODO: move to utils since it depends on numpy + scipy
#TODO: add fix_aspect_ratio option (use x_scale if Ly = 0, ...)
def read_image(mesh, filename, Lx = None, Ly = None, pos_x = None, pos_y = None):
    r"""
    Read image using pyvista and interpolate on given mesh

    :param :class:`Mesh`: Target mesh
    :param str filename:  Filename of the image
    :return :class:`torch.Tensor`: 2D tensor containing the correspoding image data

    :Examples:
      .. code::
        field = read_image(mesh, "measurement.png")
    """
    # read image data and convert to unique ids (0,1,2,...)
    image = pv.read(filename)
    data = image.get_array(image.array_names[0])
    if len(data.shape) == 2:
        data = np.prod(data, axis=1)
    data = np.unique(data, return_inverse=True)[1]
    data = data.reshape([image.dimensions[1], image.dimensions[0]]).T

    # scale and translate image
    if pos_x == None:
        pos_x = mesh.origin[0]
    if pos_y == None:
        pos_y = mesh.origin[1]
    if Lx == None:
        Lx = mesh.n[0] * mesh.dx[0]
    if Ly == None:
        Ly = mesh.n[1] * mesh.dx[1]

    x = np.linspace(0,Lx,image.dimensions[0]) + pos_x
    y = np.linspace(0,Ly,image.dimensions[1]) + pos_y
    xx, yy = np.meshgrid(x, y, indexing = "ij")
    xx = xx.reshape(-1)
    yy = yy.reshape(-1)
    data = data.reshape(-1)

    # interpolate on mesh
    x_mesh = np.arange(mesh.n[0]) * mesh.dx[0] + mesh.dx[0]/2. + mesh.origin[0]
    y_mesh = np.arange(mesh.n[1]) * mesh.dx[1] + mesh.dx[1]/2. + mesh.origin[1]
    xx_mesh, yy_mesh = np.meshgrid(x_mesh, y_mesh, indexing = "ij")

    return scipy.interpolate.griddata((xx, yy), data, (xx_mesh, yy_mesh))

def read_msh(mesh, filename):
    r"""
    Read unstructured msh meshes using pyvista

    :param str filename: Filename to be read
    :return :class:`Mesh` & dict: Mesh object and dictionary containing interpolated data tensors

    :Examples:
      .. code::
        fields = read_msh(mesh, "cylinder.msh")
    """
    unstructured_mesh = pv.read(filename)

    # read data from file
    p = unstructured_mesh.points
    fields = {}

    # interpolate files on mesh
    mesh.SpatialCoordinates()
