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
import scipy
import pyvista as pv
import os
from . import Mesh, DecoratedTensor

__all__ = ["write_vtr", "write_vti", "read_vti", "read_image", "read_mesh"]

def write_vtr(fields, filename, state = None):
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

    x = torch.hstack([state._tensor([0.]), state.dx[0].cumsum(0)]).cpu().numpy() + state.mesh.origin[0]
    y = torch.hstack([state._tensor([0.]), state.dx[1].cumsum(0)]).cpu().numpy() + state.mesh.origin[1]
    z = torch.hstack([state._tensor([0.]), state.dx[2].cumsum(0)]).cpu().numpy() + state.mesh.origin[2]

    grid = pv.RectilinearGrid(x, y, z)

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


def write_vti(fields, filename, state = None):
    r"""
    Write vti files (equidistant rectangular grid, compressed) using pyvista.

    :param fields: single torch Tensor or List/Dictionary of tensors to be written
    :type fields: :class:`Tensor`, list, dict
    :param filename: filename to be writen
    :type filename: str
    :param state: filename to be writen
    :type state: :class:`State`

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
def read_image(mesh, filename, Lx = None, Ly = None, pos_x = None, pos_y = None, fix_aspect_ratio = False):
    r"""
    Read image using pyvista and interpolate on given mesh

    :param :class:`Mesh`: Target mesh
    :param str filename:  Filename of the image
    :param float Lx: Length to which the image should be scaled (defaults to the mesh length)
    :param float Ly: Height to which the image should be scaled (defaults to the mesh height)
    :param float pos_x: x-Offest by which the image should be shifted (defaults to the mesh origin)
    :param float pos_y: y-Offest by which the image should be shifted (defaults to the mesh origin)
    :param bool fix_aspect_ratio: if True only Lx or Ly can be set. The same scale will then be applied to both dimentions.
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

    if Lx != None and Ly != None and fix_aspect_ratio == True:
        raise RuntimeError("Aspect ratio cannot be kept fix, if both Lx and Ly are provided!")
    if Ly == None:
        Ly = mesh.n[1] * mesh.dx[1]
        if fix_aspect_ratio == True:
            Lx = Ly * image.dimensions[0] / image.dimensions[1]
    if Lx == None:
        Lx = mesh.n[0] * mesh.dx[0]
        if fix_aspect_ratio == True:
            Ly = Lx * image.dimensions[1] / image.dimensions[0]

    x_image = np.linspace(0,Lx,image.dimensions[0]) + pos_x
    y_image = np.linspace(0,Ly,image.dimensions[1]) + pos_y
    xx_image, yy_image = np.meshgrid(x_image, y_image, indexing = "ij")
    xx_image = xx_image.reshape(-1)
    yy_image = yy_image.reshape(-1)
    data = data.reshape(-1)

    # interpolate on mesh
    x = np.arange(mesh.n[0]) * mesh.dx[0] + mesh.dx[0]/2. + mesh.origin[0]
    y = np.arange(mesh.n[1]) * mesh.dx[1] + mesh.dx[1]/2. + mesh.origin[1]
    xx, yy = np.meshgrid(x, y, indexing = "ij")

    return scipy.interpolate.griddata((xx_image, yy_image), data, (xx, yy), fill_value=-1)



def read_mesh(mesh, filename, scale = 1.):
    r"""
    Read unstructured msh meshes using pyvista

    :param str filename: Filename to be read
    :return :class:`Mesh` & dict: Mesh object and dictionary containing interpolated data tensors

    :Examples:
      .. code::
        fields = read_msh(mesh, "cylinder.msh")
    """
    # read image data and volume domains
    unstructured_mesh = pv.read(filename)

    # interpolate on mesh
    x = np.arange(mesh.n[0]) * mesh.dx[0] + mesh.dx[0]/2. + mesh.origin[0]
    y = np.arange(mesh.n[1]) * mesh.dx[1] + mesh.dx[1]/2. + mesh.origin[1]
    z = np.arange(mesh.n[2]) * mesh.dx[2] + mesh.dx[2]/2. + mesh.origin[2]
    points = np.stack(np.meshgrid(x, y, z, indexing = "ij"), axis=-1).reshape(-1,3) / scale

    containing_cells = unstructured_mesh.find_containing_cell(points)
    data = unstructured_mesh.get_array(0)[containing_cells]
    data[containing_cells == -1] = -1 # containing_cell == -1, if point is not included in any cell

    return data.reshape(mesh.n)
