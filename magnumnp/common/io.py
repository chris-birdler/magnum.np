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
from . import Mesh
from magnumnp.common import logging, Material

__all__ = ["write_vtr", "write_vti", "read_vti", "read_image", "read_mesh", "read_neper"]

def write_vtr(fields, filename, state = None, scale = 1.):
    if not filename.endswith(".vtr"):
        logging.warning("[write_vtr] Extention '.vtr' should be used on non-equidistant grids! (filename = '%s')" % filename)

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

    x = torch.hstack([torch.tensor([0.]), state.mesh.dx_tensor[0].cumsum(0)]).cpu().numpy() + state.mesh.origin[0]
    y = torch.hstack([torch.tensor([0.]), state.mesh.dx_tensor[1].cumsum(0)]).cpu().numpy() + state.mesh.origin[1]
    z = torch.hstack([torch.tensor([0.]), state.mesh.dx_tensor[2].cumsum(0)]).cpu().numpy() + state.mesh.origin[2]

    grid = pv.RectilinearGrid(x*scale, y*scale, z*scale)

    for name in fields:
        f = fields[name]
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


def write_vti(fields, filename, state = None, scale = 1.):
    r"""
    Write vti files (equidistant rectangular grid, compressed) using pyvista.

    :param fields: single torch Tensor or List/Dictionary of tensors to be written
    :type fields: :class:`Tensor`, list, dict
    :param filename: filename to be written
    :type filename: str
    :param state: filename to be written
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
    if str(filename)[-4:] != ".vti":
        logging.warning("[write_vti] Extention '.vti' should be used on equidistant grids!")

    dirname = os.path.dirname(filename)
    if dirname and not os.path.isdir(dirname):
        os.makedirs(dirname)

    if not (isinstance(fields, list) or isinstance(fields, dict)):
        fields = [fields]
    if not isinstance(fields, dict):
        fields = {"f%03d"%i:f for (i,f) in enumerate(fields)}

    if state is None:
        values = [fields[name] for name in fields]
        n = values[0].shape[:3]
        dx = (1., 1., 1.)
        origin = (0., 0., 0.)
    else:
        n = state.mesh.n
        dx = state.mesh.dx
        origin = state.mesh.origin

    grid = pv.ImageData(dimensions = np.array(n) + 1,
                        spacing = np.array(dx) * scale,
                        origin = np.array(origin) * scale)

    for name in fields:
        f = fields[name]
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


def read_vti(filename, scale = 1.):
    r"""
    Read vti files using pyvista

    :param str filename: Filename to be read
    :param float scale: scale with which the file was written
    :return :class:`Mesh` & dict: Mesh object and dictionary containing all data tensors


    :Examples:
      .. code::
        mesh, fields = read_vti("m0.vti")
    """
    fields = {}
    data = pv.read(filename)

    mesh = Mesh(np.array(data.dimensions)-1, np.array(data.spacing) / scale, np.array(data.origin) / scale)

    for name in data.array_names:
        f = data.get_array(name)
        vals = data.get_array(name)
        if len(vals.shape) == 1:
            dim = mesh.n
        else:
            dim = mesh.n + (vals.shape[-1],)
        t = torch.tensor(1.)
        f = torch.from_numpy(vals.reshape(dim, order="F")).to(device=t.device, dtype=t.dtype)
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

    data = scipy.interpolate.griddata((xx_image, yy_image), data, (xx, yy), fill_value=-1)
    return torch.tensor(data)



def read_mesh(mesh, filename, scale = 1.):
    r"""
    Read unstructured msh meshes using pyvista

    :param str filename: Filename to be read
    :return :class:`Mesh` & dict: Mesh object and dictionary containing interpolated data tensors

    :Examples:
      .. code::
        fields = read_mesh(mesh, "cylinder.msh")
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


def read_neper(filename, scale = 1.):
    r"""
    Read a rasterized Neper tessellation (.tesr) file.

    This function parses a Neper TESR file containing a rectilinear voxel
    representation of a tessellation. It extracts:

      • the grid dimensions and voxel spacings,
      • the per-cell (grain) IDs assigned to each voxel,
      • the per-cell group IDs (if present in the header).

    TESR files may contain both *id* and *group* lists in the header.
    The voxel data block stores only cell IDs, so group IDs are mapped
    voxel-wise using the cell → group association defined in the header.

    Parameters
    ----------
    filename : str
        Path to the `.tesr` file to be read.

    scale : float, optional
        Additional scaling factor applied to the voxel spacing
        (dx, dy, dz). Defaults to 1.

    Returns
    -------
    mesh : :class:`Mesh`
        A Mesh object describing grid resolution and voxel spacing.

    domains : torch.Tensor
        Integer tensor of shape (Nx, Ny, Nz) containing the **grain IDs**
        (cell IDs) for each voxel.

    groups : torch.Tensor
        Integer tensor of shape (Nx, Ny, Nz) containing the **group IDs**
        for each voxel. If the TESR file does not define groups, `groups`
        will contain zeros.

    Examples
    --------
    Generating a suitable Neper tessellation:

    .. code::

        # Create 1000 grains in a 500. x 250. x 40. domain
        # Assign group 1 to grains with id < 500, and group 2 otherwise
        # Rasterize to 500 x 125 x 10 voxels
        neper -T -n 1000 \
              -domain "cube(500,250,40)" \
              -group "id<500?1:2" \
              -o test \
              -format tess,tesr \
              -tesrsize 500:125:10

    Reading the generated TESR file:

    .. code::

        mesh, grain_ids, group_ids = read_neper("test.tesr")
    """
    Nx = Ny = Nz = None
    start_offset = None
    cell_ids = []
    group_ids = []

    with open(filename, "rb") as f:
        # --- Read ASCII header ---
        while True:
            pos = f.tell()
            line = f.readline().decode("ascii", errors="ignore")
            if not line:
                raise ValueError("Unexpected end of file while reading header.")

            line_s = line.strip()

            # --- Basic grid information ---
            if line_s.startswith("**general"):
                dim = int(f.readline().decode().strip())
                Nx, Ny, Nz = map(int, f.readline().decode().split())
                dx, dy, dz = tuple(scale * float(d) for d in f.readline().decode().split())

            # --- Cell definitions ---
            if line_s == "*id":
                # Next line(s) contain the cell IDs
                # They may span multiple lines until next "*"
                while True:
                    pos2 = f.tell()
                    l2 = f.readline().decode().strip()
                    if l2.startswith("*") or l2.startswith("**"):
                        f.seek(pos2)
                        break
                    cell_ids.extend(map(int, l2.split()))

            if line_s == "*group":
                # Next line(s) contain group IDs in same order as cell IDs
                while True:
                    pos2 = f.tell()
                    l2 = f.readline().decode().strip()
                    if l2.startswith("*") or l2.startswith("**"):
                        f.seek(pos2)
                        break
                    group_ids.extend(map(int, l2.split()))

            # --- Binary voxel block starts after **data / binary16 ---
            if line_s.startswith("**data"):
                fmt = f.readline().decode().strip()
                if fmt != "binary16":
                    raise ValueError(f"Unsupported TESR binary format: {fmt}")
                start_offset = f.tell()
                break

        if Nx is None:
            raise ValueError("TESR header did not contain voxel size.")

        # --- Read voxel cell ID block ---
        count = Nx * Ny * Nz
        vox_cell = np.fromfile(f, dtype="<u2", count=count)
        vox_cell = vox_cell.reshape((Nx, Ny, Nz), order="F")

    # --- Convert cell-id to group-id ---
    if group_ids:
        # cell_ids is 1-based, array index is 0-based
        lookup = np.zeros(max(cell_ids) + 1, dtype=np.int32)
        for cid, gid in zip(cell_ids, group_ids):
            lookup[cid] = gid
        vox_group = lookup[vox_cell]
    else:
        vox_group = None

    mesh = Mesh((Nx, Ny, Nz), (dx, dy, dz))
    domains = torch.tensor(vox_cell.astype(np.int64))
    groups = torch.tensor(vox_group.astype(np.int64))
    return mesh, domains, groups
