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

__all__ = ["BC", "PlaneBC", "Plane"]

class BC:
    r"""
    Boundary conditions defined on cell centers. Used for Dirichlet boundary conditions with the ``LLGWithLESolver``.

    Examples
    --------

    1. Standard usage:
        .. code:: python

            u_mask = torch.zeros(n, dtype=torch.bool)
            u_mask[-1] = 1
            
            u_vals = torch.zeros((n[1], n[2], 3))

            dirichlet_bcs = []
            dirichlet_bcs.append(BC(u_mask, u_vals))

            state.bcs = dict()
            state.bcs["ud"] = dirichlet_bcs

    2. With callable:
        .. code:: python

            u_mask = torch.zeros(n, dtype=torch.bool)
            u_mask[0] = 1
            
            def u_cond(state):
                n = state.mesh.n

                A = 1e-11
                omega = 2*np.pi*4e9
                u0 = torch.zeros((n[1], n[2], 3))
                u0[:,:,0] = A*np.cos(omega*float(state.t))
            
            return u0

            dirichlet_bcs = []
            dirichlet_bcs.append(BC(u_mask, u_cond))

            state.bcs = dict()
            state.bcs["ud"] = dirichlet_bcs

    Parameters
    ----------
    mask : :class:`torch.Tensor`
        boolean tensor that is used to set the displacement via ``state.ud[mask] = condition`` or ``state.ud[mask] = condition(state)``.
    condition : :class:`torch.Tensor` or callable
        Boundary values given by a tensor, or a function of state that returns a tensor. The shape of ``condition`` has to match ``state.ud[mask]``
    """

    def __init__(self, mask, condition):
        self._mask = mask
        
        if callable(condition):
            self._condition = condition
        elif isinstance(condition, torch.Tensor):
            self._condition = lambda state : condition
        else:
            raise Exception("BC condition has invalid type " + type(condition) + "!/n condition has to be either torch.Tensor or a function of the state that returns torch.Tensor.")

    @property
    def mask(self):
        return self._mask
    
    def condition(self, state):
        return self._condition(state)

class PlaneBC:
    r"""
    Boundary conditions defined on cell boundaries. Used for Neumann boundary conditions with the ``LLGWithLESolver``.

    Examples
    --------

    1. Standard usage:
        .. code:: python

            nx, ny, nz = state.mesh.n
            t_cond = torch.zeros((ny, nz, 3))

            neumann_bcs = []
            neumann_bcs.append(PlaneBC(Plane(0, 0, -1), t_cond))
            
            state.bcs = dict()
            state.bcs["t"] = neumann_bcs

    2. With callable:
        .. code:: python

            def t_cond(state, t0):
                nx, ny, nz = state.mesh.n

                t_in = torch.zeros((ny, nz, 3))
                t_in[:,:,2] = t0 * np.sin(float(state.t) * omega)

            return t_in

            neumann_bcs = []
            neumann_bcs.append(PlaneBC(Plane(0, 0, -1), lambda state : t_cond(state, t0)))
            
            state.bcs = dict()
            state.bcs["t"] = neumann_bcs

    Parameters
    ----------
    plane : :class:`Plane`
        The plane on which the boundary condition is defined.
    condition : :class:`torch.Tensor` or callable
        Boundary values given by a tensor, or a function of state that returns a tensor. The shape of ``condition`` needs to match the shape of ``plane``
    """

    def __init__(self, plane, condition):
        self.plane = plane

        if callable(condition):
            self._condition = condition
        elif isinstance(condition, torch.Tensor):
            self._condition = lambda state : condition
        else:
            raise Exception("BC condition has invalid type " + type(condition) + "!/n condition has to be either torch.Tensor or a function of the state that returns torch.Tensor.")

    def condition(self, state):
        return self._condition(state)
    
class Plane:
    r"""
    Defines a plane for ``PlaneBC``

    Examples
    --------

    1. A plane on the top :math:`z` boundary, with the normal pointing outwards:
        .. code:: python

            Plane(2, -1, 1)

    2. A plane normal to :math:`y` at :math:`y=\Delta y_0 + \Delta y_1`, on the first 3 layers along :math:`z`:
        .. code:: python

           Plane(1, 2, -1, trans_lim2=[0,4])

    Parameters
    ----------
    normal_dim : int
        The direction of the normal vector (0=x, 1=y, 2=z).
    normal_position_index : int
        The cell's positional index along the normal direction
    normal_sign : int 
        For ``i=normal_position_index``, ``-1`` marks the interface to the cell with index ``i-1`` , ``1`` the interface to the cell with index ``i+1``.
    trans_lim1 : list of int
        Bounds of the plane on the transversal direction with lowest index
    trans_lim2 : list of int
        Bounds of the plane on the transversal direction with highest index
    """

    def __init__(self, 
                 normal_dim, # (int) direction of the normal vector (0=x, 1=y, 2=z)
                 normal_position_index, # (int) positional index in the normal dimension
                 normal_sign, # (int), 1 or -1, direction in which the normal vector points, points away from the interiour
                 trans_lim1=[0,None], # bounds of the transversal direction with lowest index
                 trans_lim2=[0,None]): # bounds of the transversal direction with the higgest index
        
        assert isinstance(normal_dim,int)
        self.dim = normal_dim

        assert isinstance(normal_position_index,int)
        self.pos = normal_position_index

        assert isinstance(normal_sign, int)
        if (self.pos == 0):
            assert normal_sign == -1
        if (self.pos == -1):
            assert normal_sign == 1
            
        self.sign = normal_sign

        assert isinstance(trans_lim1, list) or isinstance(trans_lim1, tuple) 
        assert isinstance(trans_lim1[0], int) or isinstance(trans_lim1[0], None)
        self.trans_slice1 = slice(trans_lim1[0], trans_lim1[1])

        assert isinstance(trans_lim2, list) or isinstance(trans_lim2, tuple) 
        assert isinstance(trans_lim2[0], int) or isinstance(trans_lim2[0], None)
        self.trans_slice2 = slice(trans_lim2[0], trans_lim2[1])

        if self.dim == 0:
            self.trans_dim1 = 1
            self.trans_dim2 = 2
        elif self.dim == 1:
            self.trans_dim1 = 0
            self.trans_dim2 = 2
        elif self.dim == 2:
            self.trans_dim1 = 0
            self.trans_dim2 = 1

        self.slices = [None,None,None]
        self.slices[self.dim] = self.pos
        self.slices[self.trans_dim1] = self.trans_slice1
        self.slices[self.trans_dim2] = self.trans_slice2
        self.slices = tuple(self.slices)

    def get_mask(self, state):
        mask = torch.zeros(state.mesh.n, dtype=torch.bool)
        mask[self.slices] = True
        return mask