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

from magnumnp.common import timedmethod, constants
import torch
from .field_terms import LinearFieldTerm
import numpy as np

__all__ = ["DMIField", "InterfaceDMIField", "BulkDMIField", "D2dDMIField"]

class DMIField(LinearFieldTerm):
    r"""
    General Dzyaloshinskii-Moriya interaction

    The general expression for the DMI field can be expressed as

    .. math::
        \vec{h}^\text{dmi}(\vec{x}) = \frac{2 \, D}{\mu_0 \, M_s} \; \sum_{k=x,y,z} \vec{e}^\text{dmi}_k \times \frac{\partial \vec{m}}{\partial \vec{e}_k},

    with the DMI strength :math:`D` and the DMI vectors :math:`\vec{e}^\text{dmi}_k`, which describe which components of the gradient of :math:`\vec{m}` contribute to which component of the corresponding field. It is assumed that :math:`\vec{e}^\text{dmi}_{-k} = -\vec{e}^\text{dmi}_k`.

    The occuring gradient is discretized using central differences which finally yields

    .. math::
        \vec{h}^\text{dmi}_i = \frac{2 \, D_i}{\mu_0 \, M_{s,i}} \; \sum_{k=\pm x, \pm y,\pm z} \frac{\vec{e}^\text{dmi}_k \times \vec{m}_{i+\vec{e}_k}}{2 \, \Delta_k}.

    :param Ku: Name of the material parameter for the anisotropy constant :math:`K_\text{u}`, defaults to "Ku"
    :type Ku: str, optional
    :param Ku_axis: Name of the material parameter for the anisotropy axis :math:`\vec{e}_\text{u}`, defaults to "Ku_axis"
    :type Ku_axis: str, optional
    """
    parameters = ["D"]
    def __init__(self, dmi_vector, **kwargs):
        self._dmi_vector = dmi_vector
        super().__init__(**kwargs)

    @timedmethod
    def h(self, state):
        D = state.material[self.D].torch_tensor
        Ms = state.material["Ms"].torch_tensor
        m = state.m.torch_tensor
        h = self._h(m, D, Ms, state)
        return state.Tensor(h)

    @torch.compile
    def _h(self, m, D, Ms, state):
        h = state._zeros(state.mesh.n + (3,))
        # x
        v = state._tensor(self._dmi_vector[0]).expand(m[1:,:,:].shape)
        D_avg = torch.where(D[1:,:,:]*D[:-1,:,:] < 0,
                            torch.sqrt(torch.sqrt(-D[1:,:,:]*D[:-1,:,:])*torch.abs(D[1:,:,:]+D[:-1,:,:]) / 2.), # TODO: is sign missing?
                            2.*D[1:,:,:]*D[:-1,:,:]/(D[1:,:,:]+D[:-1,:,:]))
        h[:-1,:,:] += D_avg * torch.linalg.cross(v, m[ 1:,:,:]) / (2.*state.mesh.dx[0])
        h[ 1:,:,:] -= D_avg * torch.linalg.cross(v, m[:-1,:,:]) / (2.*state.mesh.dx[0])

        # y
        v = state._tensor(self._dmi_vector[1]).expand(m[:,1:,:].shape)
        D_avg = torch.where(D[:,1:,:]*D[:,:-1,:] < 0,
                            torch.sqrt(torch.sqrt(-D[:,1:,:]*D[:,:-1,:])*torch.abs(D[:,1:,:]+D[:,:-1,:]) / 2.),
                            2.*D[:,1:,:]*D[:,:-1,:]/(D[:,1:,:]+D[:,:-1,:]))
        h[:,:-1,:] += D_avg * torch.linalg.cross(v, m[:, 1:,:]) / (2.*state.mesh.dx[1])
        h[:, 1:,:] -= D_avg * torch.linalg.cross(v, m[:,:-1,:]) / (2.*state.mesh.dx[1])

        # z
        v = state._tensor(self._dmi_vector[2]).expand(m[:,:,1:].shape)
        D_avg = torch.where(D[:,:,1:]*D[:,:,:-1] < 0,
                            torch.sqrt(torch.sqrt(-D[:,:,1:]*D[:,:,:-1])*torch.abs(D[:,:,1:]+D[:,:,:-1]) / 2.),
                            2.*D[:,:,1:]*D[:,:,:-1]/(D[:,:,1:]+D[:,:,:-1]))
        h[:,:,:-1] += D_avg * torch.linalg.cross(v, m[:,:, 1:]) / (2.*state.mesh.dx[2])
        h[:,:, 1:] -= D_avg * torch.linalg.cross(v, m[:,:,:-1]) / (2.*state.mesh.dx[2])

        h *= 2. / (constants.mu_0 * Ms)
        h = torch.nan_to_num(h, posinf=0, neginf=0)
        return h

class InterfaceDMIField(DMIField):
    r"""
    Interface Dzyaloshinskii-Moriya interaction

    .. math::
        \vec{h}^\text{dmii}(\vec{x}) = -\frac{2 \, D_i}{\mu_0 \, M_s} \; \left[ \nabla \left(\vec{e}_z \cdot \vec{m} \right) - \left(\nabla \cdot \vec{m} \right) \, \vec{e}_z\right],

    with the DMI strength :math:`D_i` and an interface normal :math:`\vec{e}_z` in z-direction.
    The corresponding DMI vectors are :math:`\vec{e}^\text{dmi}_x = [ 0, 1, 0]`, :math:`\vec{e}^\text{dmi}_y = [-1, 0, 0]`, and :math:`\vec{e}^\text{dmi}_z = [0, 0, 0]`.

    :param Di: Name of the material parameter for the anisotropy constant :math:`D_i`, defaults to "Di"
    :type Di: str, optional
    """
    def __init__(self, Di = "Di"):
        dmi_vector = [[ 0, 1, 0], # x
                      [-1, 0, 0], # y
                      [ 0, 0, 0]] # z
        super().__init__(dmi_vector, D = Di)


class BulkDMIField(DMIField):
    r"""
    Bulk Dzyaloshinskii-Moriya interaction

    .. math::
        \vec{h}^\text{dmib}(\vec{x}) = -\frac{2 \, D_b}{\mu_0 \, M_s} \; \nabla \times \vec{m},

    with the DMI strength :math:`D_b`.
    The corresponding DMI vectors are :math:`\vec{e}^\text{dmi}_x = [1, 0, 0]`, :math:`\vec{e}^\text{dmi}_y = [0, 1, 0]`, and :math:`\vec{e}^\text{dmi}_z = [0, 0, 1]`.

    :param Db: Name of the material parameter for the anisotropy constant :math:`D_i`, defaults to "Di"
    :type Db: str, optional
    """
    def __init__(self, Db = "Db"):
        dmi_vector = [[1, 0, 0], # x
                      [0, 1, 0], # y
                      [0, 0, 1]] # z
        super().__init__(dmi_vector, D = Db)


class D2dDMIField(DMIField):
    r"""
    D2d Dzyaloshinskii-Moriya interaction

    .. math::
        \vec{h}^\text{dmib}(\vec{x}) = -\frac{2 \, D_b}{\mu_0 \, M_s} \; \nabla \times \vec{m},

    with the DMI strength :math:`D_{D2d}`.
    The corresponding DMI vectors are :math:`\vec{e}^\text{dmi}_x = [-1, 0, 0]`, :math:`\vec{e}^\text{dmi}_y = [0, 1, 0]`, and :math:`\vec{e}^\text{dmi}_z = [0, 0, 0]`.

    :param DD2d: Name of the material parameter for the anisotropy constant :math:`D_{D2d}`, defaults to "DD2d"
    :type DD2d: str, optional
    """
    def __init__(self, DD2d = "DD2d"):
        dmi_vector = [[-1, 0, 0], # x
                      [ 0, 1, 0], # y
                      [ 0, 0, 0]] # z
        super().__init__(dmi_vector, D = DD2d)
