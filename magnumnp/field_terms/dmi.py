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
    @torch.compile
    def h(self, state):
        D = state.material[self.D]
        Ms = state.material["Ms"]
        m = state.m
        dx = state.mesh.dx_tensor[0].reshape(-1,1,1,1)
        dy = state.mesh.dx_tensor[1].reshape(1,-1,1,1)
        dz = state.mesh.dx_tensor[2].reshape(1,1,-1,1)
        h = torch.zeros_like(state.m)

        # x
        if state.mesh.pbc[0] == 0:
            v = torch.tensor(self._dmi_vector[0]).expand(m[1:,:,:].shape)
            # TODO: implement for different signs (according to https://github.com/mumax/3/issues/236)
            # D_avg = torch.where(D[1:,:,:]*D[:-1,:,:] < 0, ...
            D_avg = 2.*D[1:,:,:]*D[:-1,:,:] / (D[1:,:,:]*dx[:-1,:,:,:] + D[:-1,:,:]*dx[1:,:,:,:])
            h[:-1,:,:] += D_avg * torch.linalg.cross(v, m[ 1:,:,:]) / 2.
            h[ 1:,:,:] -= D_avg * torch.linalg.cross(v, m[:-1,:,:]) / 2.
        else:
            v = torch.tensor(self._dmi_vector[0]).expand(m.shape)
            D_next = torch.roll(D, +1, dims=0)
            dx_next = torch.roll(dx, +1, dims=0)
            # TODO: implement for different signs (according to https://github.com/mumax/3/issues/236)
            # D_avg = torch.where(D_next*D < 0,
            D_avg = 2.*D_next*D / (D_next*dx + D*dx_next)
            h -= D_avg * torch.linalg.cross(v, torch.roll(state.m, +1, dims=0)) / 2.

            D_avg = torch.roll(D_avg, -1, dims=0)
            h += D_avg * torch.linalg.cross(v, torch.roll(state.m, -1, dims=0)) / 2.

        # y
        if state.mesh.pbc[1] == 0:
            v = torch.tensor(self._dmi_vector[1]).expand(m[:,1:,:].shape)
            # TODO: implement for different signs (according to https://github.com/mumax/3/issues/236)
            # D_avg = torch.where(D[:,1:,:]*D[:,:-1,:] < 0,
            D_avg = 2.*D[:,1:,:]*D[:,:-1,:] / (D[:,1:,:]*dy[:,:-1,:,:] + D[:,:-1,:]*dy[:,1:,:,:])
            h[:,:-1,:] += D_avg * torch.linalg.cross(v, m[:, 1:,:]) / 2.
            h[:, 1:,:] -= D_avg * torch.linalg.cross(v, m[:,:-1,:]) / 2.
        else:
            v = torch.tensor(self._dmi_vector[1]).expand(m.shape)
            D_next = torch.roll(D, +1, dims=1)
            dy_next = torch.roll(dx, +1, dims=1)
            # TODO: implement for different signs (according to https://github.com/mumax/3/issues/236)
            # D_avg = torch.where(D_next*D < 0,
            D_avg = 2.*D_next*D / (D_next*dy + D*dy_next)
            h -= D_avg * torch.linalg.cross(v, torch.roll(state.m, +1, dims=1)) / 2.

            D_avg = torch.roll(D_avg, -1, dims=1)
            h += D_avg * torch.linalg.cross(v, torch.roll(state.m, -1, dims=1)) / 2.

        # z
        if state.mesh.pbc[2] == 0:
            v = torch.tensor(self._dmi_vector[2]).expand(m[:,:,1:].shape)
            # TODO: implement for different signs (according to https://github.com/mumax/3/issues/236)
            # D_avg = torch.where(D[:,:,1:]*D[:,:,:-1] < 0,
            D_avg = 2.*D[:,:,1:]*D[:,:,:-1] / (D[:,:,1:]*dz[:,:,:-1,:] + D[:,:,:-1]*dz[:,:,1:,:])
            h[:,:,:-1] += D_avg * torch.linalg.cross(v, m[:,:, 1:]) / 2.
            h[:,:, 1:] -= D_avg * torch.linalg.cross(v, m[:,:,:-1]) / 2.
        else:
            v = torch.tensor(self._dmi_vector[2]).expand(m.shape)
            D_next = torch.roll(D, +1, dims=2)
            dz_next = torch.roll(dx, +1, dims=2)
            # TODO: implement for different signs (according to https://github.com/mumax/3/issues/236)
            # D_avg = torch.where(D_next*D < 0,
            D_avg = 2.*D_next*D / (D_next*dz + D*dz_next)
            h -= D_avg * torch.linalg.cross(v, torch.roll(state.m, +1, dims=2)) / 2.

            D_avg = torch.roll(D_avg, -1, dims=2)
            h += D_avg * torch.linalg.cross(v, torch.roll(state.m, -1, dims=2)) / 2.

        h *= 2. / (constants.mu_0 * Ms)
        return h.nan_to_num(posinf=0, neginf=0)



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
        dmi_vector = [[ 0., 1., 0.], # x
                      [-1., 0., 0.], # y
                      [ 0., 0., 0.]] # z
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
        dmi_vector = [[1., 0., 0.], # x
                      [0., 1., 0.], # y
                      [0., 0., 1.]] # z
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
        dmi_vector = [[-1., 0., 0.], # x
                      [ 0., 1., 0.], # y
                      [ 0., 0., 0.]] # z
        super().__init__(dmi_vector, D = DD2d)
