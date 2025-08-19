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
from magnumnp.field_terms.field_terms import LinearFieldTerm
from magnumnp.linear_elasticity import epsilon, epsilon_m, sigma
import torch

__all__ = ["LinearMagnetoElasticField", "MagnetoElasticField"]

class MagnetoElasticField(object):
    r"""
    Magnetoelastic Field

    This field term is obtained from the magnetoelastic energy:

    .. math::
        E = \int \left( \frac{1}{2}\epsilon^m :C:\epsilon^m - \epsilon :C:\epsilon^m \right) \text{d}\mathbf{x}

    References:
    - C.Y. Liang et al., Nanotechnology 25 (2014) 435701 (10pp), doi:10.1088/0957-4484/25/43/435701
    - Y.C. Shu et al., Mechanics of Materials 36 (2004) 975-997, doi:10.1016/j.mechmat.2003.04.004

    :param ud: Displacement. If neither this nor 'mechanical_strain' is set, `state.ud` is used.
    :type ud: :class:`torch.Tensor` or function, optional
    :param mechanical_strain: Strain in Voigt notation (mesh.n + (6,))
    :type mechanical_strain: :class:`torch.Tensor` or function, optional

    The parameters 'ud' and 'mechanical_strain' are mutually exclusive.
    """

    def __init__(self, ud = None, mechanical_strain = None, **kwargs):
        if (ud != None) and (mechanical_strain != None):
            raise Exception("[MagnetoElasticField] received arguments for both 'ud' and 'mechanical_strain', but these are mutually exclusive.")

        if ud != None:
            if callable(ud):
                self._get_ud = ud
            else:
                self._get_ud = lambda state : ud
        else:
            self._get_ud = self._default_get_ud

        if mechanical_strain != None:
            if callable(mechanical_strain):
                self._get_strain = mechanical_strain 
            else:
                self._get_strain = lambda state : mechanical_strain
        else:
            self._get_strain = self._default_get_strain

    def _default_get_ud(self, state):
        return state.ud
    
    def _default_get_strain(self, state):
        ud = self._get_ud(state)
        return epsilon(state, ud)

    @timedmethod
    @torch.compile
    def h(self, state):     
        Ms = state.material["Ms"]
        lambda_100 = state.material["lambda_100"][:,:,:,0]
        lambda_111 = state.material["lambda_111"][:,:,:,0]

        m = state.m
   
        eps = self._get_strain(state)
        eps_m = epsilon_m(state, m)
        sig_el = sigma(state, eps-eps_m)

        # factor 1/2 is dropped due to second derivative
        dx_eps_xx = -3.*lambda_100 * m[...,0]
        dy_eps_yy = -3.*lambda_100 * m[...,1]
        dz_eps_zz = -3.*lambda_100 * m[...,2]

        # factor 1/2 is already dropped by conversion to Voigt
        d_eps_x = -3.*lambda_111 * m[...,0]
        d_eps_y = -3.*lambda_111 * m[...,1]
        d_eps_z = -3.*lambda_111 * m[...,2]

        # factor 1/2 of the energy is droped due to product rule
        h = torch.zeros(state.mesh.n+(3,))
        h[:,:,:,0] = sig_el[:,:,:,0]*dx_eps_xx + sig_el[:,:,:,5]*d_eps_y + sig_el[:,:,:,4]*d_eps_z 
        h[:,:,:,1] = sig_el[:,:,:,1]*dy_eps_yy + sig_el[:,:,:,5]*d_eps_x + sig_el[:,:,:,3]*d_eps_z
        h[:,:,:,2] = sig_el[:,:,:,2]*dz_eps_zz + sig_el[:,:,:,4]*d_eps_x + sig_el[:,:,:,3]*d_eps_y

        h *= -1. / (constants.mu_0 * Ms)
        return h.nan_to_num(posinf=0, neginf=0)
    
    @torch.compile
    def E(self, state, domain = Ellipsis):
        m = state.m

        eps_m = epsilon_m(state, m)
        sig_m = sigma(state, eps_m)
        eps = self._get_strain(state)

        zeta = (0.5*eps_m - eps)*sig_m

        E =  zeta * state.mesh.cell_volumes
        return E[domain].sum()
    
class LinearMagnetoElasticField(LinearFieldTerm):
    r"""
    Magnetoelastic Field

    This field term is obtained only from the terms in the magnetoelastic energy
    that are quadratic in :math:`\vec{m}`:

    .. math::
        E = -\int \epsilon :C:\epsilon^m \text{d}\mathbf{x}

    :param ud: Displacement. If neither this nor 'mechanical_strain' is set, `state.ud` is used.
    :type ud: :class:`torch.Tensor` or function, optional
    :param mechanical_strain: Strain in Voigt notation (mesh.n + (6,))
    :type mechanical_strain: :class:`torch.Tensor` or function, optional
    """

    def __init__(self, ud = None, mechanical_strain = None, **kwargs):
        if (ud != None) and (mechanical_strain != None):
            raise Exception("[LinearMagnetoElasticField] received arguments for both 'ud' and 'mechanical_strain', but these are mutually exclusive.")

        if ud != None:
            if callable(ud):
                self._get_ud = ud
            else:
                self._get_ud = lambda state : ud
        else:
            self._get_ud = self._default_get_ud

        if mechanical_strain != None:
            if callable(mechanical_strain):
                self._get_strain = mechanical_strain 
            else:
                self._get_strain = lambda state : mechanical_strain
        else:
            self._get_strain = self._default_get_strain

    def _default_get_ud(self, state):
        return state.ud
    
    def _default_get_strain(self, state):
        ud = self._get_ud(state)
        return epsilon(state, ud)

    @timedmethod
    @torch.compile
    def h(self, state):     
        Ms = state.material["Ms"]
        lambda_100 = state.material["lambda_100"][...,0]
        lambda_111 = state.material["lambda_111"][...,0]
        C = state.material["C"]

        B1 = -3.*lambda_100*(C[...,0,0]-C[...,0,1]) / 2.
        B2 = -3.*lambda_111*C[...,3,3]

        m = state.m
        eps = self._get_strain(state)

        h = torch.zeros(state.mesh.n+(3,))
        # Note: eps[...,3:] already includes a factor of 2 due to Voigth notation
        h[...,0] = 2*B1*m[...,0]*eps[...,0] + B2*(eps[...,5]*m[...,1] + eps[...,4]*m[...,2])
        h[...,1] = 2*B1*m[...,1]*eps[...,1] + B2*(eps[...,5]*m[...,0] + eps[...,3]*m[...,2])
        h[...,2] = 2*B1*m[...,2]*eps[...,2] + B2*(eps[...,4]*m[...,0] + eps[...,3]*m[...,1])

        h *= -1. / (constants.mu_0 * Ms)
        return h.nan_to_num(posinf=0, neginf=0)
