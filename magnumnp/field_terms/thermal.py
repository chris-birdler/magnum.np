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
from .field_terms import FieldTerm

__all__ = ["ThermalField"]

class ThermalField(FieldTerm):
    r"""
    """
    parameters = ["T"]
    def __init__(self, domain=None, **kwargs):
        self._step = None
        self._alpha_ref = None
        self._Ms_ref = None
        super().__init__(**kwargs)

    @timedmethod
    #@torch.compile
    def h(self, state):
        if state._step != self._step: # update random field
            self._sigma = torch.normal(0., 1., size = state.m.shape, device = state.m.device, dtype = state.m.dtype)
            self._step = state._step

        # cache the time-invariant prefactor; constant material parameters return
        # the identical tensor object on every access, so an identity check is
        # sufficient (state-dependent parameters invalidate on every call)
        alpha = state.material["alpha"]
        Ms = state.material["Ms"]
        if self._alpha_ref is not alpha or self._Ms_ref is not Ms:
            self._pref = torch.sqrt(2. * alpha * constants.kb / (constants.mu_0 * Ms * constants.gamma * state.mesh.cell_volumes))
            self._alpha_ref, self._Ms_ref = alpha, Ms

        T = torch.as_tensor(state.T, dtype=state.m.dtype, device=state.m.device)
        h = self._sigma * self._pref * (T / state._dt)**0.5
        return h.nan_to_num_(posinf=0, neginf=0)
