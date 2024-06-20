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

from magnumnp.common import logging, timedmethod, constants
import torch

__all__ = ["SD_solver"]

class SD_solver(object):
    def __init__(self, terms, dm_max = 1e2, samples = 10):
        """
        This class implements the direct energy minimizing algorithm introduced in [Exl2014]_.

        .. note:: This feature is experimental.

        *Example*
          .. code:: python

            state = State(mesh)
            minimizer = Minimizer_BB([ExchangeField()])
            minimizer.minimize(state)

        *Arguments*
          terms ([:class:`LLGTerm`])
            List of LLG contributions to be considered for energy minimization
          region (:class:`str`)
            region on which the energy is minimized
          tau_min (:class:`float`)
            minimum step size
          tau_max (:class:`float`)
            maximum step size
          dm_max (:class:`float`)
            stop criterion given as supremum norm of dm/dt
          sample (:class:`int`)
            number of subsequent steps the stop criterion has to be fulfilled
        """
        self._terms = terms
        self._dm_max = dm_max
        self._samples = samples

    def _dm(self, state):
        h = sum([term.h(state) for term in self._terms])
        return torch.linalg.cross(state.m, torch.linalg.cross(state.m, h))
    
    def minimize(self, state):
        step = 0
        dm_max = 1e18
        last_dm_max = []
        energy = []
        steps = []
        tau = 1e-7

        while len(last_dm_max) < self._samples or max(last_dm_max) > self._dm_max:
            h = sum([term.h(state) for term in self._terms])
            dm = torch.linalg.cross(state.m, torch.linalg.cross(state.m, h))

            m_next = state.m - tau*dm

            # update state
            state.m = m_next
            
            dm_max = dm.max()
            last_dm_max.append(dm_max)
            if len(last_dm_max) > self._samples: last_dm_max.pop(0)

            E = sum([term.E(state) for term in self._terms])
            
            logging.info_blue("Tau: %.5g, dm_max: %.5g, E: %.5g" % (tau, dm_max, E))

            # increase step count
            step += 1
            energy.append(E)
            steps.append(step)
        return state, energy, steps

    

