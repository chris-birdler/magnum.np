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
from .ode_solvers import RKF45
from .llg import LLGSolver
import torch
from xitorch.interpolate import Interp1D

__all__ = ["StringSolver"]

class StringSolver(object):
    def __init__(self, fields, num_images = 20, fix_start = False, fix_end = False, omega = 0., interpolation = "cspline"):
        self._num_images = num_images
        self._fix_start = fix_start
        self._fix_end = fix_end
        self._omega = omega
        self._interpolation = interpolation
        self._fields = fields
        self._llg = LLGSolver(self._fields)


    def E(self, state, images):
        E = []
        for image in images:
            state.m = image
            E.append(sum([field.E(state) for field in self._fields]))
        return E


    def step(self, state, images, h = 1e-11, level = True):
        """
        Make one Step steepest descent along the gradient for each image.
        """
        for i, image in enumerate(images):
            if self._fix_start and i == 0: continue
            if self._fix_end and i == len(images): continue
            state.t = 0.0
            state.m = image
            llg = LLGSolver(self._fields)
            llg.step(state, h)
            images[i] = state.m

        if level:
            newimages = []
            limages = self.level(state, images)
            for i in range(limages.shape[-1]):
                image = limages[:,:,:,:, i]
                newimages.append(image)
            images = newimages

        return images


    def level(self, state, images):
        logging.info("[StringSolver]: Leveling Images")
        Epath = self.E(state, images)
        source_data = torch.stack(images, dim = -1)
        Emin = min(Epath)
        Emax = max(Epath)

        x_source = [0.]
        for i in range(1, (len(images))):
            w = (((Epath[i] + Epath[i-1])/2. - Emin) / (Emax - Emin) + 1.)**self._omega
            dm = torch.linalg.norm(images[i] - images[i-1])
            x_source.append(w * dm + x_source[-1])
        x_source = state._tensor(x_source)
        x_target = state._linspace(x_source[0], x_source[-1], self._num_images)

        # start interpolation
        try:
            target_data = Interp1D(x_source, source_data, method = self._interpolation)(x_target)
        except:
            target_data = Interp1D(x_source, source_data, 'linear')(x_target)

        for i in range(self._num_images):
            state.m = target_data[:,:,:,:, i]
            state.m.normalize()
            target_data[:,:,:,:,i] = state.m

        return target_data
