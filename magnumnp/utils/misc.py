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

__all__ = ["add_noise", "nsk"]

def add_noise(x, dev = 1.0, mean = 0.0):
   if torch.is_tensor(x):
        x += torch.empty_like(x).normal_(mean = mean, std = dev)
        x.normalize()

def nsk(state): # TODO: document and improve interface
    m = state.m.mean(axis=2)
    dxm = torch.stack(torch.gradient(m, spacing = state.mesh.dx[0], dim = 0), dim = -1).squeeze(-1)
    dym = torch.stack(torch.gradient(m, spacing = state.mesh.dx[1], dim = 1), dim = -1).squeeze(-1)
    return 1./(4.*pi) * (m * torch.linalg.cross(dxm, dym)).sum() * state.mesh.dx[0] * state.mesh.dx[1]


def hsl_to_rgb(h, s, l): # TODO: document and improve interface
    def hue_to_rgb(p, q, t):
        t = torch.where(t < 0, t + 1, t)
        t = torch.where(t > 1, t - 1, t)
        return torch.where(t < 1/6, p + (q - p) * 6 * t,
               torch.where(t < 1/2, q,
               torch.where(t < 2/3, p + (q - p) * (2/3 - t) * 6, p)))

    s = s.clamp(0, 1)
    l = l.clamp(0, 1)
    h = h.clamp(0, 1)

    q = l + s * torch.where(l < 0.5, l, 1 - l)
    p = 2 * l - q
    r = hue_to_rgb(p, q, h + 1/3)
    g = hue_to_rgb(p, q, h)
    b = hue_to_rgb(p, q, h - 1/3)

    return torch.stack([r, g, b], dim=-1)
