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

def nsk(x, axis = 2):
    m = self.mean(axis=axis)
    dxm = torch.stack(torch.gradient(m, spacing = dx[0], dim = 0), dim = -1).squeeze(-1)
    dym = torch.stack(torch.gradient(m, spacing = dx[1], dim = 1), dim = -1).squeeze(-1)
    nsk = 1./(4. * pi) * (m * torch.cross(dxm, dym)).sum()* dx[0] * dx[1]
    return float(nsk.detach().cpu().numpy())

