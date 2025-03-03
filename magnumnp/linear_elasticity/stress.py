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

__all__ = ["sigma"]

@torch.compile
def sigma(state, eps):
    C = state.material["C"]

    #n = state.mesh.n
    #sig = state.zeros((n[0], n[1], n[2], 6))
    #for i in range(6):
    #    for j in range(6):
    #        sig[:,:,:,i] += C[:,:,:,i,j]*eps[:,:,:,j]

    sig = (C @ eps.unsqueeze(-1)).squeeze(-1)

    return sig