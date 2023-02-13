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

import pytest
import torch
import numpy as np
from magnumnp import *
from helpers import *


@pytest.mark.parametrize("field_term", [DemagField(), DemagFieldPBC(), InterfaceDMIField(), BulkDMIField(), D2dDMIField(), ExchangeField(), ExchangeFieldPBC(), ExternalField([-24.6e-3/constants.mu_0, +4.3e-3/constants.mu_0, 0.0]), UniaxialAnisotropyField() ])
def test_material_constant(field_term):
    n  = (10, 5, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)

    state.material = {"alpha":   state.Constant([0.02]),
                      "Ms":      state.Constant([8e5]),
                      "A":       state.Constant([1.3e-11]),
                      "Ku":      state.Constant([1e5]),
                      "Ku_axis": state.Constant([0,1,0]),
                      "Di":      state.Constant([1.]),
                      "Db":      state.Constant([1.]),
                      "DD2d":    state.Constant([1.])}

    state.m = state.Constant([1,0,0])
    state.m[5:,:,:,0] = -1.0

    h = field_term.h(state)


@pytest.mark.parametrize("field_term", [DemagField(), DemagFieldPBC(), InterfaceDMIField(), BulkDMIField(), D2dDMIField(), ExchangeField(), ExchangeFieldPBC(), ExternalField([-24.6e-3/constants.mu_0, +4.3e-3/constants.mu_0, 0.0]), UniaxialAnisotropyField() ])
def test_material_tensor(field_term):
    n  = (10, 5, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)

    state.material = {"alpha":   state.Tensor([0.02]),
                      "Ms":      state.Tensor([8e5]),
                      "A":       state.Tensor([1.3e-11]),
                      "Ku":      state.Tensor([1e5]),
                      "Ku_axis": state.Tensor([0,1,0]),
                      "Di":      state.Tensor([1.]),
                      "Db":      state.Tensor([1.]),
                      "DD2d":    state.Tensor([1.])}

    state.m = state.Constant([1,0,0])
    state.m[5:,:,:,0] = -1.0

    h = field_term.h(state)


@pytest.mark.parametrize("field_term", [DemagField(), DemagFieldPBC(), InterfaceDMIField(), BulkDMIField(), D2dDMIField(), ExchangeField(), ExchangeFieldPBC(), ExternalField([-24.6e-3/constants.mu_0, +4.3e-3/constants.mu_0, 0.0]), UniaxialAnisotropyField() ])
def test_material_tensor2(field_term):
    n  = (10, 5, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)

    state.material = {"alpha":   state.Tensor(0.02),
                      "Ms":      state.Tensor(8e5),
                      "A":       state.Tensor(1.3e-11),
                      "Ku":      state.Tensor(1e5),
                      "Ku_axis": state.Tensor([0,1,0]),
                      "Di":      state.Tensor(1.),
                      "Db":      state.Tensor(1.),
                      "DD2d":    state.Tensor(1.)}

    state.m = state.Constant([1,0,0])
    state.m[5:,:,:,0] = -1.0

    h = field_term.h(state)


@pytest.mark.parametrize("field_term", [DemagField(), DemagFieldPBC(), InterfaceDMIField(), BulkDMIField(), D2dDMIField(), ExchangeField(), ExchangeFieldPBC(), ExternalField([-24.6e-3/constants.mu_0, +4.3e-3/constants.mu_0, 0.0]), UniaxialAnisotropyField() ])
def test_material_float(field_term):
    n  = (10, 5, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)

    state.material = {"alpha":   0.02,
                      "Ms":      8e5,
                      "A":       1.3e-11,
                      "Ku":      1e5,
                      "Ku_axis": [0,1,0],
                      "Di":      1.,
                      "Db":      1.,
                      "DD2d":    1.}

    state.m = state.Constant([1,0,0])
    state.m[5:,:,:,0] = -1.0

    h = field_term.h(state)
