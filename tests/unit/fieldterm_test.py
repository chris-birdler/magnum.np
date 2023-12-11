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
import pathlib
import torch
from magnumnp import *
from helpers import *


@pytest.mark.parametrize("field_term", [DemagField(), DemagFieldPBC(), InterfaceDMIField(), BulkDMIField(), D2dDMIField(), ExchangeField(), ExternalField([-24.6e-3/constants.mu_0, +4.3e-3/constants.mu_0, 0.0]), UniaxialAnisotropyField() ])
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

def test_regression():
    n  = (10, 20, 30)
    dx = (1e-9, 2e-9, 5e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"alpha":    state.Constant([0.02]),
                      "Ms":       state.Constant([8e5]),
                      "A":        state.Constant([1.3e-11]),
                      "Ku":       state.Constant([1e5]),
                      "Ku_axis":  state.Constant([0,1,0]),
                      "Kc1":      state.Constant([1e3]),
                      "Kc2":      state.Constant([1e4]),
                      "Kc_alpha": state.Constant([0.1]),
                      "Kc_beta":  state.Constant([0.2]),
                      "Kc_gamma": state.Constant([0.3]),
                      "Di":       state.Constant([1.]),
                      "Db":       state.Constant([1.]),
                      "DD2d":     state.Constant([1.])}

    x, y, z = state.SpatialCoordinate()
    state.m = torch.stack([x*y, y*z, z*x], dim=-1)
    demag        = DemagField()
    demag_pbc    = DemagFieldPBC()
    dmi_i        = InterfaceDMIField()
    dmi_b        = BulkDMIField()
    dmi_D2d      = D2dDMIField()
    exchange     = ExchangeField()
    aniso        = UniaxialAnisotropyField()
    aniso_cubic  = CubicAnisotropyField()

    m = state.m.cpu()
    h = exchange.h(state)
    h_demag        = demag.h(state).cpu()
    h_demag_pbc    = demag_pbc.h(state).cpu()
    h_dmi_i        = dmi_i.h(state).cpu()
    h_dmi_b        = dmi_b.h(state).cpu()
    h_dmi_D2d      = dmi_D2d.h(state).cpu()
    h_exchange     = exchange.h(state).cpu()
    h_aniso        = aniso.h(state).cpu()
    h_aniso_cubic  = aniso_cubic.h(state).cpu()

    m = state.m.cpu()

    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "ref" / "h_regression.vti"
    ## Uncomment to updated reference data
    #write_vti({"h_demag":h_demag,
    #           "h_demag_pbc":h_demag_pbc,
    #           "h_dmi_i":h_dmi_i,
    #           "h_dmi_b":h_dmi_b,
    #           "h_dmi_D2d":h_dmi_D2d,
    #           "h_exchange":h_exchange,
    #           "h_aniso":h_aniso,
    #           "h_aniso_cubic":h_aniso_cubic},
    #           filename)
    mesh, ref = read_vti(filename)

    torch.testing.assert_close(torch.linalg.cross(m, h_demag        / ref["h_demag"].max()),        torch.linalg.cross(m, ref["h_demag"]        / ref["h_demag"].max()),         atol=1e-15, rtol=1e-6)
    torch.testing.assert_close(torch.linalg.cross(m, h_demag_pbc    / ref["h_demag_pbc"].max()),    torch.linalg.cross(m, ref["h_demag_pbc"]    / ref["h_demag_pbc"].max()),     atol=1e-15, rtol=1e-6)
    torch.testing.assert_close(torch.linalg.cross(m, h_dmi_i        / ref["h_dmi_i"].max()),        torch.linalg.cross(m, ref["h_dmi_i"]        / ref["h_dmi_i"].max()),         atol=1e-15, rtol=1e-6)
    torch.testing.assert_close(torch.linalg.cross(m, h_dmi_b        / ref["h_dmi_b"].max()),        torch.linalg.cross(m, ref["h_dmi_b"]        / ref["h_dmi_b"].max()),         atol=1e-15, rtol=1e-6)
    torch.testing.assert_close(torch.linalg.cross(m, h_dmi_D2d      / ref["h_dmi_D2d"].max()),      torch.linalg.cross(m, ref["h_dmi_D2d"]      / ref["h_dmi_D2d"].max()),       atol=1e-15, rtol=1e-6)
    torch.testing.assert_close(torch.linalg.cross(m, h_exchange     / ref["h_exchange"].max()),     torch.linalg.cross(m, ref["h_exchange"]     / ref["h_exchange"].max()),      atol=1e-15, rtol=1e-6)
    torch.testing.assert_close(torch.linalg.cross(m, h_aniso        / ref["h_aniso"].max()),        torch.linalg.cross(m, ref["h_aniso"]        / ref["h_aniso"].max()),         atol=1e-15, rtol=1e-6)
    torch.testing.assert_close(torch.linalg.cross(m, h_aniso_cubic  / ref["h_aniso_cubic"].max()),  torch.linalg.cross(m, ref["h_aniso_cubic"]  / ref["h_aniso_cubic"].max()),   atol=1e-15, rtol=1e-6)

@pytest.mark.parametrize("field_term", [DemagFieldNonEquidistant(), ExchangeField(), ExternalField([-24.6e-3/constants.mu_0, +4.3e-3/constants.mu_0, 0.0]), UniaxialAnisotropyField() ])
def test_nonequidistant(field_term):
    n  = (10, 5, 4)
    dx2 = torch.ones(n[2]) * 5e-9
    dx2[2:] = 1.
    dx = (1e-9, 2e-9, dx2)

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

def test_energy_nonequidistant():
    n  = (10, 5, 4)
    dx2 = torch.ones(n[2]) * 5e-9
    dx2[2:] = 1e-9
    dx = (1e-9, 2e-9, dx2)

    mesh = Mesh(n, dx)
    state = State(mesh)

    state.material = {"Ms":      state.Constant([1.]),
                      "Ku":      state.Constant([1.]),
                      "Ku_axis": state.Constant([0,1,0])}

    aniso = UniaxialAnisotropyField()

    state.m = state.Constant([1,1,0])
    state.m[:,:,2:,:] = 0.
    E1 = aniso.E(state)

    state.m = state.Constant([1,1,0])
    state.m[:,:,:2,:] = 0.
    E2 = aniso.E(state)

    torch.testing.assert_close(E1, 5*E2, atol=0, rtol=1e-15)

def test_energy_domain():
    n  = (10, 5, 4)
    dx2 = torch.ones(n[2]) * 5e-9
    dx2[2:] = 1e-9
    dx = (1e-9, 2e-9, dx2)

    mesh = Mesh(n, dx)
    state = State(mesh)

    state.material = {"Ms":      state.Constant([1.]),
                      "Ku":      state.Constant([1.]),
                      "Ku_axis": state.Constant([0,1,0])}

    aniso = UniaxialAnisotropyField()

    state.m = state.Constant([1,1,0])
    domain1 = state.Constant(True, dtype=torch.bool)
    domain1[:,:,2:] = False

    E1 = aniso.E(state, domain1)
    E2 = aniso.E(state, ~domain1)

    torch.testing.assert_close(E1, 5*E2, atol=0, rtol=1e-15)

@pytest.mark.parametrize("field_term", [InterfaceDMIField(), BulkDMIField(), D2dDMIField(), ExchangeField()])
def test_pbc(field_term):
    n  = (20, 10, 5)
    dx = (1e-9, 1e-9, 1e-9)

    # no PBC
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

    x, y, z = state.SpatialCoordinate()
    state.m = torch.stack([x*y, y*z, z*x], dim=-1)
    h1 = field_term.h(state)

    # PBC
    mesh = Mesh(n, dx, pbc = (10,10,0))
    state = State(mesh)

    state.material = {"alpha":   state.Constant([0.02]),
                      "Ms":      state.Constant([8e5]),
                      "A":       state.Constant([1.3e-11]),
                      "Di":      state.Constant([1.]),
                      "Db":      state.Constant([1.]),
                      "DD2d":    state.Constant([1.])}

    x, y, z = state.SpatialCoordinate()
    state.m = torch.stack([x*y, y*z, z*x], dim=-1)
    h2 = field_term.h(state)
    diff = (h1-h2)[1:-1,1:-1,1:-1,:]

    torch.testing.assert_close(diff, torch.zeros_like(diff), atol=1e-20, rtol=1e-20)
