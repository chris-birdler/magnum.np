import pytest
import pathlib
import torch
import numpy as np
from magnumnp import *

def test_cube():
    n  = (100, 5, 5)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {"sigma": 1.}

    dirichlet_bc_nodes = state.Constant(False, dtype=bool)
    dirichlet_bc_nodes[0,:,:] = True
    dirichlet_bc_nodes[-1,:,:] = True

    state.u = state.Constant(0.)
    state.u[-1,:,:] = 1.

    ohm = OhmSolver(dirichlet_bc_nodes)
    u = ohm.u(state)
    j = ohm.j(state)
    write_vti({"u":u, "j":j}, "data/u.vti")

def test_curved():
    n  = (120, 60, 5)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx, origin = (-n[0]*dx[0]/2.,0.,0.))
    state = State(mesh)
    x, y, z = state.SpatialCoordinate()
    r_i = 30e-9
    r_a = 50e-9
    state.material["sigma"] = 5.3e5 * ((x**2+y**2 > r_i**2) & (x**2+y**2 < r_a**2))
    write_vti(state.material, "data/material.vti")

    dirichlet_bc_nodes = state.Constant(False, dtype=bool)
    dirichlet_bc_nodes[:,0,:] = True

    state.u = state.Constant(0.)
    state.u[:n[0]//2,0,:] = 1.

    ohm = OhmSolver(dirichlet_bc_nodes)
    u = ohm.u(state)
    j = ohm.j(state)
    write_vti({"u":u, "j":j}, "data/u.vti")

def test_curved_2D():
    n  = (120, 60, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx, origin = (-n[0]*dx[0]/2.,0.,0.))
    state = State(mesh)
    x, y, z = state.SpatialCoordinate()
    r_i = 30e-9
    r_a = 50e-9
    state.material["sigma"] = 1. * ((x**2+y**2 > r_i**2) & (x**2+y**2 < r_a**2))
    write_vti(state.material, "data/material.vti")

    dirichlet_bc_nodes = state.Constant(False, dtype=bool)
    dirichlet_bc_nodes[:,0,:] = True

    state.u = state.Constant(0.)
    state.u[:n[0]//2,0,:] = 1.

    ohm = OhmSolver(dirichlet_bc_nodes)
    u = ohm.u(state)
    j = ohm.j(state)
    write_vti({"u":u, "j":j}, "data/u.vti")

def test_sabri():
    n  = (1000, 1000, 1)
    dx = (5e-9, 5e-9, 1e-9)
    mesh = Mesh(n, dx, origin = (-n[0]*dx[0]/2.,-n[1]*dx[1]/2.,-n[2]*dx[2]/2.))
    state = State(mesh)
    x, y, z = state.SpatialCoordinate()

    domain = (x**2 < 1.0e-6**2) | (y**2 < 1.0e-6**2)
    state.material["sigma"] = 0.
    state.material["sigma"][domain] = 5.3e5

    write_vti(state.material, "data/material.vti")

    dirichlet_bc_nodes = state.Constant(False, dtype=bool)
    dirichlet_bc_nodes[0,:,:] = True
    dirichlet_bc_nodes[:,0,:] = True

    state.u = state.Constant(0.)
    state.u[0,:,:] = 1.

    ohm = OhmSolver(dirichlet_bc_nodes)
    u = ohm.u(state, max_iter = 1000)
    j = ohm.j(state, max_iter = 1000)
    write_vti({"u":u, "j":j}, "data/u.vti")
