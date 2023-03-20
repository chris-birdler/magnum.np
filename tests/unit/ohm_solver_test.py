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

    dirichlet_bc_u0 = state._zeros([n+2 for n in state.mesh.n])
    dirichlet_bc_nodes = state._zeros([n+2 for n in state.mesh.n], dtype=bool)
    dirichlet_bc_nodes[0,:,:] = True
    dirichlet_bc_nodes[-1,:,:] = True

    dirichlet_bc_u0[0,:,:] = 0.
    dirichlet_bc_u0[-1,:,:] = 1.

    print("test:", dirichlet_bc_nodes)
    ohm = OhmSolver(dirichlet_bc_nodes, dirichlet_bc_u0)
    u = ohm.u(state)
    write_vti(u, "data/u.vti")

def test_curved():
    n  = (120, 60, 5)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx, origin = (-n[0]*dx[0]/2.,0.,0.))
    state = State(mesh)
    x, y, z = state.SpatialCoordinate()
    r_i = 30e-9
    r_a = 50e-9
    state.material = {"sigma": 1.}
    state.material["sigma"] = (x**2+y**2 > r_i**2) & (x**2+y**2 < r_a**2)
    write_vti(state.material, "data/material.vti")
    
    dirichlet_bc_u0 = state._zeros([n+2 for n in state.mesh.n])
    dirichlet_bc_nodes = state._zeros([n+2 for n in state.mesh.n], dtype=bool)
    dirichlet_bc_nodes[:,0,:] = True

    dirichlet_bc_u0[n[0]//2:,0,:] = 0.
    dirichlet_bc_u0[:n[0]//2,0,:] = 1.

    ohm = OhmSolver(dirichlet_bc_nodes, dirichlet_bc_u0)
    u = ohm.u(state)
    write_vti(u, "data/u.vti")
