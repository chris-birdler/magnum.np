import pytest
import torch
from magnumnp import *

def test_voronoi():
    n  = (100, 100, 1)
    dx = (1e-9, 1e-9, 1e-9)
    mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2.,-n[1]*dx[1]/2.,-n[2]*dx[2]/2.))
    state = State(mesh)
    voi = Voronoi(mesh, 10)
    domains = voi.relax()
    voi.add_intergrain_phase(2)
    #state.write_vtk(domains, "data/domains.vti")
