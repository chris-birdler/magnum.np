from magnumnp import *
import torch
import numpy as np
import pathlib

def run_DMI():
    Timer.enable()
    this_dir = pathlib.Path(__file__).resolve().parent

    n  = (100, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)
    origin = (-n[0]*dx[0]/2., -n[1]*dx[1]/2., -n[2]*dx[2]/2.,)

    mesh = Mesh(n, dx, origin)
    state = State(mesh)
    state.material = {
        "alpha": 1.,
        "A": 13e-12,
        "Ms": 800e3,
        "Di": -3e-3,
        "Ku": 0.4e6,
        }
    state.material['Ku_axis'] = [0,0,1]
    x, y, z = state.SpatialCoordinate()

    write_vti(state.material, "data/material.vti", state)

    state.m = state.Constant([-0.1, 0.0, 0.9])
    state.m.normalize()

    exchange = ExchangeField()
    aniso    = UniaxialAnisotropyField()
    dmi      = InterfaceDMIField()

    llg = LLGSolver([exchange, aniso, dmi])
    logger = ScalarLogger("data/m_relax.dat", ['t', 'm'])
    while state.t <= 5e-9:
        logger << state
        llg.step(state, 1e-12)
    np.savetxt(this_dir / "data" / "m0_magnumnp.dat", torch.concat((x[:,0,0,None], state.m[:,0,0,:]), axis=1).cpu().numpy())

    Timer.print_report()

if __name__=="__main__":
    run_DMI()
