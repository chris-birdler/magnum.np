from magnumnp import *
import torch
import pathlib
from math import sin, cos, pi

def run_sp_FMR():
    Timer.enable()
    this_dir = pathlib.Path(__file__).resolve().parent

    # initialize mesh
    eps = 1e-15
    n = (24, 24, 2)
    l = (120e-9, 120e-9, 10e-9)
    dx = (l[0]/n[0], l[1]/n[1], l[2]/n[2])
    origin = (-n[0]*dx[0]/2., -n[1]*dx[1]/2., -n[2]*dx[2]/2.)
    mesh = Mesh(n, dx, origin=origin)

    # initialize state
    state = State(mesh)
    state.material = {
        "alpha": 1.,
        "Ms": 800e3,
        "A": 13e-12
        }
    state.m = state.Constant([0, 0, 1])

    #relax state
    demag    = DemagField()
    exchange = ExchangeField()
    bias     = ExternalField(80e3*state.Tensor([1., 0.715, 0]).normalize())

    llg = LLGSolver([demag, exchange, bias])
    llg.relax(state)

    #integrate
    state.t = 0.0
    state.material["alpha"] = 0.008
    bias = ExternalField(80e3*state.Tensor([1., 0.7, 0]).normalize())

    llg = LLGSolver([demag, exchange, bias])
    logger = Logger(this_dir / "data", ['t', 'm'])
    while state.t < 10e-9-eps:
        llg.step(state, 5e-12)
        logger << state

    Timer.print_report()

if __name__=="__main__":
    run_sp_FMR()
