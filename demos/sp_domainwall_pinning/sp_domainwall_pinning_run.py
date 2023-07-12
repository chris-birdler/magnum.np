from magnumnp import *
import numpy as np
import pathlib

def run_sp_domainwall_pinning():
    Timer.enable()
    this_dir = pathlib.Path(__file__).resolve().parent

    Hextmax=1.8/constants.mu_0
    Hextmin=1.4/constants.mu_0
    tfinal = 20e-9

    n  = (80, 1, 1)
    dx = (1e-9, 1e-9, 1e-9)
    origin = (-n[0]*dx[0]/2., -n[1]*dx[1]/2., -n[2]*dx[2]/2.,)
    mesh = Mesh(n, dx, origin)
    state = State(mesh)
    state.material = {"alpha": 1.}
    x, y, z = state.SpatialCoordinate()

    soft = (x < 0)
    hard = (x >= 0)

    Ms = state.Constant([0.0])
    Ms[soft] = 0.25/constants.mu_0
    Ms[hard] = 1./constants.mu_0

    Ku = state.Constant([0.0])
    Ku[soft] = 1e5
    Ku[hard] = 1e6

    A = state.Constant([0.0])
    A[soft] = 0.25e-11
    A[hard] = 1.e-11

    state.material['Ms'] = Ms
    state.material['A'] = A
    state.material['Ku'] = Ku
    state.material['Ku_axis'] = [0,1,0]

    state.m = state.Constant([np.sin(0.3), np.cos(0.3), 0.0])
    state.m[hard] = state.Tensor([0.0, -1.0, 0.0])
    state.m.normalize()

    exchange = ExchangeField()
    aniso    = UniaxialAnisotropyField()
    external = ExternalField([0, 0, 0])
#    external = ExternalField(lambda t: state.Constant([0, (Hextmax-Hextmin)*t/tfinal+Hextmin, 0]))

#    llg = LLGSolver([exchange, aniso, external])
    minimizer = MinimizerBB([exchange, aniso, external])
    logger = ScalarLogger(this_dir / "data" / "m.dat", ['t', external.h, 'm'])
    for h in state.linspace(Hextmin, Hextmax, steps=100):
        external.h = [0, h, 0]
        minimizer.minimize(state)
        logger << state

    Timer.print_report()

if __name__=="__main__":
    run_sp_domainwall_pinning()
