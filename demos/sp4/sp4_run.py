from magnumnp import *
import torch
import pathlib

def run_sp4():
    Timer.enable()
    this_dir = pathlib.Path(__file__).resolve().parent
    
    # initialize mesh
    eps = 1e-15
    n  = (100, 25, 1)
    dx = (5e-9, 5e-9, 3e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)
    
    state.material = {
        "Ms": 8e5,
        "A": state.Constant([1.3e-11]),
        "alpha": 0.02
        }
    
    # initialize field terms
    demag    = DemagField()
    exchange = ExchangeField()
    external = ExternalField([-24.6e-3/constants.mu_0,
                              +4.3e-3/constants.mu_0,
                              0.0])
    
    # initialize magnetization that relaxes into s-state
    state.m = state.Constant([0,0,0])
    state.m[1:-1,:,:,0]   = 1.0
    state.m[(-1,0),:,:,1] = 1.0
    
    # relax without external field
    llg = LLGSolver([demag, exchange])
    llg.relax(state)
    write_vti(state.m, "data/m0.vti", state)
    
    # perform integration with external field
    llg = LLGSolver([demag, exchange, external])
    logger = Logger(this_dir / "data", ['t', 'm'])
    while state.t < 1e-9-eps:
        llg.step(state, 1e-11)
        logger << state
    
    Timer.print_report()
    
if __name__ == "__main__":
    run_sp4()    
