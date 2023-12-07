from magnumnp import *
import torch
import pathlib

def run_sp5():
    Timer.enable()
    this_dir = pathlib.Path(__file__).resolve().parent

    # initialize state
    n  = (40, 40, 1)
    dx = (2.5e-9, 2.5e-9, 10e-9)
    mesh = Mesh(n, dx)
    
    state = State(mesh)
    state.material = {
        "Ms": state.Constant([8e5]),
        "A": state.Constant([1.3e-11]),
        "alpha": state.Constant([0.1]),
        "xi": state.Constant([0.05]),
        "b": state.Constant([72.17e-12])
        }
    
    # initialize magnetization that relaxes into s-state
    state.m = state.Constant([0,0,0])
    state.m[:20,:,:,1] = -1.
    state.m[20:,:,:,1] = 1.
    state.m[20,20,:,1] = 0.
    state.m[20,20,:,2] = 1.
    
    state.j = state.Tensor([1e12, 0, 0])
    
    # initialize field terms
    demag    = DemagField()
    exchange = ExchangeField()
    torque   = SpinTorqueZhangLi()
    
    # initialize sstate
    llg = LLGSolver([demag, exchange])
    llg.relax(state)
    write_vti(state.m, "data/m0.vti", state)
    
    # perform integration with spin torque
    llg = LLGSolver([demag, exchange, torque])
    logger = Logger(this_dir / "data", ['t', 'm'])
    while state.t < 5e-9:
        llg.step(state, 1e-11)
        logger << state
    
    Timer.print_report()
    
if __name__=="__main__":
    run_sp5()    
