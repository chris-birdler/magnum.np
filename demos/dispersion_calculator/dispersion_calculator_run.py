import torch
from magnumnp import *
from tqdm import tqdm
import pathlib

def run_dispersion_calculator():
    Timer.enable()
    this_dir = pathlib.Path(__file__).resolve().parent

    Nt = 500
    dt = 1e-12
    
    n = 500, 25, 1
    dx = 2e-9, 2e-9, 1e-9
    
    Ms = 8e5
    A = 1.3e-11
    alpha = 1e-8
    
    h_bias   = [804e3, 0.0, 0.0]
    h_excite = [0.0, 0.0, 400e3]
    
    
    # define mesh
    origin = (-n[0]*dx[0]/2, -n[1]*dx[1]/2, -n[2]*dx[2]/2)
    mesh = Mesh(n, dx, origin)
    
    # initialize state
    state = State(mesh, dtype = torch.float32)
    x,y,z = state.SpatialCoordinate()
    mt = state.zeros((Nt, n[0], n[1], n[2]))
    
    state.material = {
        'Ms': Ms,
        'A': A,
        'alpha': alpha
        }
    
    state.m = state.Constant(h_bias)
    state.m.normalize()
    
    # initialize spin wave excitation field
    N = n[0] / 2
    kc = N/2. * (2.*torch.pi/(n[0]*dx[0]))
    fc = 1/(2*dt)
    
    h0 = state.Constant(h_excite)
    h0 *= torch.special.sinc(kc / torch.pi * x).unsqueeze(-1)
    h0 *= torch.special.sinc(kc / torch.pi * y).unsqueeze(-1)
    
    omega = 2.*fc
    t0 = 50e-12
    ht = lambda t: torch.sinc(omega * (t-t0)) * h0
    
    # initialize energy terms
    demag    = DemagField()
    exchange = ExchangeField()
    bias     = ExternalField(state.Constant(h_bias))
    excite   = ExternalField(ht)
    
    # initialize LLG solver and relax state
    set_log_level(100)
    llg = LLGSolver([demag, exchange, bias])
    llg.relax(state)
    
    # apply spin wave excitation
    state.t = 0.0
    llg = LLGSolver([demag, exchange, bias, excite])
    logger = Logger(this_dir / "data", ['t', 'm'])
    
    for i in tqdm(range(Nt)):
        mt[i,...] = state.m[...,2]
        llg.step(state, dt)
        logger << state
    
    set_log_level(30)
    torch.save(mt, "data/mt_saved.pt")


if __name__ == "__main__":
    run_dispersion_calculator()
