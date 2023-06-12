from magnumnp import *
import torch

def run_sot():
    Timer.enable()
    
    # initialize mesh
    eps = 1e-15
    #n  = (2, 2, 2)
    n  = (1, 1, 1)
    dx = (1.0e-9, 1.0e-9, 1e-9)
    
    mesh = Mesh(n, dx)
    state = State(mesh)
    
    # initialize polarization, p, and charge current amplitude
    # thickness of thin film on which the SOT acts
    p = state.Tensor((0, -1, 0))
    je = 6.9e10
    d = n[2] * dx[2]
    Keff = 1200e3*constants.mu_0*0.4/2./constants.mu_0
    
    state.material = {
        "Ms": 1200e3,
        "A": 15e-12,
        "Ku": Keff,
        "Ku_axis": [0, 0, 1],
        "gamma": 2.211e5,
        "alpha": 0.048,
        "eta_damp": -0.1, # both eta with opposite sign as magnum.af, same as magnum.pi
        "eta_field": 0.3,
        "p": p,
        "d": d,
        "je": je}
    
    # initialize field terms
    exchange = ExchangeField()
    aniso = UniaxialAnisotropyField()
    torque = SpinOrbitTorque()
    
    # initialize magnetization that relaxes into s-state
    state.m = state.Constant([0,0,1])
    
    # relax without external field
    llg = LLGSolver([exchange, aniso])
    llg.relax(state)
    
    # perform integration with external field
    state.t = 0.
    llg = LLGSolver([exchange, torque, aniso])
    slogger = ScalarLogger("data/log.dat", ['t', 'm', torque.h])
    flogger = FieldLogger("data/fields.pvd", ['m'])
    
    while state.t < 1e-9-eps:
        slogger << state
        flogger << state
        llg.step(state, 1e-12)
    
    Timer.print_report()

if __name__=="__main__":
    run_sot()