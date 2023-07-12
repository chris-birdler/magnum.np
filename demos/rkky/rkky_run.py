# Example and benchmark problem according to
#
# "Accurate finite-difference micromagnetics of magnets including RKKY interaction",
# Suess, Dieter, et al.,  arXiv preprint arXiv:2206.11063 (2022).
#
# The analytical saturation field Hx = 5*Hk

from magnumnp import *
import torch
import pathlib

def run_rkky():
    Timer.enable()
    this_dir = pathlib.Path(__file__).resolve().parent

    A = 1e-11
    Js = 1.0
    Ms = Js/constants.mu_0
    Ku = 1e5
    tfinal = 20e-9
    J_rkky = -2.0e-3
    Hxmin = 4.8 * 2*Ku/Js
    Hxmax = 5.0 * 2*Ku/Js

    # initialize mesh
    n  = (1, 1, 400)
    dx = (2e-9, 2e-9, 2e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)

    state.material = {
            "Ms": Ms,
            "A": A,
            "Ku": Ku,
            "Ku_axis": state.Tensor([0,1,0]),
            "alpha": 1.0
            }

    domain1 = state.zeros(n, dtype=torch.bool)
    domain1[:,:,n[2]//2:] = True

    domain2 = state.zeros(n, dtype=torch.bool)
    domain2[:,:,:-n[2]//2] = True

    # Set initial magnetization
    state.m = state.Constant([0, 0, 0])
    state.m[domain1] = state.Tensor([0, -1, 0])
    state.m[domain2] = state.Tensor([0, 1, 0])

    # define field terms
    exchange1 = ExchangeField(domain1) # Two seperate exchange regions are required,
    exchange2 = ExchangeField(domain2) # so that at the RKKY interface the bulk exchange is zero
    rkky      = RKKYField(J_rkky, "z", n[2]//2-1, n[2]//2, order=2)
    aniso     = UniaxialAnisotropyField()
    zeeman    = ExternalField([0, 0, Hxmin])

    minimizer = MinimizerBB([aniso, exchange1, exchange2, rkky, zeeman])
    logger = Logger(this_dir / "data", ['t', 'm', zeeman.h], ['m'], fields_every = 100)
    for h in state.linspace(Hxmin, Hxmax, steps=100):
        zeeman.h = [0, 0, h]
        minimizer.minimize(state)
        logger << state

    Timer.print_report()

if __name__ == "__main__":
        run_rkky()
