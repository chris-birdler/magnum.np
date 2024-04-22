from magnumnp import *
import numpy as np
import pathlib

"""
Stochastic Time-Integration Demo

based on Leliaert et al., "Adaptively time stepping the stochastic Landau-Lifshitz-Gilbert equation at nonzero temperature:
                           Implementation and validation in MuMax 3", AIP advances 7, 125010 (2017)

For comparison with the paper Leliaert et al., AIP advances 7, 125010 (2017) Figure 5 we need to simulate a 2**18 (64**3) of uncoupled cells.
Since the cell are uncoupled the volume in our calculations has to be the volume of one cell and not the whole simulated material.
"""

def run_langevin():
    Timer.enable()
    this_dir = pathlib.Path(__file__).resolve().parent

    n  = (64, 64, 64) # 2**18 of uncoupled cells.
    dx = (10e-9, 10e-9, 10e-9)
    mesh = Mesh(n, dx)

    Ms = 1e6
    t_final = 2e-8

    state = State(mesh)
    state.material = {
        "Ms": state.Constant(Ms), # Js ~ 1.26 T for Ms = 1e6 A/m
        "alpha": state.Constant(0.1)
        }
    state.m = state.Constant([0,0,0])
    add_noise(state.m)
    normalize(state.m)

    external = ExternalField(state.Constant([0,0,0]))
    thermal = ThermalField()
    llg = LLGSolver([external, thermal], solver = RKF45)

    # perform stochastic integration
    for xi in [30, 91, 242, 725]:
        logger = ScalarLogger("data/log_xi%d.dat" % xi, ['t', external.h, 'm'])
        state.T = constants.mu_0 * Ms * state.mesh.cell_volumes /constants.mu_0 / (constants.kb * xi)
        print(f"Running for xi = {xi} (T = {state.T})")
        for h in np.linspace(0, 0.1, num=11):
            print(f"h = {h}")
            external.h = state.Constant([h / constants.mu_0, 0, 0])
            llg.step(state, dt = t_final)
            logger << state
    Timer.print_report()

if __name__ == "__main__":
    run_langevin()
