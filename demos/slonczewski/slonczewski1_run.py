### Slonczewski STT Extension
#
# Example taken from the Mumax3 paper: Vansteenkiste, Arne, et al. "The design and verification of MuMax3." AIP advances 4.10 (2014).
# The team compared the results of Mumax3 with the results from oommf of a 160nm x 80nm x 5nm permalloy film. Non-zero epislon prime is used.
# Implemented by Jed Cheng (jed.cheng@mag.ed.kyushu-u.ac.jp)

from magnumnp import *
import torch
import numpy as np

Timer.enable()

# initialize state
n = (64, 32, 1)
dx = (2.5e-9, 2.5e-9, 5e-9)
L = (n[0]*dx[0], n[1]*dx[1], n[2]*dx[2])
mesh = Mesh(n, dx)

state = State(mesh)
state.material = {
    "Ms": 8e5,
    "A": 1.3e-11,
    "alpha": 0.01,
    "P": 0.5669,
    "Lambda": 2,
    "epsilon_prime": 1,
    "mp": state.Tensor((np.cos(np.radians(20)), np.sin(np.radians(20)), 0)),
    "d": L[2],
    "J": -4e11,
    }
state.m = state.Constant([1,0,0])

# initialize field terms
demag    = DemagField()
exchange = ExchangeField()
torque   = SpinTorqueSlonczewski()

# perform integration with spin torque
llg = LLGSolver([demag, exchange, torque])
logger = Logger("data2", ['t', 'm'], ["m"], fields_every = 10)
while state.t < 0.5e-9:
    llg.step(state, 5e-12)
    logger << state

Timer.print_report()
