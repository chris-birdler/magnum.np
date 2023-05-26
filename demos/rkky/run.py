# Example and benchmark problem according to
#
# "Accurate finite-difference micromagnetics of magnets including RKKY interaction",
# Suess, Dieter, et al.,  arXiv preprint arXiv:2206.11063 (2022).
#
# The analytical saturation field Hx = 5*Hk

from magnumnp import *
import torch
import numpy as np

Timer.enable()
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

domain1 = state._zeros(n, dtype=torch.bool)
domain1[:,:,n[2]//2:] = True

domain2 = state._zeros(n, dtype=torch.bool)
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
zeeman    = ExternalField(TimeInterpolator(state, {0e-0:   [0, 0, Hxmin],
                                                   tfinal: [0, 0, Hxmax]}))

# integrate
llg = LLGSolver([aniso, exchange1, exchange2, rkky, zeeman], atol = 1e-6)
logger = Logger("data", ['t', 'm', zeeman.h], ['m'], fields_every = 100)
while state.t < tfinal:
    logger << state
    llg.step(state, 1e-9)

Timer.print_report()

#plot the data
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

data = np.loadtxt("data/log.dat")
ref = np.loadtxt("ref/log.dat")

Hk = 2 * 1e5 / 1.

fig, ax = plt.subplots(figsize=(10,5))
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

ax.plot(data[:,6]/Hk, data[:,3], '-', color = cycle[2], label = "magnum.np")
ax.plot(ref[:,6]/Hk, ref[:,3], '-', color = cycle[2], linewidth = 6, alpha = 0.4, label = "reference")

ax.set_xlim([4.82,5.0])
ax.set_ylim([0.99994,1.0])
ax.set_yticklabels(['0.99994', '0.99995', '0.99996', '0.99997', '0.99998', '0.99999', '1.00000'])
ax.set_title("RKKY Standard Problem")
ax.set_xlabel("External Field H$_{ext}$ / H$_k$ [1]")
ax.set_ylabel("Magnetization m")
ax.legend(ncol=3)
ax.grid()
fig.savefig("data/results.png")
