# %% [markdown]

##  Biquadratic Coupling Demo 1
#  Example taken from Ubermag
#  https://ubermag.github.io/examples/notebooks/rkky.html
# But I reduced the M_sat of a layer by a half to better demonstrate the effect of the exchange coupling.
# Implemented by Jed Cheng (jed.chengXmag.ed.kyushu-u.ac.jp, replace X with @)

# %% [markdown]
# ## Run Simulation

# %%
from magnumnp import *
import torch
import numpy as np
import pathlib
from tqdm import tqdm

set_log_level(25) # show info_green, but hide info_blue
torch.set_default_dtype(torch.float32)

Timer.enable()
try:
    this_dir = pathlib.Path(__file__).resolve().parent
except:
    this_dir = pathlib.Path().resolve()


for sigma2 in np.array([0, -0.20, -0.10, -0.02]):
    J_rkky = sigma2

    # initialize mesh
    dt = 1e-10
    n  = (20, 20, 5)
    dx = (3e-9, 3e-9, 2e-9)
    mesh = Mesh(n, dx)
    state = State(mesh)

    state.material = {"alpha": 1}

    domain1 = state.Constant(False, dtype=torch.bool)
    domain2 = state.Constant(False, dtype=torch.bool)
    domain1[:,:,0:2] = True
    domain2[:,:,3:5] = True

    # Set initial magnetization
    state.m = state.Constant([0, 0, 0])
    size_ref = state.m[domain1].shape
    state.m[domain1] = torch.tensor([-0.99, 0.01, 0])
    state.m[domain2] = torch.tensor([0.99, 0.01, 0])

    Ku = state.Constant([0.0])
    Ku[domain1] = 1e5
    Ku[domain2] = 1e5

    A = state.Constant([0.0])
    A[domain1] = 1e-12
    A[domain2] = 1e-12

    Ms = state.Constant([0.0])
    Ms[domain1] = 8e6
    Ms[domain2] = 4e6

    Ku_axis = state.Constant([1, 0, 0])
    Ku_axis[domain1] = torch.tensor([1., 0., 0.])
    Ku_axis[domain2] = torch.tensor([1., 0., 0.])

    state.material['Ms'] = Ms
    state.material['A'] = A
    state.material['Ku'] = Ku
    state.material['Ku_axis'] = Ku_axis

    exchange1 = ExchangeField(domain1)
    exchange2 = ExchangeField(domain2)
    rkky      = BiquadraticRKKYField(float(J_rkky), "z", 1, 3)
    aniso     = UniaxialAnisotropyField()
    demag     = DemagField()

    llg = LLGSolver([demag, exchange1, exchange2, aniso, rkky])
    llg.relax(state)

    external = ExternalField(TimeInterpolator(state, 
                            {0.0e-9: [0.0, 0.0, 0.0],
                             1.0e-9: [0.0, 0.0, 0.0],
                             3.5e-9: [0, 0.0, 20/constants.mu_0],
                             8.5e-9: [0, 0.0, -20/constants.mu_0],
                             13.5e-9: [0, 0.0, 20/constants.mu_0]}))

    llg = LLGSolver([demag, aniso, exchange1, exchange2, external, rkky])
    m1 = ("m1", lambda state: state.m[domain1])
    m2 = ("m2", lambda state: state.m[domain2])
    logger = ScalarLogger(f"data/m_sigma2_{sigma2}.dat", ['t', external.h, 'm', m1, m2])

    for i in tqdm(torch.arange(0, 13.5e-9, dt)):
        llg.step(state, dt)
        logger << state

Timer.print_report()
