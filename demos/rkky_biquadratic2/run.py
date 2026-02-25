# %% [markdown]

##  Biquadratic Coupling Demo 2
#  Another demonstration with an out-of-plane hysteresis. But we are using 2 CoFe layers this time. It take quite a lot of time to run this demo because we are relaxing the system at each magnetic field step. Please be patient.
#  In the previous demo, we used a time-dependent external field to compute the hysteresis. However, in this demo, the saturation magnetization was lowered to a realistic value. The time-depedent field method would create a coercivity point in the loop in this case. Therefore, we are using this apporach of relaxing at each field step.
#  Implemented by Jed Cheng (jed.chengXmag.ed.kyushu-u.ac.jp, replace X with @)

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
    state.m[domain1] = torch.tensor([0.99, 0.01, 0])
    state.m[domain2] = torch.tensor([-0.99, -0.01, 0])

    Ku = state.Constant([0.0]) 
    Ku[domain1] = 3.7e3
    Ku[domain2] = 3.7e3

    A = state.Constant([0.0]) 
    A[domain1] = 15e-12
    A[domain2] = 15e-12

    Ms = state.Constant([0.0])
    Ms[domain1] = 1300e3
    Ms[domain2] = 1300e3

    Ku_axis = state.Constant([1, 0, 0])
    Ku_axis[domain1] = torch.tensor([1., 0., 0.])
    Ku_axis[domain2] = torch.tensor([1., 0., 0.])

    state.material['Ms'] = Ms
    state.material['A'] = A
    state.material['Ku'] = Ku
    state.material['Ku_axis'] = Ku_axis

    exchange = ExchangeField()
    rkky     = BiquadraticRKKYField(float(J_rkky), "z", 1, 3)
    aniso    = UniaxialAnisotropyField()
    demag    = DemagField()
    external = ExternalField([0, 0, 0])

    llg = LLGSolver([demag, exchange, aniso, external, rkky])
    m1 = ("m1", lambda state: state.m[domain1])
    m2 = ("m2", lambda state: state.m[domain2])
    logger = ScalarLogger(f"data/m_sigma2_{sigma2}.dat", ['t', external.h, 'm', m1, m2])

    # calculate hysteresis
    B0 = 5
    BB = np.concatenate([np.linspace(-B0, B0, 20), np.linspace(B0, -B0, 20)])
    for j, B in enumerate(BB):
        external.h = state.Constant([0, 0, B/constants.mu_0])
        llg.relax(state, maxiter=100)           
        logger << state

Timer.print_report()
