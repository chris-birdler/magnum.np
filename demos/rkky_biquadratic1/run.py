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

precision = torch.float32

torch.set_default_dtype(precision)
torch.set_default_device("cuda")
Timer.enable()
try:
    this_dir = pathlib.Path(__file__).resolve().parent
except:
    this_dir = pathlib.Path().resolve()


sigma2 = np.array([0, -0.1, -0.05, -0.01])
sigma2 = sigma2 * 2

for sigma2_ in sigma2:

        Timer.enable()

        J_rkky = sigma2_

        # initialize mesh
        n  = (20, 20, 5)
        dx = (3e-9, 3e-9, 2e-9)
        mesh = Mesh(n, dx)
        state = State(mesh)

        state.material = {
                "alpha": 1
                }

        domain1 = state.Constant(False, dtype=torch.bool)
        domain2 = state.Constant(False, dtype=torch.bool)
        domain1[:,:,0:2] = True
        domain2[:,:,3:5] = True

        # Set initial magnetization
        state.m = state.Constant([0, 0, 0])
        size_ref = state.m[domain1].shape
        state.m[domain1] = torch.tensor([-0.99, 0.01, 0], dtype=precision)
        state.m[domain2] = torch.tensor([0.99, 0.01, 0], dtype=precision)


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
        Ku_axis[domain1] = torch.tensor([1, 0, 0], dtype=precision)
        Ku_axis[domain2] = torch.tensor([1, 0, 0], dtype=precision)


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


        external = ExternalField(TimeInterpolator(state, {0.0e-9: [0.0, 0.0, 0.0],
                                                        1.0e-9: [0.0, 0.0, 0.0],
                                                        3.5e-9: [0, 0.0, 20/constants.mu_0],
                                                        8.5e-9: [0, 0.0, -20/constants.mu_0],
                                                        13.5e-9: [0, 0.0, 20/constants.mu_0]}))

        llg = LLGSolver([demag, aniso, exchange1, exchange2, external, rkky])
        logger = ScalarLogger("data/m.dat", ['t', external.h, 'm'])
        eps = 1e-15

        time_np = np.zeros((int(13.5e-9/1e-10)))
        domain1_np = np.zeros((int(13.5e-9/1e-10), 3))
        domain2_np = np.zeros((int(13.5e-9/1e-10), 3))
        i = 0

        while state.t < 13.5e-9-eps:
            llg.step(state, 1e-10)
            logger << state

            domain1_m = torch.mean(state.m[domain1], dim=0)
            domain2_m = torch.mean(state.m[domain2], dim=0)

            domain1_np[i] = domain1_m.detach().cpu().numpy()
            domain2_np[i] = domain2_m.detach().cpu().numpy()
            time_np[i] = state.t
            i += 1

        data = np.loadtxt("data/m.dat")
        B = data[:,3]*4*np.pi*1e-7

        mz1 = domain1_np[:,2]
        mz2 = domain2_np[:,2]

        M_total = (mz1*2/3 + mz2*1/3)

        np.save(f"data/hysteresis_Mz_sigma2_{sigma2_}.npy", M_total)
        np.save(f"data/hysteresis_Bz_sigma2_{sigma2_}.npy", B)