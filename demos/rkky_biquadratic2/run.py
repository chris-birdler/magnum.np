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
import os 

os.makedirs("data", exist_ok=True)

set_log_level(25) # show info_green, but hide info_blue

precision = torch.float32

torch.set_default_dtype(precision)
torch.set_default_device("cuda")
Timer.enable()
try:
    this_dir = pathlib.Path(__file__).resolve().parent
except:
    this_dir = pathlib.Path().resolve()


sigma2 = np.array([0.0, -0.1, -0.05, -0.01])
sigma2 = sigma2 * 2


B_min = -5
B_max = 5
n_pts = 400
B_z = np.linspace(B_min, B_max, n_pts )
B_z2 = np.linspace(B_max, B_min, n_pts )


for i, sigma2_ in enumerate(sigma2):
        Timer.enable()

        Mz = np.zeros(2*n_pts)
        Bz = np.zeros(2*n_pts)

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
        state.m[domain1] = torch.tensor([0.99, 0.01, 0], dtype=precision)
        state.m[domain2] = torch.tensor([-0.99, -0.01, 0], dtype=precision)


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
        Ku_axis[domain1] = torch.tensor([1, 0, 0], dtype=precision)
        Ku_axis[domain2] = torch.tensor([1, 0, 0], dtype=precision)


        state.material['Ms'] = Ms
        state.material['A'] = A
        state.material['Ku'] = Ku
        state.material['Ku_axis'] = Ku_axis


        exchange = ExchangeField()
        rkky     = BiquadraticRKKYField(float(J_rkky), "z", 1, 3)
        aniso    = UniaxialAnisotropyField()
        demag    = DemagField()


        for j, B in enumerate(B_z):
                external = ExternalField([0, 0, B/constants.mu_0])
                llg = LLGSolver([demag, exchange, aniso, external, rkky])
                llg.relax(state, maxiter=100)           

                domain1_m = torch.mean(state.m[domain1], dim=0)
                domain2_m = torch.mean(state.m[domain2], dim=0)

                domain1_mz = domain1_m[2].detach().cpu().numpy()
                domain2_mz = domain2_m[2].detach().cpu().numpy()

                Mz_ = (domain1_mz + domain2_mz)/2
                Mz[j] = Mz_

                if j % 50 == 0:
                        print(f"sigma2: {sigma2_}, B: {B}, Mz: {Mz_} ({j}/{len(B_z)*2})")

        
        for j, B in enumerate(B_z2):
                external = ExternalField([0, 0, B/constants.mu_0])
                llg = LLGSolver([demag, exchange, aniso, external, rkky])
                llg.relax(state, maxiter=100)
                domain1_m = torch.mean(state.m[domain1], dim=0)
                domain2_m = torch.mean(state.m[domain2], dim=0)

                domain1_mz = domain1_m[2].detach().cpu().numpy()
                domain2_mz = domain2_m[2].detach().cpu().numpy()

                Mz_ = (domain1_mz + domain2_mz)/2
                Mz[j+n_pts] = Mz_

                if j % 50 == 0:
                        print(f"sigma2: {sigma2_}, B: {B}, Mz: {Mz_} ({j+n_pts}/{len(B_z)*2})")

        np.save(f"data/hysteresis_Mz_sigma2_{sigma2_}.npy", Mz)
        np.save(f"data/hysteresis_Bz_sigma2_{sigma2_}.npy", np.concatenate((B_z, B_z2)))