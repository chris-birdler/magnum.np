# %% [markdown]
# # strain jump over interfaces

# %% [markdown]
# ## Run Simulation

# %%
from magnumnp import *
import numpy as np
import torch
import matplotlib.pyplot as plt

set_log_level(25) # show info_green, but hide info_blue

# problem setup
u0 = 10e-11
uL = -5e-11

CI = 2.4e11
CII = 1.e11

dxL = 1e-9
dxR = 4e-9

N = 30 # number of cells
NL = 20
NR = N-NL

a = NL*dxL # interface position

cells_L = dxL*np.ones(NL)
cells_R = dxR*np.ones(NR)

dxx = np.concatenate((cells_L,cells_R))
L = np.sum(dxx)

n = (N,1,1)
dx = (dxx, dxL, dxL)

# analytical solution
A1 = CII*(uL-u0)/(CI*(L-a)+CII*a)
A2 = CI*A1/CII

D1 = u0
D2 = A1*a*(1.-CI/CII) + u0

# setup state
mesh = Mesh(n, dx)
state = State(mesh)

C_iso = C_isotropic(1,0.)
C = state.Constant(C_iso)
C[:NL] *= CI
C[NL:] *= CII

state.material = {"C" : C,
                  "lambda_100" : state.Constant(0.),
                  "lambda_111" : state.Constant(0.),
                  "A" : state.Constant(1.)}

x = mesh.SpatialCoordinate()[0][:,0,0]

state.m = state.Constant((1.,0.,0.))

state.ud = state.Constant((0.,0.,0.))
state.ud[:,0,0,0] = x*A1+D1
state.ud[NL:,0,0,0] = x[NL:]*A2+D2

# First derivative
eps0 = epsilon(state)[:,0,0,0] # with jump conditions
gx = torch.gradient(state.ud[...,0], spacing=(x,), dim=(0,))[0][:,0,0] # using torch.gradient

# displacement
ux = state.ud[:,0,0,0].cpu().numpy()


# plot data
fig, axs = plt.subplots(2,1,gridspec_kw = {'wspace':0, 'hspace':0})

for ax in axs:
    ax.set_xlim([0, L*1e9])
    ax.axvspan(0, a*1e9, color='blue', alpha=0.3)
    ax.axvspan(a*1e9, L*1e9, color='grey', alpha=0.3)

axs[0].scatter(x*1e9, ux*1e12, s=15, c="C1", label="$u_x$")
axs[1].scatter(x*1e9, 1e4*eps0, s=15, c="C1", label="$\epsilon_{xx}$")
axs[1].scatter(x*1e9, 1e4*gx, marker="x", c="k", s=11, label="torch.grad")
axs[1].plot(x*1e9, 1e4*A1*np.ones(N), ls=":", c="k", lw=0.7)
axs[1].plot(x*1e9, 1e4*A2*np.ones(N), ls="--", c="k", lw=0.7)

axs[0].set_xticklabels([])
axs[1].set_xlabel("x (nm)")

for ax in axs:
    ax.legend()

axs[0].set_ylabel("u (pm)")
axs[1].set_ylabel("Strain x 1e4")

fig.tight_layout()
plt.savefig("strain_jump.pdf")
