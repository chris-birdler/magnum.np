# %% [markdown]
# # spin wave excited by an acoustic bulk mode

# %% [markdown]
# ## Run Simulation

# %%
from magnumnp import *
import numpy as np 
from tqdm import tqdm

set_log_level(25) # show info_green, but hide info_blue

# Constants
u0 = 2e-12
L_factor = 20
f = 4e9
N_wl = 200
N_T = 450
Ms = 490e3
Aex = 8e-12
alpha = 0.05
eta=0
rho = 8900.
C11 = 2.54e11
C12 = 1.55e11
C44 = 1.23e11
lambda_100 = -4.6e-5
lambda_111 = -2.4e-5

# Derived quantities
omega = 2.*np.pi * f
c = np.sqrt(C44 / rho)
k = omega / c
wl = 2.*np.pi / k

Jex = 2.*Aex/(Ms*constants.mu_0)
h0 = 2.*np.pi*f/constants.gamma - Jex*k**2.

Lx = L_factor * wl
dx0 = wl / N_wl
dx = (dx0, dx0, dx0)
n = (int(Lx/dx[0]), 1, 1)

T = 1. / f
dt = T / N_T
print("CFL condition is %f" % (c * dt / dx0))

# setup state
mesh = Mesh(n, dx, pbc=(0,1,1))
state = State(mesh)
state.material = {
        "alpha": state.Constant(alpha),
        "Ms": state.Constant(Ms),
        "A": state.Constant(Aex),
        "rho": state.Constant(rho),
        "C": state.Constant(C_cubic(C11, C12, C44)),
        "lambda_100": state.Constant(lambda_100),
        "lambda_111": state.Constant(lambda_111),
        "eta": state.Constant(eta)
    }


# setup time integration
bias = ExternalField(state.Constant((0,0,h0)))
exchange = ExchangeField()
magEl = MagnetoElasticField()
h_terms = [bias, exchange, magEl]

llg = LLGWithLESolver(h_terms, C_sym="cubic", iteration_depth=0, boundary_nodes=3, dt=dt)

# linear elasticity boundary conditions
# ... fixed right interface
u_mask = torch.zeros(n, dtype=torch.bool)
u_mask[-1] = 1
u_cond = lambda state : torch.zeros((n[1], n[2], 3))

dirichlet_bcs = []
dirichlet_bcs.append(BC(u_mask, u_cond))

# ... left shear traction
def t0_cond(state):
    nx, ny, nz = state.mesh.n

    t0 = u0 * C44 * k

    t_in = torch.zeros((ny, nz, 3))
    t_in[:,:,2] = t0 * np.sin(float(state.t) * omega)

    return t_in

neumann_bcs = []
neumann_bcs.append(PlaneBC(Plane(0, 0, -1), t0_cond))

state.bcs = dict()
state.bcs["ud"] = dirichlet_bcs
state.bcs["t"] = neumann_bcs

# initial conditions
state.m = state.Constant((0.,0.,1.))
state.pd = state.Constant((0.,0.,0.))
state.ud = state.Constant((0.,0.,0.))

# ... initial displacement due to magnetic strain
kx_x = -0.5*lambda_100
kx_y = -0.5*lambda_100
kx_z = lambda_100
kx = kx_x + (C12/C11)*(kx_y + kx_z)

x = mesh.SpatialCoordinate()[0]
state.ud[...,0] = -((Lx-0.5*dx[0])-x)*kx

# logging
def dE_in(state):
    v_in = state.pd[0] / rho
    t_in = t0_cond(state)
    At = dx[1]*dx[2]
    return At*(t_in*v_in).sum()

scalar_log_list = ["t", # time
                   llg.E, # energy in the magnetic field terms
                   llg.U, # total potential energy
                   llg.T_el, # elastic kinetic energy
                   dE_in] # energy/time added by the boundary condition

logger = Logger("log", scalars=scalar_log_list, fields=["m", "ud", "pd"], fields_every=N_T)

# run simulation
state.t = 0.
M = int(L_factor * N_T)
for i_t in tqdm(range(M)):
    logger << state
    llg.step(state, dt)

    normalize(state.m)

logger << state
