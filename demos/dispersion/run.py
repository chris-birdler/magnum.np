# %% [markdown]
# # Calculation of Magnon Dispersion

# %% [markdown]
# ## Run Simulation

# %%
from magnumnp import *
import torch
import pathlib
from tqdm import tqdm

set_log_level(25) # show info_green, but hide info_blue
torch.set_default_dtype(torch.float32)
Timer.enable()
try:
    this_dir = pathlib.Path(__file__).resolve().parent
except:
    this_dir = pathlib.Path().resolve()

# initialize state
Nt = 500
dt = 1e-12

n = 500, 25, 1
dx = 2e-9, 2e-9, 1e-9

Ms = 8e5
A = 1.3e-11
alpha = 1e-8

h_bias   = [804e3, 0.0, 0.0]
h_excite = [0.0, 0.0, 400e3]

origin = (-n[0]*dx[0]/2, -n[1]*dx[1]/2, -n[2]*dx[2]/2)
mesh = Mesh(n, dx, origin)

state = State(mesh)
x,y,z = mesh.SpatialCoordinate()
mt = torch.zeros((Nt, n[0], n[1], n[2]))

state.material = {
    'Ms': Ms,
    'A': A,
    'alpha': alpha
    }

state.m = state.Constant(h_bias)
normalize(state.m)

# initialize spin wave excitation field
N = n[0] / 2
kc = N/2. * (2.*torch.pi/(n[0]*dx[0]))
fc = 1/(2*dt)

h0 = state.Constant(h_excite)
h0 *= Expression(torch.sinc(kc / torch.pi * x))
h0 *= Expression(torch.sinc(kc / torch.pi * y))

omega = 2.*fc
t0 = 50e-12
ht = lambda state: torch.sinc(omega * (state.t-t0)) * h0

# initialize field terms
demag    = DemagField()
exchange = ExchangeField()
bias     = ExternalField(h_bias)
excite   = ExternalField(ht)

# initialize LLG solver and relax state
llg = LLGSolver([demag, exchange, bias])
llg.relax(state)

# apply spin wave excitation
state.t = 0.0
llg = LLGSolver([demag, exchange, bias, excite])
logger = Logger(this_dir / "data", ['t', 'm'])

for i in tqdm(range(Nt)):
    mt[i,...] = state.m[...,2]
    llg.step(state, dt)
    logger << state

torch.save(mt, "data/mt_saved.pt")
Timer.print_report()
