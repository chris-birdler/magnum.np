from pathlib import Path

from magnumnp import *
import pyvista as pv
import torch
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt


base_dir = Path(__file__).resolve().parent
data_dir = base_dir / "data"
data_dir.mkdir(parents=True, exist_ok=True)

# initialize state
dt = 5e-12
n = (24, 24, 2)
l = (120e-9, 120e-9, 10e-9)
dx = (l[0]/n[0], l[1]/n[1], l[2]/n[2])
origin = (-n[0]*dx[0]/2., -n[1]*dx[1]/2., -n[2]*dx[2]/2.)
mesh = Mesh(n, dx, origin=origin)
state = State(mesh)

state.material = {
    "alpha": 0.008,
    "Ms": 800e3, # A/m
    "A": 13e-12, # J/m
    "Di": -3e-3, # J/m^2
    "Ku": 0.5e6, # J/m^3
    "Ku_axis": [0,0,1]
    }

domain1 = torch.zeros(n, dtype=torch.bool)
domain1[:,:, 0] = True
domain2 = torch.zeros(n, dtype=torch.bool)
domain2[:,:, 1] = True

state.m = state.Constant([0.,0.,0.])
state.m[domain1] = torch.tensor([0., 0., -1.])
state.m[domain2] = torch.tensor([0., 0., +1.])

# Néel skyrmion field
radius = 10e-9
x, y, z = mesh.SpatialCoordinate()
r   = torch.hypot(x, y)
phi = torch.atan2(y, x)
core = r <= radius
s = torch.pi * (r / radius)

mx = torch.cos(phi) * torch.sin(s)
my = torch.sin(phi) * torch.sin(s)
mz = torch.cos(s)
sky = torch.stack([mx, my, mz], dim=-1)

# apply to layers
mask_bot = core & domain1
mask_top = core & domain2

state.m[mask_bot] =  sky[mask_bot]   # bottom core up
state.m[mask_top] = -sky[mask_top]   # top core down

# define field terms
exchange_b = ExchangeField(domain1)
exchange_t = ExchangeField(domain2)
dmi        = InterfaceDMIField()
aniso      = UniaxialAnisotropyField()
rkky       = RKKYField(-3e-4, "z", 0, 1) # adding antiferromagnetic coupling between the two layers
bias       = ExternalField(state.Constant([65076.68505349, 46529.82981324, 0.0])) # static field bias

# calculate groundstate
try:
    mesh0, fields0 = read_vti(str(data_dir / "m0.vti"))
    state.m[...] = fields0["m0"]
except FileNotFoundError:
    with Timer("Calculate Groundstate"):
        minimizer = MinimizerBB([exchange_b, exchange_t, dmi, aniso, rkky, bias])
        minimizer.minimize(state, maxiter=5000, dm_tol=1e-4)
        state.write_vtk({"m0":state.m}, str(data_dir / "m0.vti"))
m0 = state.m.clone()

# new equilibrium under the modified bias
bias_new = ExternalField(state.Constant([65538.55364152, 45876.98754907, 0.0]))
try:
    mesh1, fields1 = read_vti(str(data_dir / "m1.vti"))
    m1 = fields1["m1"]
except FileNotFoundError:
    with Timer("Calculate Deviated Groundstate"):
        state.m = m0.clone()
        minimizer = MinimizerBB([exchange_b, exchange_t, dmi, aniso, rkky, bias_new])
        minimizer.minimize(state, maxiter=5000, dm_tol=1e-4)
        state.write_vtk({"m1": state.m}, str(data_dir / "m1.vti"))
        m1 = state.m.clone()
state.m = m0.clone()
# initial deviation corresponds to the difference between the old and new equilibrium states
delta_m = m0 - m1

# Ring-Down method
tt = torch.arange(0, 50e-9, dt)
Nt = len(tt)

try:
    stored = torch.load(str(data_dir / "ringdown_50ns.pt"), map_location=state.device)
    data4d = stored['data4d']
except FileNotFoundError:
    with Timer("Ring-Down Method "):
        llg = LLGSolver([exchange_b, exchange_t, dmi, aniso, rkky, bias_new], atol=1e-10, rtol=1e-10)
        logger = Logger(str(data_dir), ['t', 'm'], [])

        data4d = torch.zeros((Nt,)+state.m.shape)
        for i, t in enumerate(tt):
            data4d[i,...] = state.m
            llg.step(state, dt)
            logger << state
        torch.save({"data4d":data4d}, str(data_dir / "ringdown_50ns.pt"))

freq = np.fft.rfftfreq(Nt, d=dt)
freq_axis = freq[1:]
m_fft = np.fft.rfft(data4d - m0[None,...], axis=0)

# convert PSD from average-per-cell to volume integral to match the eigenmode spectrum
num_cells = np.prod(n)
cell_volume = np.prod(dx)
power = (np.abs(m_fft)**2).mean(axis=(1,2,3)).sum(axis=-1) * (num_cells * cell_volume)
peaks = scipy.signal.find_peaks(power, prominence=1e-30)[0]

# EigenSolver method
with Timer("EigenSolver"):
    try:
        res = EigenResult.load(state, str(data_dir / "eigen.pt"))
    except Exception:
        state.m = m0
        eigen = EigenSolver(state, [exchange_b, exchange_t, rkky, aniso, dmi], [bias])
        res = eigen.solve(k=20)
        res.store(str(data_dir / "eigen.pt"))
        res.save_evecs3D(str(data_dir / "evecs.pvd"))

h_excite = bias_new.h(state) - bias.h(state)
spectrum = res.spectrum(2*np.pi*freq_axis, h_excite)

# Modal projections handled by EigenResult helpers
modal_freq, modal_power = res.modal_projection_psd(delta_m, tt, volume_scale=num_cells * cell_volume)
simple_modal_power = res.simple_modal_projection(delta_m, 2*np.pi*freq_axis, volume_scale=cell_volume)

fig, ax = plt.subplots(figsize=(15,10))
##np.savez("data/ringdown_10ns.npz", f=freq_axis, p = power[1:])
#ref = np.load("data/ringdown_10ns.npz")
#ax.plot(ref["f"] * 1e-9, ref["p"], "k--", label="PSD(RingDown) 10ns", linewidth=2.0)
ax.plot(freq_axis * 1e-9, power[1:], label="PSD(RingDown)", linewidth=2.0)
ax.plot(modal_freq[1:] * 1e-9, modal_power[1:], color="green", linewidth=2.0, label="PSD(Modal projection)")
ax.plot(freq_axis * 1e-9, spectrum, color="red", linewidth=2.0, label="PSD(Harmonic drive)")
ax.plot(freq_axis * 1e-9, simple_modal_power, "--", color="purple", linewidth=2.0, label="PSD(Simple modal projection)")

ax.scatter(freq[peaks] * 1e-9, power[peaks], color="red", label="Peaks")
ax.set_xlim([0, 50])
ax.set_yscale("log")
ax.set_xlabel("Frequency [GHz]")
ax.set_ylabel("PSD [arb.]")
ax.set_title("Spatially Resolved PSD")

freq_eig = res.freq * 1e-9
tick_labels = [f"{f:.5f}" for f in freq_eig]
for p in peaks[:12]:
    x_val = freq[p] * 1e-9      # GHz
    y_val = power[p]
    ax.text(x_val, y_val, f"{freq[p]*1e-9:.2f}", rotation=45, ha='left',va='bottom')
ax.set_xticks(np.arange(0,50,5))
ax.tick_params(axis='both', direction='in', length=6, width=1.2)
ax.grid()
ax.legend(loc='upper right')

fig.savefig(base_dir / "result.png")
