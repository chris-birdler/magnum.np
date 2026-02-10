from magnumnp import *
import torch
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt

### initialize state
Timer.enable()
dt = 2e-12
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

# Neel skyrmion field
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
rkky       = RKKYField(-3e-4, "z", 0, 1)
bias       = ExternalField(state.Constant([65076.68, 46529.82, 0.0]))
excite     = ExternalField(state.Constant([461.87, -652.84, 0.0]))


### calculate groundstate
with Timer("Calculate Groundstate"):
    try:
        mesh0, fields0 = read_vti("data/m0.vti")
        m0 = fields0["m0"]
    except FileNotFoundError:
        minimizer = MinimizerBB([exchange_b, exchange_t, dmi, aniso, rkky, bias])
        minimizer.minimize(state, maxiter=5000, dm_tol=1e-4)
        state.write_vtk({"m0":state.m}, "data/m0.vti")
        m0 = state.m.clone()


### calculate groundstate with excitation
with Timer("Calculate Groundstate with Excitation"):
    try:
        mesh1, fields1 = read_vti("data/m1.vti")
        m1 = fields1["m1"]
    except FileNotFoundError:
        minimizer = MinimizerBB([exchange_b, exchange_t, dmi, aniso, rkky, bias, excite])
        minimizer.minimize(state, maxiter=5000, dm_tol=1e-4)
        state.write_vtk({"m1": state.m}, "data/m1.vti")
        m1 = state.m.clone()

delta_m = m0 - m1
state.m = m0


### use time-domain simulation with sinc excitation
with Timer("Time-domain simulation with sinc excitation"):
    tt = torch.arange(0, 10e-9, dt)
    freq = np.fft.rfftfreq(len(tt), d=dt)
    f_max = 100e9  # Hz, maximum frequency of interest
    t_pulse = 5e-11 # pulse center time (at t=0)

    try:
        stored = torch.load("data/sinc.pt", map_location=state.device)
        data4d_sinc = stored['data4d_sinc']
    except FileNotFoundError:
        # Store field and response history
        data4d_sinc = torch.zeros((len(tt),) + state.m.shape)
        bias_pulse = ExternalField(lambda state: excite.h(state) * np.sinc(2 * f_max * (state.t-t_pulse)))
        llg_sinc = LLGSolver([exchange_b, exchange_t, dmi, aniso, rkky, bias, bias_pulse], atol=1e-10, rtol=1e-10)

        for i, t in enumerate(tt):
            data4d_sinc[i, ...] = state.m
            llg_sinc.step(state, dt)
        torch.save({"data4d_sinc":data4d_sinc}, "data/sinc.pt")

    delta_m_sinc = data4d_sinc - m0[None, ...]
    m_fft_sinc = torch.fft.rfft(delta_m_sinc, dim=0)
    H0 = 1 / (2 * f_max * dt) # m_fft includes H0 factor from sinc spectrum, so divide by H0² to get |δm̂|²
    power_sinc = (state.material["Ms"]**2 * torch.abs(m_fft_sinc)**2).mean(axis=(1,2,3)).sum(axis=-1) / H0**2 / 2
    peaks = scipy.signal.find_peaks(power_sinc.numpy())[0]


### calculate eigenmodes and spectrum
with Timer("Caculate Eigenmodes"):
    try:
        res = EigenResult.load(state, "data/eigen.pt")
    except Exception:
        state.m = m0
        eigen = EigenSolver(state, [exchange_b, exchange_t, rkky, aniso, dmi], [bias])
        res = eigen.solve(k=20)
        res.store("data/eigen.pt")
        res.save_evecs3D("data/evecs.pvd")

    spectrum = res.spectrum(2*np.pi*freq, excite.h(state))
    projection = res.projection(2*np.pi*freq, delta_m)


### plot results
fig, ax = plt.subplots(figsize=(15,10))
ax.plot(freq * 1e-9, spectrum, "--", linewidth=2.0, label="Eigenmode (spectrum)")
ax.plot(freq * 1e-9, projection, "-.", linewidth=2.0, label="Eigenmode (projection)")
ax.plot(freq * 1e-9, power_sinc, "k--", linewidth=2.0, label="Time-Domain (sinc excitation)")
ax.scatter(freq[peaks] * 1e-9, power_sinc[peaks], color="red")

ax.set_xlim([0, 50])
ax.set_ylim([1e5, 1e9])
ax.set_yscale("log")
ax.set_xlabel("Frequency [GHz]")
ax.set_ylabel("Volume-Averaged Magnetization Power Spectrum\n" + r"$P(\omega) = \frac{1}{2} \langle M_s^2 |\delta\hat{m}|^2 \rangle$ [A$^2$/m$^2$]")

freq_eig = res.freq * 1e-9
tick_labels = [f"{f:.5f}" for f in freq_eig]
for p in peaks[:12]:
    x_val = freq[p] * 1e-9 # GHz
    y_val = power_sinc[p]
    ax.text(x_val, y_val, f"{freq[p]*1e-9:.1f}", rotation=45, ha='left',va='bottom')
ax.set_xticks(np.arange(0,50,5))
ax.tick_params(axis='both', direction='in', length=6, width=1.2)
ax.grid()
ax.legend(loc='upper right')

fig.savefig("data/result_bilayer.png")
Timer.print_report()
