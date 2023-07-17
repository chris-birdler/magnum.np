from magnumnp import *
import matplotlib.pyplot as plt
import numpy as np

Timer.enable()

lex = 5.71e-9
Js = 1.
n  = (40, 40, 1)
dx = (2.5e-9, 2.5e-9, 3e-9)

mesh = Mesh(n, dx)
state = State(mesh)
state.material = {
        "A":lex**2*Js**2/(2.*constants.mu_0),
        "Ms":1./constants.mu_0,
        }
demag    = DemagField()
exchange = ExchangeField()
external = ExternalField([0.,0.,1.2/constants.mu_0])

state.m = state.Constant([0.,0.,1.])
eigen = EigenSolver(state, [demag, exchange], [external])
res = eigen.solve(k=20, tol=1e-6)
print("evals[GHz]:", res.omega.numpy()/2./torch.pi*1e-9)
res.save_evecs3D("data/evecs.vti")

# plot the results
ref_magnumpi = np.loadtxt("ref/thinfilm_magnumpi.dat")
ref_daquino = np.loadtxt("ref/thinfilm_daquino.dat")

fig, ax = plt.subplots()
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

ax.plot(ref_daquino, '^--', fillstyle = "none", color = "black", alpha = 0.7, linewidth=0.8, label = "d'Aquino (FD)")
ax.plot(ref_magnumpi, 'o', mfc = 'none', color = cycle[1], label = "magnum.pi (FE)")
ax.plot(res.omega.numpy()/2./torch.pi*1e-9, '+', color = cycle[0], label = "magnum.np (FD)")

ax.set_xlim([0,4])
ax.set_ylim([7.5,15])
ax.set_title("Saturated Thin-Film")
ax.set_xlabel("Eigenmode Index [1]")
ax.set_ylabel("Frequency f [GHz]")
ax.legend(loc=2)
ax.grid()
fig.savefig("data/results_thinfilm.png")

Timer.print_report()
