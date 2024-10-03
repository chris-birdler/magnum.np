import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data/log.dat")
ref = np.loadtxt("data/ref.dat")

fig, ax = plt.subplots(figsize=(10, 5))
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
ax.semilogx(data[:,0], data[:,1]*1e9, 'k--', label = "analytic", alpha = 0.5)

ax.semilogx(data[:,0], data[:,2]*1e9, '-', color = cycle[0], label = "magnum.np (int B dA)")
ax.semilogx(ref[:,0], ref[:,2]*1e9, '-', color = cycle[0], linewidth = 6, alpha = 0.4)

ax.semilogx(data[:,0], data[:,3]*1e9, '-', color = cycle[1], label = "magnum.np (int j*A dV)")
ax.semilogx(ref[:,0], ref[:,3]*1e9, '-', color = cycle[1], linewidth = 6, alpha = 0.4)

ax.set_xlabel("Aspect Ratio $W/H$")
ax.set_ylabel("Self-Inductivity $L[nH]$")
ax.legend()
ax.grid()
fig.savefig("data/results.png")
