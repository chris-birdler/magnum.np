import numpy as np
import matplotlib.pyplot as plt
import itertools

### plot eval
data = np.loadtxt("timings.dat")
ref = np.loadtxt("ref/timings.dat")
mumax = np.loadtxt("ref/mumax.dat")

fig, ax = plt.subplots(figsize=(15, 5))
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

ax.loglog(data[:,0], data[:,3]/1e6, color = cycle[0], label = "magnum.np")
ax.loglog(data[:,0], ref[:,3]/1e6, '-', color = cycle[0], linewidth = 6, alpha = 0.4, label = "magnum.np(ref)")
ax.loglog(mumax[:,0], mumax[:,1]/1e6, 'k--', label = "MuMax3")

ax.set_xlabel("Number of Cells N")
ax.set_ylabel(r"Throughput [million cells / s]")
ax.legend()
ax.grid()
fig.savefig("timings.png")
