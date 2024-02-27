import numpy as np
import matplotlib.pyplot as plt
import itertools

with open("timings_eval.dat", 'r') as file:
    names = file.readline().strip().split()

### plot eval
data = np.loadtxt("timings_eval.dat")
ref = np.loadtxt("ref/timings_eval.dat")

fig, ax = plt.subplots(3,1,figsize=(15, 12))
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
for axis in ax:
    axis.set_prop_cycle(marker=['o', '+', 'x', '*', '<', '>', '^', 'v', 'X', 's', 'd', 'D'])

for i in range(1, data.shape[1]):
    ax[i%3].loglog(data[:,0]**3, data[:,i], color = cycle[(i-1)//3], label = names[i])
for i in range(1, ref.shape[1]):
    ax[i%3].loglog(ref[:,0]**3, ref[:,i], '-', color = cycle[(i-1)//3], linewidth = 6, alpha = 0.4)

for axis in ax:
    axis.set_xlabel("Number of Cells N")
    axis.set_ylabel(r"$t_{eval}$ [s]")
    axis.legend(ncol=3)
    axis.grid()
fig.savefig("timings_eval.png")



### plot setup
data = np.loadtxt("timings_setup.dat")
ref = np.loadtxt("ref/timings_setup.dat")

fig, ax = plt.subplots(3,1,figsize=(15, 12))
for axis in ax:
    axis.set_prop_cycle(marker=['o', '+', 'x', '*', '<', '>', '^', 'v', 'X', 's', 'd', 'D'])

for i in range(1, data.shape[1]):
    ax[i%3].loglog(data[:,0]**3, data[:,i], color = cycle[(i-1)//3], label = names[i])
for i in range(1, ref.shape[1]):
    ax[i%3].loglog(ref[:,0]**3, ref[:,i], '-', color = cycle[(i-1)//3], linewidth = 6, alpha = 0.4)

for axis in ax:
    axis.set_xlabel("Number of Cells N")
    axis.set_ylabel(r"$t_{setup}$ [s]")
    axis.legend(ncol=3)
    axis.grid()
fig.savefig("timings_setup.png")
