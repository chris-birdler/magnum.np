import numpy as np
import matplotlib.pyplot as plt
import itertools

with open("timings_eval.dat", 'r') as file:
    names = file.readline().strip().split()

### plot eval
data = np.loadtxt("timings_eval.dat")
ref = np.loadtxt("ref/timings_eval.dat")
ref_mumax = np.loadtxt("ref/timings_eval_mumax.dat")

fig, ax = plt.subplots(4,1,figsize=(15, 12))
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
for axis in ax:
    axis.set_prop_cycle(marker=['o', '+', 'x', '*', '<', '>', '^', 'v', 'X', 's', 'd', 'D'])

for n, i in enumerate([1, 2, 3]): #  'DemagField.h', 'DemagFieldPBC.h', 'OerstedField.h'
    ax[0].loglog(data[:,0]**3, data[:,i], color = cycle[n], label = names[i])
    ax[0].loglog(ref[:,0]**3, ref[:,i], '-', color = cycle[n], linewidth = 6, alpha = 0.4)
ax[0].loglog(ref_mumax[:,0]**3, ref_mumax[:,1]/10000, "--", color = "black", label = "DemagField (MuMax3)")

for n, i in enumerate([4, 7, 8, 9]): # 'ExchangeField.h', 'InterfaceDMIField.h', 'BulkDMIField.h', 'D2dDMIField.h'
    ax[1].loglog(data[:,0]**3, data[:,i], color = cycle[n], label = names[i])
    ax[1].loglog(ref[:,0]**3, ref[:,i], '-', color = cycle[n], linewidth = 6, alpha = 0.2)

for n, i in enumerate([11, 12, 13]): # 'SpinTorqueSlonczewski.h', 'SpinOrbitTorque.h', 'SpinTorqueZhangLi.h'
    ax[2].loglog(data[:,0]**3, data[:,i], color = cycle[n], label = names[i])
    ax[2].loglog(ref[:,0]**3, ref[:,i], '-', color = cycle[n], linewidth = 6, alpha = 0.4)

for n, i in enumerate([5, 6, 10]): # 'UniaxialAnisotropyField.h', 'CubicAnisotropyField.h', 'ExternalField.h'
    ax[3].loglog(data[:,0]**3, data[:,i], color = cycle[n], label = names[i])
    ax[3].loglog(ref[:,0]**3, ref[:,i], '-', color = cycle[n], linewidth = 6, alpha = 0.4)

ax[-1].set_xlabel("Number of Cells N")
for axis in ax:
    axis.set_ylabel(r"$t_{eval}$ [s]")
    axis.legend(ncol=3)
    axis.grid()
fig.savefig("timings_eval.png")



### plot setup
data = np.loadtxt("timings_setup.dat")
ref = np.loadtxt("ref/timings_setup.dat")

fig, ax = plt.subplots(4,1,figsize=(15, 12))
for axis in ax:
    axis.set_prop_cycle(marker=['o', '+', 'x', '*', '<', '>', '^', 'v', 'X', 's', 'd', 'D'])

for n, i in enumerate([1, 2, 3]): #  'DemagField.h', 'DemagFieldPBC.h', 'OerstedField.h'
    ax[0].loglog(data[:,0]**3, data[:,i], color = cycle[n], label = names[i])
    ax[0].loglog(ref[:,0]**3, ref[:,i], '-', color = cycle[n], linewidth = 6, alpha = 0.4)

for n, i in enumerate([4, 7, 8, 9]): # 'ExchangeField.h', 'InterfaceDMIField.h', 'BulkDMIField.h', 'D2dDMIField.h'
    ax[1].loglog(data[:,0]**3, data[:,i], color = cycle[n], label = names[i])
    ax[1].loglog(ref[:,0]**3, ref[:,i], '-', color = cycle[n], linewidth = 6, alpha = 0.2)

for n, i in enumerate([11, 12, 13]): # 'SpinTorqueSlonczewski.h', 'SpinOrbitTorque.h', 'SpinTorqueZhangLi.h'
    ax[2].loglog(data[:,0]**3, data[:,i], color = cycle[n], label = names[i])
    ax[2].loglog(ref[:,0]**3, ref[:,i], '-', color = cycle[n], linewidth = 6, alpha = 0.4)

for n, i in enumerate([5, 6, 10]): # 'UniaxialAnisotropyField.h', 'CubicAnisotropyField.h', 'ExternalField.h'
    ax[3].loglog(data[:,0]**3, data[:,i], color = cycle[n], label = names[i])
    ax[3].loglog(ref[:,0]**3, ref[:,i], '-', color = cycle[n], linewidth = 6, alpha = 0.4)

ax[-1].set_xlabel("Number of Cells N")
for axis in ax:
    axis.set_ylabel(r"$t_{setup}$ [s]")
    axis.legend(ncol=3)
    axis.grid()
fig.savefig("timings_setup.png")
