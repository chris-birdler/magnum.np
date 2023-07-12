import numpy as np
import matplotlib.pyplot as plt

magnumnp = np.loadtxt("data/m0_magnumnp.dat")
#magnumpi = np.loadtxt("data/m0_magnumpi.dat")
finmag = np.loadtxt("data/m0_finmag.dat")

ref_x = np.loadtxt("ref/interfaceDMI_mx.dat")
ref_z = np.loadtxt("ref/interfaceDMI_mz.dat")

cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

fig, (ax1,ax2) = plt.subplots(nrows=2, sharex=True)
ax1.plot(ref_x[:,0], ref_x[:,1], "--", color=cycle[0], alpha = 0.7, label = r"reference")
ax1.plot(magnumnp[1:-1,0]*1e9, magnumnp[1:-1,1], "x", color=cycle[1], alpha = 0.7, label = "magnum.np")
ax1.plot(finmag[::2,0]-50., finmag[::2,1], "o", mfc='none', color=cycle[2], alpha = 0.7,  label = "finmag")
ax1.set_ylabel("Magnetization $m_x$")

ax2.plot(ref_z[:, 0], ref_z[:,1], "--", color=cycle[0], alpha = 0.7, label = r"reference")
ax2.plot(magnumnp[1:-1:2,0]*1e9, magnumnp[1:-1:2,3], "x", color=cycle[1], alpha = 0.7, label = "magnum.np")
ax2.plot(finmag[::2,0]-50., finmag[::2,3], "o", mfc='none', color=cycle[2], alpha = 0.7, label = "finmag")
ax2.set_ylabel("Magnetization $m_z$")

for ax in ax1, ax2:
    ax.legend()
    ax.grid()
ax2.set_xlabel("Position x [nm]")
fig.savefig("data/results.png")