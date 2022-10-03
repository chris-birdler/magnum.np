import numpy as np
import matplotlib.pyplot as plt

magnumnp = np.loadtxt("data/m0_magnumnp.dat")
#magnumpi = np.loadtxt("data/m0_magnumpi.dat")
finmag = np.loadtxt("data/m0_finmag.dat")

ref_x = np.loadtxt("data/interfaceDMI_mx.dat")
ref_z = np.loadtxt("data/interfaceDMI_mz.dat")

fig, (ax1,ax2) = plt.subplots(nrows=2, sharex=True)
ax1.plot(ref_x[:,0], ref_x[:,1], label = r"reference")
ax1.plot(magnumnp[1:-1,0]*1e9, magnumnp[1:-1,1], label = "magnum.np")
#ax1.plot(magnumpi[:,0]-50., magnumpi[:,1], label = "magnum.pi")
ax1.plot(finmag[:,0]-50., finmag[:,1], label = "finmag")
ax1.set_ylabel("Magnetization m_x")

ax2.plot(ref_z[:, 0], ref_z[:,1], label = r"reference")
ax2.plot(magnumnp[1:-1,0]*1e9, magnumnp[1:-1,3], label = "magnum.np")
#ax2.plot(magnumpi[:,0]-50., magnumpi[:,3], label = "magnum.pi")
ax2.plot(finmag[:,0]-50., finmag[:,3], label = "finmag")
ax2.set_ylabel("Magnetization m_z")

for ax in ax1, ax2:
    ax.legend()
    ax.grid()
ax2.set_xlabel("Position x [nm]")
fig.savefig("data/results.png")
