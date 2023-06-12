import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data/log.dat")
data_ref = np.loadtxt("ref/m_ref.dat")

t, mx, my, mz = data.T
t_ref, mx_ref, my_ref, mz_ref = data_ref.T

sy = abs((np.fft.fft(my)**2.))
sy_ref = abs((np.fft.fft(my_ref)**2.))
f = np.fft.fftfreq(len(t), t[1]-t[0])*1e-9
f_ref = np.fft.fftfreq(len(t_ref), t_ref[1]-t_ref[0])*1e-9


fig, (ax0, ax1) = plt.subplots(2,1)
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

ax0.plot(t*1e9, my, color=cycle[0])
ax0.plot(t_ref*1e9, my_ref, '-', color=cycle[0], linewidth=6, alpha=0.4)
ax0.set_xlabel("Time (ns)")
ax0.set_ylabel("$m_y (A/m)$")
ax0.set_xlim(0, 2.5)
ax0.set_ylim(0.580, 0.594)
ax0.grid()

ax1.plot(f[2: -len(f)//2], sy[2: -len(f)//2], '-', color=cycle[1])
ax1.plot(f_ref[2: -len(f_ref)//2], sy_ref[2: -len(f_ref)//2], '-', color=cycle[1], linewidth=6, alpha=0.4)
ax1.set_yscale("log")
ax1.set_xlabel("Frequency (GHz)")
ax1.set_ylabel("$m_y (A/m)$")
ax1.set_xlim(0, 20)
ax1.grid()
fig.tight_layout()
fig.savefig("data/results.png")