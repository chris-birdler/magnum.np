# %% [markdown]
# ## Plot Results

# %%
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data/m.dat")
ref = np.loadtxt("ref/m_ref.dat")

fig, ax = plt.subplots(figsize=(15,5))
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

ax.plot(data[:,1]*4*np.pi*1e-7, data[:,4], '-', color = cycle[2], label = "magnum.np")
ax.plot(ref[:,1]*4*np.pi*1e-7, ref[:,4], '-', color = cycle[2], linewidth = 6, alpha = 0.4, label = "magnum.af")

ax.set_xlim([-0.1,0.1])
ax.set_ylim([-1,1])
ax.set_title("SMC Demo")
ax.set_xlabel("External Field $\mu$$_0$ H$^{ext}$$_x$ [T]")
ax.set_ylabel("Reduced Magnetization m [1]")
ax.legend(ncol=3)
ax.grid()
fig.savefig("data/results.png")
