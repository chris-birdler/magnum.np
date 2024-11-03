# %% [markdown]
# ## Plot Results

# %%
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(15,5))
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

for i, xi in enumerate([30, 91, 242, 725]):
    try:
        data = np.loadtxt("data/log_xi%d.dat" % xi)
        h = 4*np.pi*1e-7*data[:,1]
        ax.plot(h, data[:,4], 'x', color = cycle[i], label = f"$\\xi = {xi}$")
    except:
        pass

    h_ref = np.linspace(0.0001, 0.1, num=100)
    ax.plot(h_ref, np.cosh(xi*h_ref)/np.sinh(xi*h_ref)-1/(xi*h_ref), '-', color = cycle[i])


ax.set_xlim([0,0.1])
ax.set_ylim([0,1])
ax.set_xlabel("External Field $\mu_0 H^{ext}_z [T]$")
ax.set_ylabel("Average Magnetization $m_z$")
ax.legend()
ax.grid()
fig.savefig("data/results.png")
