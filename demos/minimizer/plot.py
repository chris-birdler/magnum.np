# %% [markdown]
# ## Plot Results

# %%
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data/log.dat")

fig, ax = plt.subplots(figsize=(15,5))
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

ax.plot(4*np.pi*1e-7*1e3*data[:,2], data[:,5]*7./4., '-', color = cycle[2], label = "magnum.np")
ax.axvline(x = -2646, linestyle = '--', color = 'black', alpha = 0.2)
ax.axvline(x = -2753, linestyle = '--', color = 'black', alpha = 0.2)
ax.axvline(x = -2845, linestyle = '--', color = 'black', alpha = 0.2)
ax.axvline(x = -2848, linestyle = '--', color = 'black', alpha = 0.2)

ax.set_title("Minimizer Standard Problem")
ax.set_xlim([-2900,-2600])
ax.set_xlabel("External Field $\mu_0 H_{ext}$ [mT]")
ax.set_ylabel("Magnetization m")
ax.grid()
fig.savefig("data/results.png")
