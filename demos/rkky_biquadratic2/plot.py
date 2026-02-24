# %% [markdown]
# ## Plot Results

# %%
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

fig, ax = plt.subplots(figsize=(10, 10))
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

for i, sigma in enumerate(np.array([0, -0.1, -0.05, -0.01])):
    ref_Bz = np.load(f"ref/oommf_hysteresis_Bz_sigma2_{sigma}.npy")/1000
    ref_Mz = np.load(f"ref/oommf_hysteresis_mz_sigma2_{sigma}.npy")

    data = np.loadtxt(f"data/m_sigma2_{sigma*2}.dat")
    B = data[:,3]*4*np.pi*1e-7
    M = data[:,9]/2 + data[:,12]/2

    ax.plot(ref_Bz, ref_Mz, '-', color = cycle[i], linewidth = 6, alpha = 0.4, label = f"reference (sigma2 = {sigma*2})")
    ax.plot(B, M, '-', color = cycle[i], label = f"magnum.np (sigma2 = {sigma * 2})")

ax.set_xlabel("Bz (T)")
ax.set_ylabel("Mz")
ax.legend()
ax.grid()
fig.savefig("data/results.png", dpi=300)
