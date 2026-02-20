# %% [markdown]
# ## Plot Results

# %%
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

fig, ax = plt.subplots(figsize=(10, 10))
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

sigma2 = np.array([0, -0.1, -0.05, -0.01])


for i in range(1):
    sigma2_ = sigma2[i]
    
    ref_Bz = np.load(f"ref/oommf_hysteresis_Bz_sigma2_{sigma2_}.npy")/1000
    ref_Mz = np.load(f"ref/oommf_hysteresis_mz_sigma2_{sigma2_}.npy")

    data_Bz = np.load(f"data/hysteresis_Bz_sigma2_{sigma2_*2}.npy")
    data_Mz = np.load(f"data/hysteresis_Mz_sigma2_{sigma2_*2}.npy")

    ax.plot(ref_Bz, ref_Mz, '-', color = cycle[i], label = f"magnum.np (sigma2 = {sigma2_*2})")
    ax.plot(data_Bz, data_Mz, '-', color = cycle[i], linewidth = 6, alpha = 0.4, label = f"reference (sigma2 = {sigma2_})")



ax.set_xlabel("Bz (T)")
ax.set_ylabel("Mz")
ax.legend(ncol = 2)
ax.grid()
fig.savefig("data/results.png", dpi=300)
