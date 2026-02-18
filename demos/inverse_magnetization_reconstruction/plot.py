# %% [markdown]
# ## Plot Results

# %%
import torch
import matplotlib.pyplot as plt

data = torch.load("data/results.pt")
m_true = data["m_true"].cpu()
m_0 = data["m_0"].cpu()
m_rec_1 = data["m_rec_1"].cpu()
m_rec = data["m_rec"].cpu()
h_true = data["h_true"].cpu()
h_0 = data["h_0"].cpu()
h_rec_1 = data["h_rec_1"].cpu()
h_rec = data["h_rec"].cpu()

# %% [markdown]
# ### Magnetization Comparison (Ground Truth / Starting Guess / After 1 Step / Final)

# %%
comp_labels = ["$m_x$", "$m_y$", "$m_z$", r"$\mathbf{H}^\mathrm{dem} \cdot \mathbf{n}_\mathrm{NV}$ (A/m)"]
row_labels = [
    "Ground Truth",
    "Starting Guess",
    "After 1 Step",
    "Final",
]
row_data = [m_true, m_0, m_rec_1, m_rec]
row_fields = [h_true, h_0, h_rec_1, h_rec]

fig, axes = plt.subplots(4, 4, figsize=(20, 16))

for row, (label, m, h) in enumerate(zip(row_labels, row_data, row_fields)):
    # Plot magnetization components
    for i in range(3):
        im = axes[row, i].imshow(m[:, :, 0, i].T, origin="lower", cmap="RdBu_r", vmin=-1, vmax=1)
        # Only show component labels in the top row
        if row == 0:
            axes[row, i].set_title(comp_labels[i], fontsize=12)
        # Add row labels on the left
        if i == 0:
            axes[row, i].set_ylabel(label, fontsize=11, fontweight='bold')
        axes[row, i].set_xticks([])
        axes[row, i].set_yticks([])
        cb = plt.colorbar(im, ax=axes[row, i], fraction=0.046)
        cb.set_ticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    
    # Plot H-field projection
    h_vlim = h.abs().max().item()
    im = axes[row, 3].imshow(h.T, origin="lower", cmap="coolwarm", vmin=-h_vlim, vmax=h_vlim)
    if row == 0:
        axes[row, 3].set_title(comp_labels[3], fontsize=12)
    axes[row, 3].set_xticks([])
    axes[row, 3].set_yticks([])
    plt.colorbar(im, ax=axes[row, 3], fraction=0.046)

fig.suptitle("Magnetization Reconstruction Results", fontsize=16, fontweight='bold')
fig.tight_layout()
fig.savefig("data/results.png", dpi=150)

# %% [markdown]
# ### Pointwise Error

# %%
m_err = (m_rec - m_true).norm(dim=-1)[:, :, 0]
fig2, ax = plt.subplots(figsize=(6, 5))
im = ax.imshow(m_err.T, origin="lower", cmap="hot")
ax.set_title("Pointwise $\|\\mathbf{m}_{rec} - \\mathbf{m}_{true}\|$")
ax.set_xticks([])
ax.set_yticks([])
plt.colorbar(im, ax=ax, fraction=0.046)
fig2.tight_layout()
fig2.savefig("data/error_map.png", dpi=150)
