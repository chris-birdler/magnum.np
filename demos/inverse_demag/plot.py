# %% [markdown]
# ## Plot Results

# %%
import torch
import numpy as np
import matplotlib.pyplot as plt

fig, (ax1,ax2) = plt.subplots(1,2,figsize=(15,5))
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

n, dx, grad = torch.load("data/grad.pt")
n = n.cpu().numpy()
dx = dx.cpu().numpy()
grad = torch.log(grad[:,:,0,0]).cpu().numpy()
ref_grad = torch.log(torch.load("ref/grad.pt")[:,:,0,0]).cpu().numpy()
x = np.linspace(-n[0]*dx[0]/2.*1e9, n[0]*dx[0]/2.*1e9, num = 200)

ax1.imshow(grad, origin = "lower", extent=(-n[0]*dx[0]/2.*1e9,n[0]*dx[0]/2.*1e9,0*1e9,n[1]*dx[1]*1e9))
ax1.plot(x, np.abs(x)/np.sqrt(2), "k--")
ax1.set_title("Current Simulation")

ax2.imshow(ref_grad, origin = "lower", extent=(-n[0]*dx[0]/2.*1e9,n[0]*dx[0]/2.*1e9,0*1e9,n[1]*dx[1]*1e9))
ax2.plot(x, np.abs(x)/np.sqrt(2), "k--")
ax2.set_title("Reference Simulation")

fig.savefig("data/results.png")
