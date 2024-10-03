import torch
import numpy as np
import matplotlib.pyplot as plt

mt = torch.load("data/mt_saved.pt").cpu()
k = 0.2
f = 100

dt = 1e-12
n = 500, 25, 1
dx = 2e-9, 2e-9, 1e-9

# plot dispersion
mfft = torch.fft.fftn(mt, dim = (0,1))
ff = torch.fft.fftfreq(mfft.shape[0], dt)
kk = torch.fft.fftfreq(mfft.shape[1], dx[0])*2*torch.pi

disp = mfft[:,:,int(n[1]//2),0].abs()
cbmax_bv = 10.*np.log10(disp**2).max()
cbmin_bv = 10.*np.log10(disp**2).min()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize = [10, 4])
ax1.set_title('Dispersion')
ax1.pcolormesh(torch.fft.fftshift(kk)/1e9, torch.fft.fftshift(ff)/1e9, 10.*np.log(torch.fft.fftshift(disp)**2.), vmin = cbmin_bv, vmax = cbmax_bv, cmap = "viridis", shading = "auto")
dot, = ax1.plot(k, f, 'r.')
ax1.set_xlabel(r'$k_{x}\,\mathrm{(rad/nm)}$', fontsize=10)
ax1.set_ylabel('f (GHz)', fontsize=10)
ax1.set_ylim([0, 500])
ax1.set_xlim([-0.5, 0.5])

# plot mode profile
def get_mi(f,k):
    # Determine ki and fi indizes
    fi = (ff-f*1e9).abs().argmin()
    ki = (kk-k*1e9).abs().argmin()

    mfft_i = torch.zeros_like(mfft[fi,...])
    mfft_i[ki,...] = mfft[fi,ki,...]
    mi = torch.fft.ifftn(mfft_i, dim = 0)
    return mi[:,:,0].real

mi = get_mi(f,k)
mi = mi / mi.max()

ax2.set_title('Mode Profile')
mode = ax2.imshow(mi)
ax2.axes.xaxis.set_visible(False)
ax2.axes.yaxis.set_visible(False)

fig.canvas.header_visible = False
fig.canvas.footer_visible = False

fig.tight_layout()
fig.savefig("data/result.png")
