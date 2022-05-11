from fd import *
import numpy as np
import torch
import os

CUDA_DEVICE = os.environ.get('CUDA_DEVICE', '0')
cuda = torch.device(f"cuda:{CUDA_DEVICE}" if torch.cuda.is_available() else "cpu")

# initialize mesh
N = 31
n = (N, N, N)
dx = (10e-9, 10e-9, 10e-9)
mesh = Mesh(n, dx)

# initialize field terms
oersted = OerstedField(mesh)

# initialize current
j = torch.zeros(n + (3,), dtype=torch.float64, device=cuda)
j[N//2,N//2,:,0] = 0.
j[N//2,N//2,:,1] = 0.
j[N//2,N//2,:,2] = 1.

h = oersted.h(0., j)
write_vtr(j, "data/j", mesh)
write_vtr(h, "data/h", mesh)

x = np.arange(-N//2+1,N//2+1)*dx[0]
h_sim = h[N//2,:,N//2,0]
h_analytic = -dx[0]*dx[1]/(2.*np.pi*x)

np.savetxt("data/result.txt", np.vstack((x, h_sim, h_analytic)).T)
