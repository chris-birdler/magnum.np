from fd import *
import numpy as np

# initialize mesh
N = 31
n = (N, N, N)
dx = (10e-9, 10e-9, 10e-9)
mesh = Mesh(n, dx)

# initialize field terms
oersted = OerstedField(mesh)

# initialize current
j = np.zeros(n + (3,))
j[N//2,N//2,:,:] = [0, 0, 1]

h = oersted.h(0., j)
write_vtr(j, "data/j", mesh)
write_vtr(h, "data/h", mesh)

x = np.arange(-N//2+1,N//2+1)*dx[0]
h_sim = h[N//2,:,N//2,0]
h_analytic = -dx[0]*dx[1]/(2.*np.pi*x)

np.savetxt("data/result.txt", np.vstack((x, h_sim, h_analytic)).T)
