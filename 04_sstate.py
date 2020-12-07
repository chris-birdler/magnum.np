from fd import *
import numpy as np

# initialize mesh
n  = (100, 25, 1)
dx = (5e-9, 5e-9, 3e-9)
mesh = Mesh(n, dx)

# initialize material
material = {
        "Ms": 8e5,
        "A": 1.3e-11
        }

# initialize field terms
demag    = DemagField(mesh, material)
exchange = ExchangeField(mesh, material)

# initialize magnetization that relaxes into s-state
m0 = np.zeros(n + (3,))
m0[:,:,:,:] = [1./np.sqrt(2), 0, 1./np.sqrt(2.)]

# minimize energy
write_vtr(m0, "m_start", mesh)
minimizer = Minimizer([demag, exchange])
m = minimizer.minimize(m0, 1e-3, 1e-4)
write_vtr(m, "m_relaxed", mesh)
