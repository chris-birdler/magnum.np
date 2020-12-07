from fd import *
import numpy as np

# initialize mesh
n  = (100, 25, 1)
dx = (5e-9, 5e-9, 3e-9)
mesh = Mesh(n, dx)

# initialize material
material = {
        "Ms": 8e5,
        "A": 1.3e-11,
        "K": 1e5,
        "K_axis": (0,0,1)
        }

# initialize anisotropy and exchange field
aniso    = AnisotropyField(mesh, material)
exchange = ExchangeField(mesh, material)

# initialize magnetization that relaxe
m0 = np.zeros(n + (3,))
m0[:50,:,:,2] = +1.0
m0[50:,:,:,2] = -1.0
m0[50,:,:,:]  = [0.0, 1.0, 0.0]

# minimize energy
write_vtr(m0, "m_start", mesh)
minimizer = Minimizer([aniso, exchange])
m = minimizer.minimize(m0, 1e-6, 1e-4)
write_vtr(m, "m_relaxed", mesh)
