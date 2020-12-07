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

# initialize anisotropy field
aniso    = AnisotropyField(mesh, material)

# initialize random magnetization
m0 = np.random.rand(n[0], n[1], n[2], 3) - 0.5
m0 = m0 / np.linalg.norm(m0, axis=3).repeat(3).reshape(m0.shape)

# minimize energy
write_vtr(m0, "m_start", mesh)
minimizer = Minimizer([aniso])
m = minimizer.minimize(m0, 1e-4, 1e-4)
write_vtr(m, "m_relaxed", mesh)
