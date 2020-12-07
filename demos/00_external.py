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
        "K": 1e4,
        "K_axis": (0,0,1)
        }

# initialize external field pointing in y-direction
external = ExternalField(mesh, material, (0,1,0))

# initialize magnetization pointing in x=direction
m0 = np.zeros(n + (3,))
m0[:,:,:,0] = 1.0

# minimize energy
write_vtr(m0, "m_start", mesh)
minimizer = Minimizer([external])
m = minimizer.minimize(m0, 3e-2, 1e-4)
write_vtr(m, "m_relaxed", mesh)
