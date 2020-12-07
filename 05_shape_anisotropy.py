from fd import *
import numpy as np
from scipy import constants

# initialize mesh
n  = (10, 10, 1)
dx = (10e-9, 10e-9, 2e-9)
mesh = Mesh(n, dx)

# initialize material
material = {
        "Ms": 8e5,
        "K": - constants.mu_0 * 8e5**2 / 2.,
        "K_axis": (0,0,1)
        }

# initialize field terms
demag = DemagField(mesh, material)
aniso = AnisotropyField(mesh, material)

m = np.zeros(n + (3,))
for theta in range(181):
    m[:,:,:,:] = [0, np.sin(np.radians(theta)), np.cos(np.radians(theta))]
    print "%d %g %g" % (theta, demag.E(0, m), aniso.E(0, m))
