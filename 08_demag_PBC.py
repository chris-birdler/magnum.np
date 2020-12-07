from .. import *
import numpy as np

# initialize mesh
n  = (10, 10, 1)
dx = (10e-9, 10e-9, 2e-9)
mesh = Mesh(n, dx)

# initialize material
material = {
        "Ms": 8e5,
        }

# initialize field terms
demag = DemagField(mesh, material)

# initialize magnetization that relaxes into s-state
m = np.zeros(n + (3,))
m[:,:,:,:] = [0, 0, 1]

h = demag.h(0., m)
write_vtr(h, "data/h", mesh)
