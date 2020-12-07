from fd import *
import numpy as np

# initialize mesh
n  = (5, 7, 9)
dx = (10e-9, 10e-9, 10e-9)
mesh = Mesh(n, dx)

# initialize material
material = {
        "Ms": 8e5,
        }

# initialize field terms
demag = DemagFieldPBC(mesh, material)

# initialize magnetization that relaxes into s-state
m = np.zeros(n + (3,))
m[:,:,:,:] = [1, 0, 0]

h = demag.h(0., m)
write_vtr(h, "data/h", mesh)
