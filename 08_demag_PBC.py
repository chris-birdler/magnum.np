from fd import *
import numpy as np

# initialize mesh
n  = (10, 10, 10)
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
m[3:7,3:7,3:7,:] = [-1, 0, 0]
#m[4:6,:,:,:] = [-1, 0, 0]

h, u, div = demag.h(0., m)
write_vtr(m, "data/m", mesh)
write_vtr(h, "data/h", mesh)
write_vtr(u, "data/u", mesh)
write_vtr(div, "data/div", mesh)
