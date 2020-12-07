from fd import *
import numpy as np
from scipy import constants

# initialize mesh
n  = (100, 25, 1)
dx = (5e-9, 5e-9, 3e-9)
mesh = Mesh(n, dx)

# initialize material
material = {
        "Ms": 8e5,
        "A": 1.3e-11,
        "gamma": 2.211e5,
        "alpha": 0.02
        }

# initialize field terms
demag    = DemagField(mesh, material)
exchange = ExchangeField(mesh, material)
external = ExternalField(mesh, material, [-24.6e-3/constants.mu_0,
                                          +4.3e-3/constants.mu_0,
                                          0.0])

# initialize magnetization that relaxes into s-state
m0 = np.zeros(n + (3,))
m0[1:-1,:,:,0]   = 1.0
m0[(-1,0),:,:,1] = 1.0

# initialize sstate
minimizer = Minimizer([demag, exchange])
m = minimizer.minimize(m0, 1e-2, 1e-4)

# perform integration with external field
llg = LLG([demag, exchange, external], material, m)
def E(m):
    return demag.E(material, m) + exchange.E(material, m) + external.E(material, m)

i = 0
with open('sp4.dat', 'w') as f:
    while llg.t < 1e-9:
        if i % 10 == 0:
            write_vtr(m, "sp4/m_%04d" % (i/10))
        f.write("%g %g %g %g %g\n" % ((llg.t,) + tuple(np.mean(m, axis=(0,1,2))) + (E(m),)))

        m = llg.step(1e-12)
        i += 1
