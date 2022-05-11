from fd import *
import numpy as np
from scipy import constants

# initialize mesh
n  = (40, 40, 1)
dx = (2.5e-9, 2.5e-9, 10e-9)
mesh = Mesh(n, dx)

# initialize material
material = {
        "Ms": 8e5,
        "A": 1.3e-11,
        "gamma": 2.211e5,
        "alpha": 0.1,
        "xi": 0.05,
        "b": 72.17e-12
        }

# initialize field terms
demag    = DemagField(mesh, material)
exchange = ExchangeField(mesh, material)
torque   = SpinTorque(mesh, material, [1e12, 0, 0])

# initialize magnetization that relaxes into s-state
m0 = np.zeros(n + (3,))
m0[:20,:,:,1]   = -1.0
m0[20:,:,:,1]   = 1.0
m0[20,20,:,1]   = 0.
m0[20,20,:,2]   = 1.


# initialize sstate
minimizer = Minimizer([demag, exchange])
m = minimizer.minimize(m0, 1e-2, 1e-4)
#m = m0

# perform integration with external field
llg = LLG([demag, exchange, torque], material, m)
i = 0
with open('sp5.dat', 'w') as f:
    while llg.t < 10e-9:
        if i % 10 == 0:
            write_vtr(m, "sp5/m_%04d" % (i/10))
        f.write("%g %g %g %g\n" % ((llg.t,) + tuple(np.mean(m, axis=(0,1,2)))))

        m = llg.step(1e-12)
        i += 1
