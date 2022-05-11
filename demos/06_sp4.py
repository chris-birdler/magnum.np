from magnumnp import *
from scipy import constants
import torch
import os

CUDA_DEVICE = os.environ.get('CUDA_DEVICE', '0')
cuda = torch.device(f"cuda:{CUDA_DEVICE}" if torch.cuda.is_available() else "cpu")

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
m0 = torch.zeros(n + (3,), dtype=torch.float64, device=cuda)
m0[1:-1,:,:,0]   = 1.0
m0[(-1,0),:,:,1] = 1.0

# initialize sstate
minimizer = Minimizer([demag, exchange])
#m = minimizer.minimize(m0, 1e-2, 1e-4)
m = minimizer.minimize_event(m0)
write_vtr(m, "data/sp4_m0")

# perform integration with external field
llg = LLGSolver([demag, exchange, external], material, m)
def E(m):
    return demag.E(material, m) + exchange.E(material, m) + external.E(material, m)

res = llg.solve(1e-9, 1e-12)
t = torch.arange(0,1e-9,1e-12)
with open('data/sp4.dat', 'w') as f:
    for i in range(res.shape[0]):
        m = res[i,:,:,:]
        if i % 10 == 0:
            write_vtr(m, "data/sp4_m_%04d" % (i/10))
        f.write("%g %g %g %g %g\n" % ((t[i],) + tuple(torch.mean(m, axis=(0,1,2))) + (E(m),)))
