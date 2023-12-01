from magnumnp import *
import torch
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('N', type=int)
args = parser.parse_args()
print("N:", args.N)


Timer.enable()

# initialize mesh
n  = (args.N, args.N, args.N)
dx = (1e-9, 1e-9, 1e-9)
mesh = Mesh(n, dx)

# initialize material
state = State(mesh, dtype=torch.float32)
state.material = {
        "Ms": 8e5,
        "A": state.Constant([1.3e-11]),
        "gamma": 2.211e5,
        "Ku": 0.4e6,
        "Ku_axis": (0,0,1),
        "Di": state.Constant([3e-3]),
        "Db": state.Constant([3e-3]),
        "DD2d": state.Constant([3e-3])
        }


demag    = DemagField()
exchange = ExchangeField()
aniso    = UniaxialAnisotropyField()
dmii     = InterfaceDMIField()
dmib     = BulkDMIField()
dmiD2d   = D2dDMIField()
external = ExternalField((0,0,1))

state.m = state.Constant((0,0,0))
state.m[1:-1,:,:,0]   = 1.0
state.m[(-1,0),:,:,1] = 1.0

with Timer("first"):
    h1 = demag.h(state)
    h2 = exchange.h(state)
    h3 = aniso.h(state)
    h4 = dmii.h(state)
    h5 = dmib.h(state)
    h6 = dmiD2d.h(state)
    h7 = external.h(state)

with Timer("warmup"):
    for i in range(1000):
        state.m *= 1.01
        h1 = demag.h(state)
        h2 = exchange.h(state)
        h3 = aniso.h(state)
        h4 = dmii.h(state)
        h5 = dmib.h(state)
        h6 = dmiD2d.h(state)
        h7 = external.h(state)

with Timer("measure"):
    for i in range(10000):
        state.m *= 1.01
        h1 = demag.h(state)
        h2 = exchange.h(state)
        h3 = aniso.h(state)
        h4 = dmii.h(state)
        h5 = dmib.h(state)
        h6 = dmiD2d.h(state)
        h7 = external.h(state)

Timer.print_report()
