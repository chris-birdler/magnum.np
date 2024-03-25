from magnumnp import *
import torch
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('e', type=int)
args = parser.parse_args()

N = 2**args.e
print("N:", N)

Timer.enable()

# initialize mesh
n  = (N, N, 1)
dx = (4e-9, 4e-9, 4e-9)
mesh = Mesh(n, dx)

# initialize material
state = State(mesh, dtype=torch.float32)
state.material = {
        "Ms": 800e3,
        "A": 13e-12,
        "alpha": 0.01,
        }

demag    = DemagField()
exchange = ExchangeField()
external = ExternalField([0,0.01/constants.mu_0,0])

state.m = state.Constant([1,0,0])

llg = LLGSolver([demag, exchange, external])

with Timer("warmup"):
    llg.step(state, 1e-12)
with Timer("measure"):
    llg.step(state, 1e-11)

Timer.print_report()

t_total = Timer._timers['###measure']['total_time']
n_eval = Timer._timers['###measure###LLGSolver.step###DemagField.h']['calls']

print("NN=", N*N, "n_eval=", n_eval, "t_total=", t_total, "throuhput=", N*N*n_eval/t_total)
