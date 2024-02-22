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
state = State(mesh)
state.material = {
        "Ms": 8e5,
        "A": state.Constant([1.3e-11]),
        "gamma": 2.211e5,
        "Ku": 0.4e6,
        "Ku_axis": (0,0,1),
        "Di": state.Constant([3e-3]),
        "Db": state.Constant([3e-3]),
        "DD2d": state.Constant([3e-3]),
        "P": 0.5669,
        "Lambda": 2,
        "epsilon_prime": 1,
        "mp": state.Tensor((0.940, 0.342, 0)),
        "p": [0, -1, 0],
        "d": dx[2],
        "J": -4e11,
        "je": 6.9e10,
        "xi": 0.05,
        "b": 72.17e-12,
        "eta_damp": -0.1,
        "eta_field": 0.3,
        }
state.j = state.Tensor((1e12, 0, 0))


demag          = DemagField()
exchange       = ExchangeField()
aniso          = UniaxialAnisotropyField()
dmii           = InterfaceDMIField()
dmib           = BulkDMIField()
dmiD2d         = D2dDMIField()
external       = ExternalField((0,0,1))
st_slonczewski = SpinTorqueSlonczewski()
st_sot         = SpinOrbitTorque()
st_zhangli     = SpinTorqueZhangLi()


def benchmark_fieldterm(term, NN = 10000):
    state.m = state.Constant((0,0,0))
    state.m[1:-1,:,:,0]   = 1.0
    state.m[(-1,0),:,:,1] = 1.0

    with Timer("setup"):
        h = term.h(state)

    with Timer("warmup"):
        for i in range(NN // 10):
            state.m *= 1.01
            h = term.h(state)

    with Timer("measure"):
        for i in range(NN):
            state.m *= 1.01
            h = term.h(state)

benchmark_fieldterm(demag, NN = 1000)
benchmark_fieldterm(exchange)
benchmark_fieldterm(aniso)
benchmark_fieldterm(dmii)
benchmark_fieldterm(dmib)
benchmark_fieldterm(dmiD2d)
benchmark_fieldterm(external)
benchmark_fieldterm(st_slonczewski)
benchmark_fieldterm(st_sot)
benchmark_fieldterm(st_zhangli)

Timer.print_report()
names = [Timer._timers[t]["name"] for t in Timer._timers.keys() if t.startswith("###measure###")]
t_setup = [Timer._timers[t]["total_time"] / Timer._timers[t]["calls"] for t in Timer._timers.keys() if t.startswith("###setup###")]
t_eval = [Timer._timers[t]["total_time"] / Timer._timers[t]["calls"] for t in Timer._timers.keys() if t.startswith("###measure###")]

print("#names:   N= ", args.N, *names, sep='\t')
print("t_setup: N= ", args.N, *t_setup, sep='\t')
print("t_eval:  N= ", args.N, *t_eval, sep='\t')
