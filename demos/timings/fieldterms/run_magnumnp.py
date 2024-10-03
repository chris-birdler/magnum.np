from magnumnp import *
import torch
import argparse

torch.set_default_dtype(torch.float32)

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
        "Ms": state.Constant(8e5),
        "A": state.Constant(1.3e-11),
        "Ku": state.Constant(0.4e6),
        "Ku_axis": state.Constant([0,0,1]),
        "Kc_alpha": state.Constant(0.),
        "Kc_beta": state.Constant(0.),
        "Kc_gamma": state.Constant(0.),
        "Kc1": state.Constant(1e3),
        "Kc2": state.Constant(0.),
        "Di": state.Constant(3e-3),
        "Db": state.Constant(3e-3),
        "DD2d": state.Constant(3e-3),
        "P": state.Constant(0.5669),
        "Lambda": state.Constant(2),
        "epsilon_prime": state.Constant(1),
        "mp": state.Constant([0.940, 0.342, 0]),
        "p": state.Constant([0, -1, 0]),
        "d": state.Constant(dx[2]),
        "J": state.Constant(-4e11),
        "je": state.Constant(6.9e10),
        "xi": state.Constant(0.05),
        "b": state.Constant(72.17e-12),
        "eta_damp": state.Constant(-0.1),
        "eta_field": state.Constant(0.3),
        }
state.j = state.Constant([1e12, 0, 0])

demag          = DemagField()
demagPBC       = DemagFieldPBC()
oersted        = OerstedField()
exchange       = ExchangeField()
aniso          = UniaxialAnisotropyField()
caniso         = CubicAnisotropyField()
dmii           = InterfaceDMIField()
dmib           = BulkDMIField()
dmiD2d         = D2dDMIField()
external       = ExternalField(state.Constant([0,0,1]))
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
benchmark_fieldterm(demagPBC, NN = 1000)
benchmark_fieldterm(oersted, NN = 1000)
benchmark_fieldterm(exchange)
benchmark_fieldterm(aniso)
benchmark_fieldterm(caniso)
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
