from magnumnp import *
import torch

Timer.enable()

# initialize mesh
n  = (40, 40, 1)
dx = (2.5e-9, 2.5e-9, 10e-9)
mesh = Mesh(n, dx)

# initialize material
material = {
        "Ms": 8e5,
        "A": 1.3e-11,
        "alpha": 0.1,
        "xi": 0.05,
        "b": 72.17e-12
        }

state = State(mesh, material)

# initialize magnetization that relaxes into s-state
state.m = state.zeros(n + (3,))
state.m[:20,:,:,1] = -1.
state.m[20:,:,:,1] = 1.
state.m[20,20,:,1] = 0.
state.m[20,20,:,2] = 1.

state.j = state.tensor((1e12, 0, 0))

# initialize field terms
demag    = DemagField(state)
exchange = ExchangeField(state)
torque   = SpinTorqueZhangLi(state)

# initialize sstate
minimizer = Minimizer(state, [demag, exchange])
minimizer.minimize_event()

# perform integration with external field
llg = LLGSolver(state, [demag, exchange, torque])
i = 0
with open('data/m.dat', 'w') as f:
    while state.t < 5e-9:
        if i % 10 == 0:
            write_vti(state.m, "data/m_%04d.vti" % (i/10), state)
        f.write("%g %g %g %g\n" % ((state.t,) + tuple(torch.mean(state.m, axis=(0,1,2)))))
        f.flush()

        m = llg.step(1e-10)
        i += 1

Timer.print_report()
