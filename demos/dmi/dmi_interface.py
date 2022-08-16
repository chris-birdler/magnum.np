from magnumnp import *
import torch

Timer.enable()

# initialize mesh
eps = 1e-15
n  = (50, 50, 1)
dx = (2e-9, 2e-9, 2e-9)
origin = (-n[0]*dx[0]/2., -n[1]*dx[1]/2., -n[2]*dx[2]/2.)
mesh = Mesh(n, dx, origin=origin)

# initialize material
state = State(mesh, {})
state._material = {
        "alpha": 1.,
        "A":  13e-12, # [J/m]
        "Di": 3.0e-3, # [J/m2]
        "Ku": 0.4e6,   # [J/m3]
        "Ku_axis": [0.0, 0.0, 1.0]
        }

XX, YY, ZZ = state.SpatialCoordinates()
R = 50e-9
disk = XX**2. + YY**2. < R**2.
Ms = state.zeros(n + (1,))
Ms[disk] = 0.86e6 # [A/m]
state._material["Ms"] = Ms

# initialize field terms
exchange = ExchangeField(state)
aniso = UniaxialAnisotropyField(state)
dmi = InterfaceDMIField(state)

# initialize magnetization that relaxes into s-state
state.m = state.zeros(n + (3,))
state.m[disk] = state.tensor((0.,0.,1.))
state.m[n[0]//2-5:n[0]//2+5, n[1]//2-5:n[1]//2+5, :, 2] = -1.0
write_vti(state.m, "data/m_init.vti", state)

# perform integration with external field
llg = LLGSolver(state, [exchange, aniso, dmi])
i = 0
with open('data/m_interface.dat', 'w') as f:
    while state.t < 5e-9-eps:
        llg.step(1e-11)
        f.write("%g %g %g %g\n" % ((state.t,) + tuple(torch.mean(state.m, axis=(0,1,2)))))
        f.flush()

        if i % 10 == 0:
            write_vti(state.m, "data/m_interface_%.5d.vti" % (i//10))
        i += 1

Timer.print_report()
