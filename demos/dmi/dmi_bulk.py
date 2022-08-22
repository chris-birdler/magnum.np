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
state = State(mesh)
state.material = {
        "alpha": 1.,
        }

XX, YY, ZZ = state.SpatialCoordinates()
R = 50e-9
disk = XX**2. + YY**2. < R**2.

Ms = state.Constant([0])
Ms[disk] = 0.86e6 # [A/m]
state.material["Ms"] = Ms

A = state.Constant([0])
A[disk] = 13e-12 # [pJ/m]
state.material["A"] = A

Db = state.Constant([0])
Db[disk] = 3e-3 # [mJ/m²]
state.material["Db"] = Db

K = state.Constant([0])
K[disk] = 0.4e6 # [J/m³]
state.material["Ku"] = K


K_axis = state.Constant([0,0,1.0])
K_axis[~disk] = 0.0
state.material["Ku_axis"] = K_axis

# initialize field terms
exchange = ExchangeField()
aniso = UniaxialAnisotropyField()
dmi = BulkDMIField()

# initialize magnetization that relaxes into s-state
state.m = state.Constant([0.0, 0.0, 0.0]) 
state.m[disk] = state.Tensor((0.,0.,-1.))
state.m[n[0]//2-5:n[0]//2+5, n[1]//2-5:n[1]//2+5, :, 2] = 1.0
state.m.normalize()
write_vti(state.m, "data/m_init.vti", state)

# relax structure without external fields
llg = LLGSolver([exchange, aniso, dmi])
i = 0
with open('data/log.dat', 'w') as f:
    while state.t < 5e-9-eps:
        llg.step(state, 1e-11)
        f.write("%g %g %g %g\n" % ((state.t,) + tuple(state.m.avg())))
        f.flush()

        if i % 10 == 0:
            write_vti(state.m, "data/m_%.5d.vti" % (i//10))
        i += 1

Timer.print_report()
