from magnumnp import *
from torch import sqrt
import matplotlib.pyplot as plt
import numpy as np

Timer.enable()

lex = 5.71e-9
Js = 1.
n  = (50, 50, 1)
dx = (5e-9, 5e-9, 20e-9)

mesh = Mesh(n, dx, origin=(-n[0]*dx[0]/2.,-n[1]*dx[1]/2.,-n[2]*dx[2]/2.))
state = State(mesh)
state.material = {
        "A": 13e-12,
        "Ms":1./constants.mu_0,
        }

x, y, z = mesh.SpatialCoordinate()
r = n[0]*dx[0]/2.
disk = x**2 + y**2 < r**2
state.m = Expression([-y,x,0*z])
normalize(state.m)

state.material["A"][~disk] = 0.
state.material["Ms"][~disk] = 0.
write_vti(state.material, "data/material.vti", state)

demag    = DemagField()
exchange = ExchangeField()

try:
    mesh0, fields0 = read_vti("data/m0_vortex_cylindric.vti")
    state.m[...] = fields0["m0"]
except:
    with Timer("Calculate Groundstate"):
        llg = LLGSolver([demag, exchange])
        logger = ScalarLogger("data/m0_vortex_cylindric.dat", ['t', 'm'])
        while state.t < 5e-9:
            llg.step(state, 1e-11, alpha=1.)
            logger << state
        write_vti({"m0":state.m}, "data/m0_vortex_cylindric.vti", state)

with Timer("Calculate Eigenvectors"):
    try:
        res = EigenResult.load(state, "data/eigen_vortex_cylindric.pt")
    except:
        eigen = EigenSolver(state, [demag, exchange], [], domain = disk)
        res = eigen.solve(k=20, tol=1e-6)
        res.store("data/eigen_vortex_cylindric.pt")

with Timer("Store evecs"):
    print("evals[GHz]:", res.freq.numpy()*1e-9)
    res.save_evecs3D("data/evecs.pvd")

with Timer("Plot Results"):
    ref_magnumpi = np.loadtxt("ref/vortex_cylindric_magnumpi.dat")
    ref_daquino = np.loadtxt("ref/vortex_cylindric_daquino.dat")

    fig, ax = plt.subplots()
    cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    ax.plot(ref_daquino[1], '^', fillstyle = "none", color = "black", alpha = 0.7, linewidth=0.8, label = "d'Aquino (FD)")
    ax.plot(ref_magnumpi[2,1:].T, 'o', mfc = 'none', color = cycle[1], label = "magnum.pi (FEM)")
    ax.plot(res.freq.numpy()*1e-9, '+-', color = cycle[0], label = "magnum.np (FD)")

    ax.set_xlim([-0.1,5.1])
    #ax.set_ylim([7.5,15])
    ax.set_title("Vortex")
    ax.set_xlabel("Eigenmode Index [1]")
    ax.set_ylabel("Frequency f [GHz]")
    ax.legend(loc=2)
    ax.grid()
    fig.savefig("data/results_vortex_cylindric.png")

Timer.print_report()
