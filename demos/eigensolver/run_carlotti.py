from magnumnp import *
import matplotlib.pyplot as plt
import numpy as np

Timer.enable()


n = (200, 100, 1)
dx = (1e-9, 1e-9, 5e-9)

mesh = Mesh(n, dx, origin = (-n[0]*dx[0]/2., -n[1]*dx[1]/2., 0.0))
state = State(mesh, scale = 1e9)

Ms = 800e3
state.material = {
        "A": 13e-12,
        "Ms": Ms,
        "alpha": 0.01
        }

x, y, z = mesh.SpatialCoordinate()
a = 200e-9 / 2.
b = 100e-9 / 2.
magnetic = (x/a)**2. + (y/b)**2. <= 1.

state.m = state.Constant([1.,0.,0.])
#state.m[~magnetic] = 0.
state.material["Ms"][~magnetic] = 0.
state.material["A"][~magnetic] = 0.

demag    = DemagField()
exchange = ExchangeField()

# calculate groundstate
try:
    mesh0, fields0 = read_vti("data/m0_carlotti.vti")
    state.m[...] = fields0["m0"]
except:
    with Timer("Calculate Groundstate"):
        minimizer = MinimizerBB([demag, exchange])
        minimizer.minimize(state)
        state.write_vtk({"m0":state.m}, "data/m0_carlotti.vti")

# calculate eigenmodes
with Timer("Calculate Eigenvectors"):
    try:
        res = EigenResult.load(state, "data/eigen_carlotti.pt")
    except:
        eigen = EigenSolver(state, [demag, exchange], [], domain = magnetic)
        res = eigen.solve(k=20, tol=1e-6)
        res.store("data/eigen_carlotti.pt")

with Timer("Store evecs"):
    print("evals[GHz]:", res.freq.cpu().numpy()*1e-9)
    res.save_evecs3D("data/evecs_carlotti.pvd")

with Timer("Plot Absorbtion"):
    ref = np.loadtxt("ref/carlotti.dat")
    freq = np.arange(0.05e9, 30e9, 0.05e9)
    h_excite = state.Constant([0.,0.5e-3/constants.mu_0,0.])
    absorption = res.absorption(2*np.pi*freq, h_excite)#, magnetic)

    V = magnetic.sum() * mesh.cell_volumes
    P0 = constants.mu_0 * Ms**2 * V * constants.gamma * Ms

    fig, ax = plt.subplots(figsize=(8,8))
    ax.plot(freq * 1e-9, absorption / P0, color="red", linewidth=2.0)
    #ax2 = ax.twinx()
    #ax2.set_yscale("log")
    #ax2.plot(ref[:,0], ref[:,1], "k--")

    ax.plot(ref[:,0], ref[:,1], "k--")
    ax.set_yscale("log")
    ax.set_xlabel("Frequency [GHz]")
    ax.set_ylabel("PSD [$P_0$]")
    ax.grid()
    fig.savefig("result_carlotti.png")

Timer.print_report()
