from magnumnp import *
import matplotlib.pyplot as plt
import numpy as np

Timer.enable()

n = (400, 20, 1)
dx = (2.5e-9, 2.5e-9, 1e-9)

mesh = Mesh(n, dx)
state = State(mesh)
state.material = {
        "A": 13e-12,
        "Ms": 800e3,
        }

state.m = state.Constant([1.,0.,0.])

demag    = DemagField()
exchange = ExchangeField()
external = ExternalField([804e3,0.0,0.0])

# calculate groundstate
try:
    mesh0, fields0 = read_vti("data/m0_dispersion.vti")
    state.m[...] = fields0["m0"]
except:
    with Timer("Calculate Groundstate"):
        llg = LLGSolver([demag, exchange, external])
        logger = ScalarLogger("data/m0_dispersion.dat", ['t', 'm'])
        while state.t < 5e-9:
            llg.step(state, 1e-11, alpha=1.)
            logger << state
        write_vti({"m0":state.m}, "data/m0_dispersion.vti", state)

# calculate eigenmodes
with Timer("Calculate Eigenvectors"):
    try:
        res = EigenResult.load(state, "data/eigen.pt")
    except:
        eigen = EigenSolver(state, [demag, exchange], [external])
        res = eigen.solve(k=200, tol=1e-6)
        res.store("data/eigen.pt")

with Timer("Store evecs"):
    print("evals[GHz]:", res.freq.cpu().numpy()*1e-9)
    res.save_evecs3D("data/evecs.pvd")

with Timer("Calculate Dispersion"):
    points = lambda vvv: vvv[:,10,0,2,:]
    omega = res.omega.cpu().numpy()
    kk, ww, mm = res.dispersion(points, dx[0])

    fig, ax = plt.subplots()
    cbmax = np.log10(mm*2.).max()
    cbmin = np.log10(mm*2.).min()
    ax.pcolormesh(kk*1e-9, ww*1e-9/2./np.pi, np.log10(np.abs(mm)**2.), cmap = "viridis") #, vmax = cbmax, vmin = cbmin, shading = "auto")

    ref1 = np.loadtxt("ref/disp1.dat")
    ref2 = np.loadtxt("ref/disp2.dat")
    ref3 = np.loadtxt("ref/disp3.dat")
    ref4 = np.loadtxt("ref/disp4.dat")
    ax.plot(ref1[:,0], ref1[:,1], 'k--', linewidth=2)
    ax.plot(ref2[:,0], ref2[:,1], 'k--', linewidth=2)
    ax.plot(ref3[:,0], ref3[:,1], 'k--', linewidth=2)
    ax.plot(ref4[:,0], ref4[:,1], 'k--', linewidth=2)

    ax.set_xlim([-0.3, 0.3])
    ax.set_ylim([40, 80])
    ax.set_xlabel('$k_{x}\,\mathrm{[rad/nm]}$')
    ax.set_ylabel('Frequency f [GHz]')
    fig.savefig("data/dispersion.png")

Timer.print_report()

