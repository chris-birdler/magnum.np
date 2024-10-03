from magnumnp import *
from numpy import pi, sqrt, log, logspace
import pathlib

def run_self_induction():
    Timer.enable()
    this_dir = pathlib.Path(__file__).resolve().parent
    set_log_level(100)

    with open(this_dir / "data" / "log.dat", "w") as fd:
        print("#aspect  L_analytic  L_sim1  L_sim2", file=fd)

        for aspect in logspace(-1, 1, num=11):
            dx = (20e-9, 20e-9, 100e-9)
            n  = (int(aspect*500), 500, 1)
            L = (n[0]*dx[0], n[1]*dx[1], n[2]*dx[2])

            origin = (-L[0]/2., -L[1]/2., -L[2]/2.)
            mesh = Mesh(n, dx, origin = origin)
            state = State(mesh)

            j0 = 1e8
            state.j = state.Constant([0.,0.,0.])
            state.j[0,:,:,1] = -j0
            state.j[-1,:,:,1] = j0
            state.j[:,0,:,0] = j0
            state.j[:,-1,:,0] = -j0

            meas = state.Constant(False, dtype=torch.bool)
            meas[1:-1,1:-1,:] = 1

            # sim1 (int B dA)
            oersted = OerstedField()
            h = oersted.h(state)
            psi = constants.mu_0 * h[meas][:,2].sum() * dx[0] * dx[1]
            i = j0 * dx[0] * dx[2]
            L1 = psi / i

            # sim2 (int A*j dV)
            j_scale = 1/i
            state.j *= j_scale
            A = VectorPotential().A(state)
            psi2 = (A*state.j).sum() * dx[0] * dx[1] * dx[2] * constants.mu_0
            L2 = psi2

            # analytic solution
            H = L[0]
            W = L[1]
            d = dx[0]
            L0 = constants.mu_0 / pi * (-2.*(W+H) + 2.*sqrt(H**2.+W**2.) - H*log( (H+sqrt(H**2.+W**2.)) /W) - W*log( (W+sqrt(H**2.+W**2.)) /H) + H*log(2.*H/(d/2.)) + W*log(2.*W/(d/2.)))

            print(aspect, L0, L1.cpu().numpy(), L2.cpu().numpy(), file=fd)

    Timer.print_report()

if __name__=="__main__":
    run_self_induction()
