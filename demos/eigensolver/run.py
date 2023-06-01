from magnumnp import *
from torch import sqrt



n = (20, 20, 1)
dx = (10e-9, 10e-9, 1e-9)
origin = (-1./2.*n[0]*dx[0], -1./2.*n[1]*dx[1], -1./2.*n[2]*dx[2])
mesh = Mesh(n, dx, origin)
state = State(mesh)
x,y,z = state.SpatialCoordinate()
r = x**2 + y**2
state.m = state.Tensor(torch.stack([y/r, -x/r, 1e-9/r], dim = -1))
state.m.normalize()
write_vti(state.m, "m0.vti")

state.material = {'Ms': 800e3, 'A': 13e-12, 'Ku': -5e3, 'Ku_axis': (0,0,1)}

class EigenSolverTest(object):
    def __init__(self, state, linear_terms, constant_terms):
        #self._linear_terms = linear_terms
        #self._constant_terms = constant_terms
        #h = sum([term.h(state) for term in self._terms])
 
        self._state = state
 
        ez = state.Constant([1e-15,0.,1.])
        self._e1 = torch.linalg.cross(ez, self._state.m)
        self._e1 = self._e1 / torch.linalg.norm(self._e1, axis=3).unsqueeze(dim = -1)
        self._e0 = -torch.linalg.cross(self._e1, self._state.m)
        self._e0 = self._e0 / torch.linalg.norm(self._e0, axis=3).unsqueeze(dim = -1)
 
        print(self._e0.shape, self._e1.shape)
        write_vti(self._e0, "e0.vti")

eig = EigenSolverTest(state, [], [])
exchange = ExchangeField()
demag = DemagField()
external = ExternalField([10e3, 0, 0])
eigs = EigenSolver(state, [demag, exchange], [external])
res = eigs.solve(k = 5)


print(res[0])
print(res[1].shape)
