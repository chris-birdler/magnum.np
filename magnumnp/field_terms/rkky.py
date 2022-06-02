from magnumnp.common import timedmethod, constants
import torch

__all__ = ["RKKYField"]

class RKKYField(object):
    def __init__(self, state, J_rkky, dir, id1, id2):
        self._state = state
        self._mesh = state._mesh

        self._J_rkky = J_rkky
        self._Ms = state._material["Ms"]
        self._dir = dir
        self._id1 = id1
        self._id2 = id2

        self._h = state.zeros(self._mesh.n + (3,))

    @timedmethod
    def h(self, t, m):
        h1 = m2*J_rkky 
        h2 = m1*J_rkky 

    def E(self, t, m):
        E = (m[:,:,(self._id2,),:] * m[:,:,(self._id1,),:]).sum()
        E *= -self._mesh.dx[0]*self._mesh.dx[1]* self._J_rkky
        return E
