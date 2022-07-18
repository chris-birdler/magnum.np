from magnumnp.common import timedmethod, constants
import torch

__all__ = ["RKKYField"]

class RKKYField(object):
    r"""
    Effective field contribution corresponding to the RKKY interaction between two interface layers.

    .. math::

        E_\text{RKKY} = \int_\Gamma J_\text{RKKY} \, \vec{m}_1 \codt \vec{m}_2 d\vec{A}

    :param state: current state
    :type state: State

    :Example:

      .. code::

        # create state with named domains from mesh
        state = State(mesh)

        # create domains as bool arrays, e.g:
        domain1 = state.zeros(n, dtype=torch.bool)
        domain1[n[0]//2:,:,:] = True

        domain2 = state.zeros(n, dtype=torch.bool)
        domain2[:-n[0]//2:,:,:] = True

        # rotate magnetization within one subdomain
        state.m[domain1] = state.tensor([np.cos(phi), np.sin(phi), 0])

        # without interface layer, two seperate exchange fields need to be defined
        exchange1 = ExchangeField(state, Aex1, domain1)
        exchange2 = ExchangeField(state, Aex2, domain2)
        rkky = RKKYField(state, domain1, domain2)
    """
    def __init__(self, state, J_rkky, dir, id1, id2, order = 0):
        self._state = state
        self._mesh = state._mesh

        self._J_rkky = J_rkky
        self._Ms = state._material["Ms"]
        self._dir = dir #TODO: dir is ignored
        self._id1 = min(id1,id2)
        self._id2 = max(id1,id2)

        self._order = order
        self._h = state.zeros(self._mesh.n + (3,))

    @timedmethod
    def h(self, t, m):
        m1 = m[:,:,(self._id1,),:]
        m2 = m[:,:,(self._id2,),:]
        if self._order == 1:
            m1 += 0.5 * (m[:,:,(self._id1,),:] - m[:,:,(self._id1-1,),:])
            m2 += 0.5 * (m[:,:,(self._id2,),:] - m[:,:,(self._id2+1,),:])
        elif self._order == 2:
            m1 += 0.25 * (3*m[:,:,(self._id1,),:] - 4*m[:,:,(self._id1-1,),:] + m[:,:,(self._id1-2,),:])
            m2 += 0.25 * (3*m[:,:,(self._id2,),:] - 4*m[:,:,(self._id2+1,),:] + m[:,:,(self._id2+2,),:])

        self._h[:,:,(self._id1,),:] = self._J_rkky * (m2 - (m1*m2).sum(axis = 3) * m1)
        self._h[:,:,(self._id2,),:] = self._J_rkky * (m1 - (m1*m2).sum(axis = 3) * m2)

        self._h /= constants.mu_0 * self._Ms * self._mesh.dx[2]
        return torch.nan_to_num(self._h)

    def E(self, t, m):
        m1 = m[:,:,(self._id1,),:]
        m2 = m[:,:,(self._id2,),:]
        if self._order == 1:
            m1 += 0.5 * ( m[:,:,(self._id1,),:] - m[:,:,(self._id1-1,),:])
            m2 += 0.5 * (-m[:,:,(self._id2,),:] + m[:,:,(self._id2+1,),:])
        elif self._order == 2:
            m1 += 0.25 * (3*m[:,:,(self._id1,),:] - 4*m[:,:,(self._id1-1,),:] + m[:,:,(self._id1-2,),:])
            m2 += 0.25 * (3*m[:,:,(self._id2,),:] - 4*m[:,:,(self._id2+1,),:] + m[:,:,(self._id2+2,),:])

        E = (m2*m1).sum()
        E *= -self._mesh.dx[0] * self._mesh.dx[1] * self._J_rkky
        return E


