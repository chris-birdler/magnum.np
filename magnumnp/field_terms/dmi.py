from magnumnp.common import timedmethod, constants
import torch
import numpy as np

__all__ = ["InterfaceDMIField", "BulkDMIField", "D2dDMIField"]

class DMIField(object):
    def __init__(self, dmi_vector, D = "D"):
        self._dmi_vector = dmi_vector
        self._D = D

    @timedmethod
    def h(self, state):
        h = state._zeros(state.mesh.n + (3,))
        D = state.material[self._D]

        full = slice(None, None)
        current = (slice(None, -1), full, full)
        next = (slice(1, None), full, full)

        for dim in range(3):
            v = state.Tensor(self._dmi_vector[dim]).expand(state.m[next].shape)
            D_avg = torch.where(D[next]*D[current] < 0,
                                torch.sqrt(torch.sqrt(-D[next]*D[current])*torch.abs(D[next]+D[current]) / 2.),
                                2.*D[next]*D[current]/(D[next]+D[current]))
            h[current] += D_avg * torch.linalg.cross(v, state.m[next]   ) / (2.*state.mesh.dx[dim])
            h[next]    -= D_avg * torch.linalg.cross(v, state.m[current]) / (2.*state.mesh.dx[dim])

            # rotate dimension
            current = current[-1:] + current[:-1]
            next = next[-1:] + next[:-1]

        h *= 2. / (constants.mu_0 * state.material["Ms"])
        h = torch.nan_to_num(h, posinf=0, neginf=0)
        return state.Tensor(h)

    def E(self, state):
        return -0.5 * constants.mu_0 * state.mesh.cell_volume * torch.sum(state.material["Ms"] * m * self.h(state))


class InterfaceDMIField(DMIField):
    def __init__(self, Di = "Di"):
        dmi_vector = [[ 0, 1, 0], # x
                      [-1, 0, 0], # y
                      [ 0, 0, 0]] # z
        super().__init__(dmi_vector, Di)


class BulkDMIField(DMIField):
    def __init__(self, Db = "Db"):
        dmi_vector = [[1, 0, 0], # x
                      [0, 1, 0], # y
                      [0, 0, 1]] # z
        super().__init__(dmi_vector, Db)


class D2dDMIField(DMIField):
    def __init__(self, DD2d = "DD2d"):
        dmi_vector = [[-1, 0, 0], # x
                      [ 0, 1, 0], # y
                      [ 0, 0, 0]] # z
        super().__init__(dmi_vector, DD2d)
