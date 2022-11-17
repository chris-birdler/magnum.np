from magnumnp.common import constants
import torch

__all__ = ["FieldTerm", "LinearFieldTerm"]

class FieldTerm(object):
    parameters = []

    def __init__(self):
        for key in parameters:
            self.setattr(key, key)

class LinearFieldTerm(FieldTerm):
    def E(self, state):
        return -0.5 * constants.mu_0 * state.mesh.cell_volume * torch.sum(state.material["Ms"] * state.m * self.h(state))
