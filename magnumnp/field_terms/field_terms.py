from magnumnp.common import constants
import torch

__all__ = ["FieldTerm", "LinearFieldTerm"]

class FieldTerm(object):
    parameters = []

    def __init__(self, **kwargs):
        unknown_params = set(kwargs.keys()) - set(self.parameters)
        if unknown_params:
            raise Warning("Got unknown parameters '%s'. Ignoring!" % unknown_params)
        params = {key:key for key in self.parameters}
        params.update(kwargs)
        for key, value in params.items():
            setattr(self, key, value)

class LinearFieldTerm(FieldTerm):
    def E(self, state):
        return -0.5 * constants.mu_0 * state.mesh.cell_volume * torch.sum(state.material["Ms"] * state.m * self.h(state))
