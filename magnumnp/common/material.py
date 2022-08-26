import torch

__all__ = ["Material"]

class Material(object):
    def __init__(self, state):
        self._state = state
        self._material = {}

    def __getitem__(self, key):
        return self._material[key](self._state.t)

    def __setitem__(self, key, value):
        self._material[key] = self._state.Tensor(value)
