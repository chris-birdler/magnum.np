import torch

__all__ = ["Material"]

class Material(dict):
    def __init__(self, state):
        self._state = state

    def __getitem__(self, key):
        return super().__getitem__(key)(self._state.t)

    def __setitem__(self, key, value):
        value = self._state.Tensor(value)
        super().__setitem__(key, value)
