import numpy as np
import torch
from scipy import interpolate

__all__ = ["TimeInterpolator"]

class TimeInterpolator(object):
    def __init__(self, state, points):
        self._tp = state.Tensor(list(points.keys()))
        self._fp = state._tensor(list(points.values()))
        self._state = state

    def __call__(self, t):
        i = torch.searchsorted(self._tp, t) # upper index
        i = torch.clamp(i, min=1, max=len(self._tp)-1) # extrapolate on bounds
        tp = self._tp
        fp = self._fp
        return fp[i-1] + (t-tp[i-1]) / (tp[i]-tp[i-1]) * (fp[i] - fp[i-1])

    @property
    def final_time(self):
        return self._tp[-1]
