import numpy as np
from scipy import interpolate

__all__ = ["TimeInterpolator"]

class TimeInterpolator(object):
    def __init__(self, state, points):
        xp = np.array(list(points.keys()))
        fp = np.array(list(points.values()))
        self._f = interpolate.interp1d(xp, fp, axis=0, fill_value="extrapolate")
        self._final_time = xp[-1]
        self._state = state

    def __call__(self, t):
        result = self._f(t)
        return self._state.Tensor(result)

    @property
    def final_time(self):
        return self._final_time
