import numpy as np
from scipy import interpolate

__all__ = ["TimeInterpolator"]

class TimeInterpolator(object):
    def __init__(self, points):
        xp = np.array(list(points.keys()))
        fp = np.array(list(points.values()))
        self._f = interpolate.interp1d(xp, fp, axis=0, fill_value="extrapolate")
        self._final_time = xp[-1]

    def __call__(self, t):
        result = self._f(t)
        if result.shape == ():
            return float(result)
        else:
            return result

    @property
    def final_time(self):
        return self._final_time
