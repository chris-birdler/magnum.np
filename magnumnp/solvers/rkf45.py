import torch

__all__ = ["RKF45"]

#Runge-Kutta-Fehlberg method with stepsize control
class RKF45(object):
    def __init__(self, f, y, t, dt=1e-15):
        self._f = f
        self._y = y
        self._t = t
        self._dt = torch.DoubleTensor([dt])
        self._order = 4

        # Numerical Recipies 3rd Edition suggests these values:
        self._headroom = 0.9
        self._maxstep = torch.Tensor([1e-11])
        self._minscale = 0.2
        self._maxscale = 10.
        self._atol = 1e-5

    def _try_step(self):
        t, dt, y, f = self._t, self._dt, self._y, self._f
        k1 = dt * f(t,              y)
        k2 = dt * f(t +  1./ 4.*dt, y +      1./ 4.*k1)
        k3 = dt * f(t +  3./ 8.*dt, y +      3./32.*k1 +      9./32.*k2)
        k4 = dt * f(t + 12./13.*dt, y + 1932./2197.*k1 - 7200./2197.*k2 + 7296./2197.*k3)
        k5 = dt * f(t +      1.*dt, y +   439./216.*k1 -          8.*k2 + 3680./ 513.*k3 -  845./4104.*k4)
        k6 = dt * f(t +  1./ 2.*dt, y -    8. / 27.*k1 +          2.*k2 - 3544./2565.*k3 + 1859./4104.*k4 - 11./40.*k5)

        dy = 16./135.*k1 + 6656./12825.*k3 + 28561./56430.*k4 - 9./50.*k5 + 2./55.*k6
        rk_error = dy - (25./216.*k1 + 1408./2565.*k3 + 2197./4104.*k4 - 1./5.*k5)
        return y+dy, t+dt, rk_error

    def _optimal_stepsize(self, rk_error):
        norm = torch.linalg.norm(rk_error.flatten() / self._atol, torch.inf)

        if norm > 1.1:
            # decrease step, no more than factor of 5, but a fraction S more
            # than scaling suggests (for better accuracy)
            r = self._headroom / torch.pow(norm, 1.0/self._order)
            if (r < self._minscale):
                r = self._minscale
        elif norm < 0.5:
            # increase step, but no more than by a factor of 5
            r = self._headroom / torch.pow(norm, 1.0/(self._order+1.0));
            if r > self._maxscale: # increase no more than factor of 5
                r = self._maxscale
            if r < 1.: # don't allow any decrease caused by S<1
                r = 1.
        else: # no change
            r = 1.

        dt_opt = self._dt * r
        if dt_opt > self._maxstep:
            dt_opt = self._maxstep
        return dt_opt

    def step(self, dt):
        t0, t1 = self._t, self._t + dt
        while self._t < t1:
            _y1, _t1, err = self._try_step()
            dt_opt = self._optimal_stepsize(err)
            if self._dt > dt_opt or self._dt > t1 - self._t:
                # step size was too large, retry with optimal stepsize
                self._dt = min(dt_opt, t1 - self._t)
                print("REVERT step: %g, new step size: %g, t1-_t: %g" % (self._dt, dt_opt, t1 - self._t))
            else:
                # accept step, adapt stepsize for next step
                self._y = _y1
                self._t = _t1
                self._dt = dt_opt
                print("ACCEPT step: %g, new step size: %g, time: %g" % (self._dt, dt_opt, self._t))
