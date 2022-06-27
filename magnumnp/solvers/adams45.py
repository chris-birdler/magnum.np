import torch
from magnumnp.common import logging

__all__ = ["Adams45"]

#Adams-Bashforth method with stepsize control
class Adams45(object):
    def __init__(self, f, y, t, dt=1e-15):
        self._f = f
        self._f1 = f(t,y)
        self._f2 = self._f1 # TODO: maybe use RKF45 for better initialization
        self._f3 = self._f1

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
        dt, t, y, f = self._dt, self._t, self._y, self._f
        f0 = f(t, y)
        f1, f2, f3 = self._f1, self._f2, self._f3
        
        dy3 = dt/12. * (23.*f0 - 16.*f1 +  5.*f2)  
        dy4 = dt/24. * (55.*f0 - 59.*f1 + 37.*f2 - 9.*f3)
        rk_error = dy4-dy3
        return y+dy4, t+dt, rk_error, f0

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
        t0, t_new = self._t, self._t + dt
        while self._t < t_new:
            _y_new, _t_new,  err, f0 = self._try_step()
            dt_opt = self._optimal_stepsize(err)
            if self._dt > dt_opt or self._dt > t_new - self._t:
                # step size was too large, retry with optimal stepsize
                self._dt = min(dt_opt, t_new - self._t)
                logging.debug("REVERT step: %g, new step size: %g, time: %g" % (self._dt, dt_opt, self._t))
            else:
                # accept step, adapt stepsize for next step
                self._y = _y_new
                self._t = _t_new
                self._dt = dt_opt
                self._f1, self._f2, self.f3 = f0, self._f1, self._f2
                logging.debug("ACCEPT step: %g, new step size: %g, time: %g" % (self._dt, dt_opt, self._t))
