
#Runge-Kutta-Fehlberg method with stepsize control
class RKF45(object)
    def __init__(self, f, t, x, dt0=1.):
        self.f = f
        self.t = t
        self.x = x
        self.dt = dt0
        self.order = 4 # 

        # Numerical Recipies 3rd Edition suggests these values:
        self.beta = 0.4 / 5.0
        self.alpha = 0.2 - 0.75 * beta
        self.headroom = 0.9
        self.minscale = 0.2
        self.maxscale = 10.

    def _try_step(self):
        dt, x, f = self.dt, self.x, self.f
        k1 = dt * f(t,              x)
        k2 = dt * f(t +  1./ 4.*dt, x +      1./ 4.*k1)
        k3 = dt * f(t +  3./ 8.*dt, x +      3./32.*k1 +      9./32.*k2)
        k4 = dt * f(t + 12./13.*dt, x + 1932./2197.*k1 - 7200./2197.*k2 + 7296./2197.*k3)
        k5 = dt * f(t +      1.*dt, x +   439./216.*k1 -          8.*k2 + 3680./ 513.*k3 -  845./4104.*k4)
        k6 = dt * f(t +  1./ 2.*dt, x -    8. / 27.*k1 +          2.*k2 - 3544./2565.*k3 + 1859./4104.*k4 - 11./40.*k5)
        sumbk = 16./135.*k1 + 6656./12825.*k3 + 28561./56430.*k4 - 9./50.*k5 + 2./55.*k6
        rk_error = sumbk - (25./216.*k1 + 1408./2565.*k3 + 2197./4104.*k4 - 1./5.*k5)
        scale = self.atol + self.rtol * torch.abs(torch.max(self.x, self.x + sumbk)) # scale from magnumaf
        return (rk_error/scale).max()

    def _accept(self, dt):

    def step(self, dt):
        # shortcuts
        y0 = y.detach()

        while True:
            # Try a step from state.t to state.t+h_try
            dydt = self._try_step(state)
#
#            # Step size control (was h too big?)
#            # calculate the minimal acceptable step size
            accept, h_new = adjust_stepsize(self.tab.order, y, y_err, dydt)
#            if accept:
#                # done -> exit loop
#                break
#            else:
#                # oh, tried step size was too large.
#                y.assign(y0)  # reverse last step
#                h_try = h_new  # try again with new (smaller) h.
#                continue  # need to retry -> redo loop
#
#        # But: Don't overshoot past t_max!
#        if state.t + h_try > t_max:
#            h_try = t_max - state.t   # make h_try smaller.
#            y.assign(y0)              # reverse last step
#            self.apply(state, h_try)  # assume that a smaller step size is o.k.
#
#        # Update state
#        state.t += h_try
#        state.h = h_try
#        state.__runge_kutta_next_h = h_new
#        state.step += 1
#        state.flush_cache()
#        state.finish_step()
#        return state



#    // Numerical Recipies 3rd Edition suggests these values:
#    const double beta_ = 0.4 / 5.0;
#    const double alpha_ = 0.2 - 0.75 * beta_;
#    const double headroom_ = 0.9;
#    const double minscale_ = 0.2;
#    const double maxscale_ = 10.;

