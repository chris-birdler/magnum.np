import pytest
import types
import torch
from magnumnp import RKF45


class OldRKF45(RKF45):
    """step() as before the change: it can try a step past the end of the interval, which is then rejected."""
    def step(self, t, x, dt, rtol = None, atol = None, **kwargs):
        remaining = float(dt)
        while remaining > 0.:
            _x1, _t1, err = self._try_step(t, x, **kwargs)
            if atol == None:
                atol = self._atol
            dt_opt = self._optimal_stepsize(x, err, atol, rtol or self._rtol)
            if self._dt > dt_opt or self._dt > remaining:
                self._dt = min(dt_opt, remaining)
            else:
                remaining -= self._dt
                x = _x1
                t = _t1
                self._dt = dt_opt
                kwargs["state"]._step += 1
        return t, x


def run(solver_cls, n_intervals, interval):
    calls = []
    def f(t, x, **kwargs):
        calls.append((float(t), solver._dt))
        return -3e10 * x + 1e10 * torch.sin(7e10 * t) * torch.ones_like(x)
    solver = solver_cls(f, dt = 1e-13, atol = 1e-6, rtol = 1e-6)
    state = types.SimpleNamespace(_dt = None, _step = 0)
    t = torch.tensor(0., dtype = torch.float64)
    x = torch.linspace(0.5, 1.5, 7, dtype = torch.float64)
    ends = []
    for i in range(n_intervals):
        t, x = solver.step(t, x, interval, state = state)
        ends.append(float(t))
    return x, calls, ends


def test_no_try_past_interval_end():
    interval = 1e-10
    x, calls, ends = run(RKF45, 20, interval)
    # each try (6 stages) starts at the time of its first stage call; it must not pass the end of its interval
    for k in range(0, len(calls), 6):
        t0, dt = calls[k]
        end = (int(t0 / interval + 1e-9) + 1) * interval
        assert t0 + dt <= end * (1 + 1e-12)


def test_same_result_fewer_tries():
    # smooth problem: the step past the end would be accepted by the error control, so the old code retries with
    # exactly the remaining time; the new code tries it at once and gives the same result bitwise
    x_new, calls_new, ends_new = run(RKF45, 20, 1e-10)
    x_old, calls_old, ends_old = run(OldRKF45, 20, 1e-10)
    assert torch.equal(x_new, x_old)
    assert ends_new == ends_old
    assert len(calls_new) < len(calls_old)


class OldTryRKF45(RKF45):
    """_try_step and the error norm as before the compiled combinations (audit item A3)."""
    def _try_step(self, t, x, **kwargs):
        f, dt = self._f, self._dt
        k1 = dt * f(t,              x, **kwargs)
        k2 = dt * f(t +  1./ 4.*dt, x +      1./ 4.*k1, **kwargs)
        k3 = dt * f(t +  3./ 8.*dt, x +      3./32.*k1 +      9./32.*k2, **kwargs)
        k4 = dt * f(t + 12./13.*dt, x + 1932./2197.*k1 - 7200./2197.*k2 + 7296./2197.*k3, **kwargs)
        k5 = dt * f(t +      1.*dt, x +   439./216.*k1 -          8.*k2 + 3680./ 513.*k3 -  845./4104.*k4, **kwargs)
        k6 = dt * f(t +  1./ 2.*dt, x -    8. / 27.*k1 +          2.*k2 - 3544./2565.*k3 + 1859./4104.*k4 - 11./40.*k5, **kwargs)
        dx = 16./135.*k1 + 6656./12825.*k3 + 28561./56430.*k4 - 9./50.*k5 + 2./55.*k6
        rk_error = dx - (25./216.*k1 + 1408./2565.*k3 + 2197./4104.*k4 - 1./5.*k5)
        return x+dx, t+dt, rk_error

    def _optimal_stepsize(self, x, err, atol, rtol):
        self.norms.append(float(torch.linalg.norm((err / (atol + rtol*x.abs())).flatten(), torch.inf)))
        return super()._optimal_stepsize(x, err, atol, rtol)


class NewTryRKF45(RKF45):
    def _optimal_stepsize(self, x, err, atol, rtol):
        self.norms.append(float(torch.linalg.norm((err / (atol + rtol*x.abs())).flatten(), torch.inf)))
        return super()._optimal_stepsize(x, err, atol, rtol)


@pytest.mark.parametrize("dtype, tol", [(torch.float64, 1e-12), (torch.float32, 1e-5)])
def test_compiled_try_against_reference(dtype, tol):
    """One try and 1000 steps (PLAN 11, S1): the compiled combinations give x1 and the error within round-off, and
    in double precision the same sequence of accepted and rejected tries."""
    g = torch.Generator().manual_seed(3)
    B = torch.randn(3, 3, generator=g, dtype=torch.float64)
    A = ((B - B.T) * 1e10 - 1e9 * torch.eye(3, dtype=torch.float64)).to(dtype)   # rotation and damping
    x0 = torch.randn(4, 4, 4, 3, generator=g, dtype=torch.float64).to(dtype)
    f = lambda t, x, **kw: x @ A.T + 1e10 * torch.sin(7e10 * t)
    state = types.SimpleNamespace(_dt = None, _step = 0)
    t0 = torch.tensor(0., dtype = dtype)
    new, old = NewTryRKF45(f, dt = 1e-12, atol = 1e-6, rtol = 1e-6), OldTryRKF45(f, dt = 1e-12, atol = 1e-6, rtol = 1e-6)
    (x1n, _, en), (x1o, _, eo) = new._try_step(t0, x0, state = state), old._try_step(t0, x0, state = state)
    assert float((x1n - x1o).abs().max() / x1o.abs().max()) < tol
    assert float((en - eo).abs().max() / eo.abs().max()) < tol * 1e3   # the error is a small difference of terms
    if dtype == torch.float64:
        new.norms, old.norms = [], []
        tn, xn, to, xo = t0, x0, t0, x0
        for i in range(500):
            tn, xn = new.step(tn, xn, 2e-11, state = state)
            to, xo = old.step(to, xo, 2e-11, state = state)
        accept = lambda norms: [v <= 1.1 for v in norms]
        assert len(new.norms) == len(old.norms) and accept(new.norms) == accept(old.norms)
        assert len(new.norms) >= 1000
        assert float((xn - xo).abs().max() / xo.abs().max()) < 1e-9
