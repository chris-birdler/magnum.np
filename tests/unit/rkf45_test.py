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
